"""Synthetic local functional tests; no database or clinical fixture is shipped."""
import contextlib
import hashlib
import importlib.util
import io
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
from unittest.mock import patch
from types import SimpleNamespace

ROOT = Path(__file__).resolve().parents[1]
os.environ['PYTHONDONTWRITEBYTECODE'] = '1'
sys.path.insert(0, str(ROOT))
import mycowes
from workflow.reference_profiles import MODES, bindings, validate


class Functionality(unittest.TestCase):
    @staticmethod
    def mock_complete_engine(command, env):
        """Transaction fixture only; this does not run or validate native analysis."""
        attempt = Path(command[command.index('--run-root') + 1])
        (attempt/'run_manifest.tsv').write_text('status\tCOMPLETE\n')
        for mode in ('PRIMARY', 'STRICT_DF2'):
            summary = attempt/'07_breadth'/mode/'07_species_dedup_dispersion_summary.tsv'
            summary.parent.mkdir(parents=True)
            summary.write_text('species\tdispersion_tier\nunit_fixture\tweak\n')
        return subprocess.CompletedProcess(command, 0)

    def make_assets(self, root, mode=MODES[0]):
        values = bindings(root, mode)
        for key, value in values.items():
            if key.endswith(('_FASTA', '_MMI', '_TAXONOMY', '_IDENTITY', '_MANIFEST', '_BED', '_DF2')):
                p = Path(value); p.parent.mkdir(parents=True, exist_ok=True)
                p.write_text('synthetic-test-only\n')
        Path(values['MYCOWES_MICROBIAL_IDENTITY']).write_text('profile_id\t' + values['MYCOWES_MICROBIAL_PROFILE_ID'] + '\n')
        return values

    def test_default_is_one_representative(self):
        self.assertEqual(bindings('x'), bindings('x', MODES[0]))
        self.assertIn('ONE_REP_TARGET', bindings('x')['MYCOWES_JOINT_FASTA'])

    def test_profile_switch_only_changes_microbe_assets(self):
        one, multi = bindings('x'), bindings('x', MODES[1])
        self.assertNotEqual(one['MYCOWES_JOINT_FASTA'], multi['MYCOWES_JOINT_FASTA'])
        self.assertNotEqual(one['MYCOWES_JOINT_MMI'], multi['MYCOWES_JOINT_MMI'])
        for key in ('MYCOWES_TARGET_FASTA', 'MYCOWES_RISK_BED', 'MYCOWES_GLOBAL_CAPTURE_DF2'):
            self.assertEqual(one[key], multi[key])
        self.assertFalse(any('HOST' in k or k.startswith(('CORE_', 'ANCHOR_', 'STRONG_', 'SUPPORTIVE_', 'ADAPTER_')) for k in one))

    def test_both_profile_loads(self):
        with tempfile.TemporaryDirectory() as t:
            for mode in MODES:
                validate(self.make_assets(t, mode))

    def test_missing_reference_does_not_fall_back(self):
        with tempfile.TemporaryDirectory() as t:
            self.make_assets(t, MODES[1])
            with self.assertRaisesRegex(ValueError, 'Missing'):
                validate(bindings(t))

    def test_empty_asset_rejected(self):
        with tempfile.TemporaryDirectory() as t:
            v = self.make_assets(t);Path(v['MYCOWES_JOINT_MMI']).write_text('')
            with self.assertRaisesRegex(ValueError, 'empty'):validate(v)

    def test_mismatched_identity_rejected(self):
        with tempfile.TemporaryDirectory() as t:
            v = self.make_assets(t);Path(v['MYCOWES_MICROBIAL_IDENTITY']).write_text('profile_id\twrong\n')
            with self.assertRaisesRegex(ValueError, 'identity mismatch'):validate(v)

    def test_unknown_profile_rejected(self):
        with self.assertRaises(ValueError):bindings('x', 'typo')

    def test_quality_profiles_and_boundary_logic(self):
        spec=importlib.util.spec_from_file_location('risk',ROOT/'workflow/annotate_fragment_risk_v0.1.py')
        risk=importlib.util.module_from_spec(spec);spec.loader.exec_module(risk)
        row={'pair_support_class':'paired_target_consistent','R1_primary_aligned_bases':'60','R2_primary_aligned_bases':'60',
             'R1_primary_MAPQ':'5','R2_primary_MAPQ':'5','maximum_primary_NM_rate':'0.12'}
        self.assertEqual(risk.core_filter_reason(row,*mycowes.DIFFICULTIES['simple']),'PASS')
        self.assertNotEqual(risk.core_filter_reason(row,*mycowes.DIFFICULTIES['normal']),'PASS')
        row.update(R1_primary_aligned_bases='100',R2_primary_aligned_bases='100',R1_primary_MAPQ='20',R2_primary_MAPQ='20',maximum_primary_NM_rate='0.08')
        self.assertEqual(risk.core_filter_reason(row,100),'PASS')
        self.assertEqual(risk.core_filter_reason(row,*mycowes.DIFFICULTIES['normal']),'PASS')
        row['maximum_primary_NM_rate']='0.0801'
        self.assertNotEqual(risk.core_filter_reason(row,*mycowes.DIFFICULTIES['normal']),'PASS')

    def test_explicit_preset_resolution_and_conflicts(self):
        base=dict(preset=None,difficulty=None,reference_mode=None,min_aligned_bp=None,min_mapq=None,max_nm_rate=None)
        resolve=lambda **kw:mycowes.resolve_run_settings(SimpleNamespace(**dict(base,**kw)))
        self.assertEqual(resolve(),('standard',MODES[0],(100,20,.08)))
        self.assertEqual(resolve(preset='relax2'),('relax2',MODES[0],(60,5,.12)))
        self.assertEqual(resolve(difficulty='hard'),('standard',MODES[1],(100,20,.08)))
        self.assertEqual(resolve(preset='relax2',difficulty='simple')[2],(60,5,.12))
        for difficulty,preset in [('simple','standard'),('normal','relax2'),('hard','relax2')]:
            with self.assertRaisesRegex(ValueError,'conflicts'):resolve(preset=preset,difficulty=difficulty)
        for name,value in [('min_aligned_bp',0),('min_aligned_bp',-1),('min_mapq',-1),('min_mapq',256),
                           ('max_nm_rate',-.01),('max_nm_rate',float('nan')),('max_nm_rate',float('inf'))]:
            with self.assertRaises(ValueError):resolve(**{name:value})
        self.assertEqual(resolve(min_aligned_bp=90,min_mapq=0,max_nm_rate=0)[2],(90,0,0))

    def test_cli_preset_and_overrides_have_recorded_effective_values(self):
        cases=[([], 'standard',(100,20,.08)),(['--preset','relax2'],'relax2',(60,5,.12)),
               (['--preset','relax2','--min-aligned-bp','90','--min-mapq','12','--max-nm-rate','0.09'],'relax2',(90,12,.09))]
        with tempfile.TemporaryDirectory() as t:
            root=Path(t);host=root/'host.mmi';host.write_text('test');self.make_assets(root)
            r1,r2=root/'R1.fq',root/'R2.fq'
            for r in (r1,r2):r.write_text('@test\nACGT\n+\nIIII\n')
            for i,(flags,preset,thresholds) in enumerate(cases):
                out=root/f'out{i}'
                args=['run','--sample','unit','--r1',str(r1),'--r2',str(r2),'--out',str(out),'--reference-root',str(root),*flags]
                with patch.dict(os.environ,{'MYCOWES_HOST_MMI':str(host),'CORE_MIN_ALIGNED_BP':'7','CORE_MIN_MAPQ':'1','CORE_MAX_NM_RATE':'0.99'}), \
                     patch.object(mycowes.subprocess,'run',side_effect=self.mock_complete_engine) as run:
                    self.assertEqual(mycowes.main(args),0)
                env=run.call_args.kwargs['env'];selection=json.loads((out/'reference_selection.json').read_text())
                self.assertEqual(tuple(env[k] for k in ('CORE_MIN_ALIGNED_BP','CORE_MIN_MAPQ','CORE_MAX_NM_RATE')),tuple(map(str,thresholds)))
                self.assertEqual(selection['preset'],preset);self.assertEqual(selection['effective_min_aligned_bp'],thresholds[0])
                self.assertEqual(selection['quality_thresholds'],dict(zip(('min_aligned_bp','min_mapq','max_nm_rate'),thresholds)))
                self.assertEqual(set(selection['config_files_sha256']),{'default.env','relax2.env'})
                self.assertTrue(run.call_args.args[0][run.call_args.args[0].index('--config')+1].endswith('relax2.env' if preset=='relax2' else 'default.env'))

    def test_conflicting_reference_mode_fails_before_output(self):
        with contextlib.redirect_stderr(io.StringIO()),self.assertRaises(SystemExit):
            mycowes.main(['run','--sample','test','--r1','missing','--r2','missing','--out','missing','--difficulty','hard','--reference-mode',MODES[0]])

    def test_manuscript_preset_and_thread_arguments_reach_engine(self):
        with tempfile.TemporaryDirectory() as t:
            p=Path(t);host=p/'host.mmi';host.write_text('test');self.make_assets(p,MODES[1])
            for mate in ('r1','r2'):(p/mate).write_text('@test\nACGT\n+\nIIII\n')
            args=['run','--sample','unit','--r1',str(p/'r1'),'--r2',str(p/'r2'),
                  '--out',str(p/'out'),'--reference-root',str(p),'--preset','manuscript',
                  '--threads','2','--samtools-threads','1']
            with patch.dict(os.environ,{'MYCOWES_HOST_MMI':str(host)}), patch.object(mycowes.subprocess,'run',side_effect=self.mock_complete_engine) as engine:
                self.assertEqual(mycowes.main(args),0)
            env=engine.call_args.kwargs['env']
            self.assertEqual(env['CORE_MAPQ_POLICY'],'assembly-top-ties-gap10')
            self.assertEqual(env['MYCOWES_RISK_POLICY'],'low-complexity-nonrisk180-pair')
            self.assertEqual((env['THREADS'],env['SAMTOOLS_THREADS']),('2','1'))
            selection=json.loads((p/'out/reference_selection.json').read_text())
            self.assertEqual(selection['mode'],MODES[1])
            self.assertEqual(selection['software_version'],'MycoWES 1.0.0')

    def test_invalid_thread_arguments_do_not_start_engine(self):
        with tempfile.TemporaryDirectory() as t:
            p=Path(t);self.make_assets(p);host=p/'host.mmi';host.write_text('test')
            for mate in ('r1','r2'):(p/mate).write_text('@test\nACGT\n+\nIIII\n')
            base=['run','--sample','unit','--r1',str(p/'r1'),'--r2',str(p/'r2'),
                  '--out',str(p/'out'),'--reference-root',str(p)]
            for flags in (['--threads','0'],['--samtools-threads','-1']):
                with patch.dict(os.environ,{'MYCOWES_HOST_MMI':str(host)}), patch.object(mycowes.subprocess,'run') as engine, contextlib.redirect_stderr(io.StringIO()), self.assertRaises(SystemExit):
                    mycowes.main(base+flags)
                engine.assert_not_called()
                self.assertFalse((p/'out').exists())

    def test_cli_help_and_version(self):
        for args in (['--help'], ['--version'], ['run', '--help'], ['reference', '--help']):
            result = subprocess.run([sys.executable, str(ROOT/'mycowes.py'), *args], capture_output=True, text=True)
            self.assertEqual(result.returncode, 0, result.stderr)

    def test_cli_reference_is_same_resolver(self):
        with contextlib.redirect_stdout(io.StringIO()) as out:
            self.assertEqual(mycowes.main(['reference', '--reference-root', 'folder with spaces']), 0)
        self.assertEqual(json.loads(out.getvalue())['assets'], bindings('folder with spaces'))

    def test_run_binds_both_modes_without_changing_host(self):
        with tempfile.TemporaryDirectory() as t:
            p=Path(t);host=p/'host.mmi';host.write_text('test');r1=p/'R1.fq';r2=p/'R2.fq'
            for r in (r1,r2):r.write_text('@test\nACGT\n+\nIIII\n')
            for difficulty in ('simple','normal','hard'):
                mode=mycowes.REFERENCE_MODES[difficulty]
                self.make_assets(p, mode)
                args=['run','--sample','synthetic','--r1',str(r1),'--r2',str(r2),'--out',str(p/('output_'+difficulty)), '--reference-root',str(p),'--difficulty',difficulty]
                with patch.dict(os.environ, {'MYCOWES_HOST_MMI':str(host)}), patch.object(mycowes.subprocess,'run',side_effect=self.mock_complete_engine) as run:
                    self.assertEqual(mycowes.main(args),0)
                    call=run.call_args;self.assertEqual(call.kwargs['env']['MYCOWES_HOST_MMI'],str(host))
                    self.assertEqual(call.kwargs['env']['MYCOWES_JOINT_MMI'],bindings(p,mode)['MYCOWES_JOINT_MMI'])
                    for key,val in zip(('CORE_MIN_ALIGNED_BP','CORE_MIN_MAPQ','CORE_MAX_NM_RATE'),mycowes.DIFFICULTIES[difficulty]):
                        self.assertEqual(call.kwargs['env'][key],str(val))
                    self.assertNotIn('--host-depleted',call.args[0])
                    failed_args=args.copy()
                    failed_args[failed_args.index('--out')+1]=str(p/('failed_'+difficulty))
                    run.side_effect=lambda command,env: subprocess.CompletedProcess(command,17)
                    self.assertEqual(mycowes.main(failed_args),17)

    def test_reusing_output_with_other_profile_rejected(self):
        with tempfile.TemporaryDirectory() as t:
            p=Path(t);host=p/'h';host.write_text('test')
            for m in MODES:self.make_assets(p,m)
            for name in ('r1','r2'):(p/name).write_text('test')
            args=['run','--sample','synthetic','--r1',str(p/'r1'),'--r2',str(p/'r2'),'--out',str(p/'out'),'--reference-root',str(p)]
            with patch.dict(os.environ, {'MYCOWES_HOST_MMI':str(host)}), patch.object(mycowes.subprocess,'run',side_effect=self.mock_complete_engine):
                mycowes.main(args)
                with contextlib.redirect_stderr(io.StringIO()),self.assertRaises(SystemExit):mycowes.main(args+['--difficulty','hard'])

    def test_no_reference_or_analysis_payload(self):
        forbidden={'.fa','.fasta','.bed','.mmi','.bam','.fastq','.fq','.xlsx','.png','.pdf','.docx','.log','.pyc'}
        for p in ROOT.rglob('*'):
            self.assertNotIn(p.suffix.lower(), forbidden, str(p))

    def test_core_pair_dedup_synthetic(self):
        with tempfile.TemporaryDirectory() as t:
            p=Path(t)
            for name in ('r1','r2'):
                (p/name).write_text('@a\nACGT\n+\nIIII\n@b\nACGT\n+\nIIII\n')
            result=subprocess.run([sys.executable,str(ROOT/'workflow/exact_pair_dedup_v0.2.py'),*[str(p/n) for n in ('r1','r2','o1','o2','groups')]],capture_output=True,text=True)
            self.assertEqual(result.returncode,0,result.stderr)
            self.assertEqual((p/'o1').read_text().count('@'),1)
            self.assertIn('a\t2\ta,b',(p/'groups').read_text())


if __name__=='__main__':unittest.main()
