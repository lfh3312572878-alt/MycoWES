#!/usr/bin/env python3
"""Copy the package with its default8/32/8/4 execution profile unchanged."""
from pathlib import Path
import argparse,hashlib,json,shutil,sys
SOURCE=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(SOURCE))
from mycowes import CPU_PROFILES
def engine(root):
 files=[root/'mycowes.py',*sorted((root/'config').glob('*.env')),root/'VERSION',*sorted((root/'workflow').glob('*.py')),*sorted((root/'workflow').glob('*.sh'))]
 return hashlib.sha256(b''.join(p.read_bytes() for p in files)).hexdigest()
def main():
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--out',type=Path,required=True);a=p.parse_args();out=a.out.resolve()
 if out==SOURCE or SOURCE in out.parents or out.exists():p.error('Choose a new directory outside the source package.')
 before=engine(SOURCE)
 if before!=json.loads((SOURCE/'PACKAGE_IDENTITY.json').read_text())['canonical_engine_sha256']:p.error('Engine differs from this package identity.')
 shutil.copytree(SOURCE,out)
 assert engine(SOURCE)==before==engine(out)
 (out/'RESOURCE_EXECUTION_PROFILE.json').write_text(json.dumps(dict(canonical_engine_sha256=before,execution_engine_sha256=engine(out),cpu_profile='benchmark',threads=CPU_PROFILES['benchmark'],anchor_mapper_threads=2,scientific_rules_changed=False),indent=2)+'\n')
 files=sorted(f for f in out.rglob('*') if f.is_file() and f.name!='SHA256SUMS')
 (out/'SHA256SUMS').write_text(''.join(hashlib.sha256(f.read_bytes()).hexdigest()+'  '+f.relative_to(out).as_posix()+'\n' for f in files))
 print('Created execution profile:',out)
if __name__=='__main__':main()
