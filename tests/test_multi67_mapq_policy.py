import unittest,importlib.util
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
def load(n,p):
 s=importlib.util.spec_from_file_location(n,p);m=importlib.util.module_from_spec(s);s.loader.exec_module(m);return m
parser=load('tie_parser',ROOT/'workflow/parse_fragment_evidence_v0.1.py')
risk=load('tie_quality',ROOT/'workflow/annotate_fragment_risk_v0.1.py')
def hit(asm,**kw):
 r=dict(read_end='R1',mapped=True,supplementary=False,reference_class='TARGET',species_taxid='1',assembly=asm,AS=300,reference=asm+'_contig',position=10,reference_end=160,reverse=False)
 return dict(r,**kw)
class TestAssemblyTies(unittest.TestCase):
 def test_top_gap10_boundaries(self):
  a,b=hit('a'),hit('b')
  candidate=lambda other: parser.assembly_tie_only([a,b,other],a,True)
  self.assertTrue(candidate(hit('c',AS=290,species_taxid='2')))
  self.assertTrue(candidate(hit('c',AS=290,reference_class='COMPETITOR')))
  self.assertTrue(candidate(hit('a',AS=290,position=20)))
  for other in [hit('c',AS=291,species_taxid='2'),hit('c',AS=300,species_taxid='2'),
                hit('a',position=20),hit('c',AS=290,supplementary=True),
                hit('c',AS=None),hit('c',species_taxid=''),hit('c',assembly='')]:
   self.assertFalse(candidate(other))
  self.assertFalse(parser.assembly_tie_only([a,hit('b',AS=299)],a,True))
  self.assertFalse(parser.assembly_tie_only([a,b,hit('c',AS=301)],a,True))
  self.assertFalse(parser.assembly_tie_only([a]+[b]*100,a,True))
 def test_top_gap10_quality_guard(self):
  r=dict(pair_support_class='paired_target_consistent',R1_primary_aligned_bases='150',R2_primary_aligned_bases='150',R1_primary_MAPQ='0',R2_primary_MAPQ='0',maximum_primary_NM_rate='0',R1_assembly_top_tie_gap10='1',R2_assembly_top_tie_gap10='1',R1_primary_species='s',R2_primary_species='s',R1_primary_class='TARGET',R2_primary_class='TARGET',R1_best_target_species='s',R2_best_target_species='s')
  self.assertEqual(risk.core_filter_reason(r,100,mapq_policy='assembly-top-ties-gap10'),'PASS')
  for key,value in [('R2_assembly_top_tie_gap10','0'),('R1_primary_aligned_bases','99'),('maximum_primary_NM_rate','0.08001'),('R2_primary_species','x'),('R2_best_target_species','s;x'),('R2_primary_class','COMPETITOR')]:
   self.assertNotEqual(risk.core_filter_reason(dict(r,**{key:value}),100,mapq_policy='assembly-top-ties-gap10'),'PASS')
 def test_boundaries(self):
  a=hit('a');b=hit('b')
  self.assertTrue(parser.assembly_tie_only([a,b],a))
  for bad in [hit('b',species_taxid='2'),hit('b',reference_class='COMPETITOR'),hit('b',assembly=''),hit('b',supplementary=True),hit('b',AS=299)]:
   self.assertFalse(parser.assembly_tie_only([a,bad],a))
  self.assertFalse(parser.assembly_tie_only([a,b,hit('a',position=20)],a))
  self.assertFalse(parser.assembly_tie_only([a,b,hit('c',species_taxid='2',AS=100)],a))
  self.assertFalse(parser.assembly_tie_only([a]+[b]*100,a))
 def test_gate(self):
  r=dict(pair_support_class='paired_target_consistent',R1_primary_aligned_bases='150',R2_primary_aligned_bases='150',R1_primary_MAPQ='0',R2_primary_MAPQ='0',maximum_primary_NM_rate='0',R1_assembly_tie_only='1',R2_assembly_tie_only='1',R1_primary_species='species1',R2_primary_species='species1',R1_primary_class='TARGET',R2_primary_class='TARGET',R1_best_target_species='species1',R2_best_target_species='species1')
  self.assertEqual(risk.core_filter_reason(r,100),'MAPQ_LT20')
  self.assertEqual(risk.core_filter_reason(r,100,mapq_policy='assembly-ties-only'),'PASS')
  for key,value in [('R2_assembly_tie_only','0'),('R1_primary_aligned_bases','99'),('maximum_primary_NM_rate','0.09'),('R2_primary_species','species2'),('R1_best_target_species','species1;species2'),('R2_primary_class','COMPETITOR')]:
   self.assertNotEqual(risk.core_filter_reason(dict(r,**{key:value}),100,mapq_policy='assembly-ties-only'),'PASS')
if __name__=='__main__':unittest.main()
