"""Execution-only: preseed byte-identical read-only FASTA, then run unchanged anchor code."""
from pathlib import Path
import json,sys,runpy,os
ROOT=Path(__file__).resolve().parents[2]
cache=Path(os.environ.get('MYCOWES_ANCHOR_CACHE',str(ROOT/'anchor_cache')))
out=Path(sys.argv[sys.argv.index('--out')+1]);tmp=out/'tmp_assembly_selection'
records=json.loads((cache/'MANIFEST.json').read_text())['files'] if (cache/'MANIFEST.json').is_file() else []
for r in records:
 p=Path(r['path']);st=p.stat()
 assert [st.st_size,st.st_mtime_ns,st.st_ino]==r['identity'],'Cached reference replaced'
 d=tmp/r['safe_species'];d.mkdir(parents=True,exist_ok=True)
 link=d/(r['assembly']+'.fa');assert not link.exists(),'Unexpected pre-existing per-run reference'
 link.symlink_to(p)
runpy.run_path(str(Path(__file__).with_name('evidence_anchor_breadth_v0.2.py')),run_name='__main__')
