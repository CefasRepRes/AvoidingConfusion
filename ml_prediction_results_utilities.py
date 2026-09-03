#!/usr/bin/env python3
from pathlib import Path
import json,re
import pandas as pd
META={"timestamp","run_name","blob_name"}
def load_json_payload(s,**_):
 if isinstance(s,(dict,list)): return s
 return json.loads(Path(s).read_text())
def load_inference_rows_from_source(s,**_):
 p=load_json_payload(s); p=p.get('rows',p); rows=[]
 for r in p:
  rows.append({'timestamp':pd.to_datetime(r['timestamp'],utc=True),'run_name':r.get('run_name',str(s)),'blob_name':r.get('blob_name',str(s)),**{k:int(v) for k,v in r['class_counts'].items()}})
 return rows
def build_dataframe(rows):
 d=pd.DataFrame(rows); cs=sorted(c for c in d if c not in META)
 for c in cs: d[c]=d[c].fillna(0).astype(int)
 d=d.sort_values(['timestamp','run_name','blob_name']).reset_index(drop=True); return d[['timestamp','run_name','blob_name']+cs],d
def safe_part(x): return re.sub(r'[^A-Za-z0-9._-]+','_',x).strip('_') or 'value'
def write_long_counts_csv(d,p):
 cs=[c for c in d if c not in META]; d.melt(id_vars=['timestamp','run_name','blob_name'],value_vars=cs,var_name='class',value_name='count').to_csv(p,index=False)
def write_long_uncertainty_csv(d,p):
 names=sorted({c[:-15] for c in d if c.endswith('_corrected_mean')}); rows=[]
 for _,r in d.iterrows():
  for n in names: rows.append({'timestamp':r.timestamp,'run_name':r.run_name,'blob_name':r.blob_name,'class':n,**{m:r.get(f'{n}_{m}') for m in ('corrected_mean','corrected_median','corrected_lower','corrected_upper')}})
 pd.DataFrame(rows).to_csv(p,index=False)
