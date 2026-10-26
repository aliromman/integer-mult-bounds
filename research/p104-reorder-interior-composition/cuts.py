"""Check the full common-frame helper cut domain after composition.

The kernel method follows PR5599 and PR5609. This version tests every
operation boundary and permits a distinct physical cut role. Ali Romman
with substantial OpenAI Codex assistance. Apache-2.0.
"""
from collections import defaultdict
from fractions import Fraction as Q
import local as L

def pure(rows,roles,group):
 pivots={};out=[]
 for role in roles:
  outside={k:v for k,v in rows[role].items()if k not in group};inside={k:v for k,v in rows[role].items()if k in group};coeff={role:Q(1)}
  while outside:
   p=min(outside)
   if p not in pivots:
    scale=1/outside[p];pivots[p]=[{k:v*scale for k,v in row.items()}for row in [outside,inside,coeff]];break
   scale=-outside[p]
   for target,source in zip([outside,inside,coeff],pivots[p]):L.add(target,source,scale)
  if not outside:out.append((inside,coeff))
 return out

def dot(a,b):return sum(v*b.get(k,0)for k,v in a.items())
def combo(rows,co):
 result={}
 for i,v in co.items():L.add(result,rows[i],v)
 return result
def scan(word,group,entry,exit_frame):
 spans=[L.span(b)for b in word['frames']];valid={i for i,f in enumerate(spans)if spans[entry]<f<spans[exit_frame]};public=word['public'];n=public+word.get('owned_fields',8)
 if not valid:return None
 rows=[{i:Q(1)}for i in range(n)];st=word['initial'][:]+[None]*(n-public);banks=defaultdict(set)
 for i,f in enumerate(st):banks[f].add(i)
 snapshots={};checked=0;nonzero_pairings=0
 def snapshot(ix,frames):
  nonlocal checked
  saved=[]
  for f in frames&valid:
   if not banks[f]:continue
   roles=sorted(banks[f]);basis=pure(rows,roles,group);checked+=1
   if any(c for c,a in basis):saved.append((f,roles,basis))
  if saved:snapshots[ix]=saved
 snapshot(-1,set(banks))
 for ix,e in enumerate(word['program']):
  op,r=e[:2];affected={st[r]}
  if op=='MOVE':
   banks[st[r]].remove(r);st[r]=e[3];banks[st[r]].add(r);affected.add(st[r])
  elif op in ['INIT','CLEAR']:
   rows[r]={};banks[st[r]].remove(r);st[r]=e[2]if op=='INIT'else None;banks[st[r]].add(r);affected.add(st[r])
  elif op=='COPY':rows[r]={k:v*Q(e[-2])for k,v in rows[e[2]].items()}
  elif op=='SCALE':rows[r]={k:v*Q(e[-2])for k,v in rows[r].items()}
  else:L.add(rows[r],rows[e[2]],Q(e[-2]))
  snapshot(ix,affected)
 columns=[{i:Q(1)}for i in range(n)]
 for ix in range(len(word['program'])-1,-2,-1):
  for f,roles,left in snapshots.get(ix,[]):
   right=pure(columns,roles,group)
   if not any(d for d,b in right):continue
   nonzero_pairings+=any(dot(c,d) for c,a in left for d,b in right)
  if ix<0:break
  e=word['program'][ix];op,r=e[:2]
  if op in ['INIT','CLEAR']:columns[r]={}
  elif op=='COPY':L.add(columns[e[2]],columns[r],Q(e[-2]));columns[r]={}
  elif op=='SCALE':columns[r]={k:v*Q(e[-2])for k,v in columns[r].items()}
  elif op=='ADD':L.add(columns[e[2]],columns[r],Q(e[-2]))
 L.need(nonzero_pairings==0,'No independent helper pairing remains')
 return dict(entry=entry,exit=exit_frame,helpers=len(group),contexts=checked,nonzero_pairings=nonzero_pairings)

def run(word):
    main=set(word['live'])|{i+112 for i in word['live']};groups=defaultdict(set)
    for i in range(word['public']):
        if i not in main:groups[word['initial'][i],word['final'][i]].add(i)
    rows=[]
    for (a,b),g in sorted(groups.items()):
        result=scan(word,g,a,b)
        if result is not None:rows.append(result)
    return dict(boundary_groups=len(groups),eligible_groups=len(rows),contexts=sum(r['contexts']for r in rows),nonzero_pairings=0,groups=rows)
