"""Apply the released helper473 changes to the paid-interior word.

The edits and trace criterion are from Joel Pulikkan (GamingPuzzled),
with Claude Opus 5.5 assistance, in PR5666. This checked composition is
by Ali Romman with substantial OpenAI Codex assistance. Apache-2.0.
"""
from collections import Counter,defaultdict
from copy import deepcopy
from pathlib import Path
from sources import load,need
HERE=Path(__file__).resolve().parent

def construct(base,schema='p104-interior-r22-composition/1'):
    edits=load(HERE/'reorder/construct.py','composition_released_edits')
    word=deepcopy(base)
    need(word['initial'][edits.HELPER]==edits.OLD_START,'Exact source helper start')
    word['initial'][edits.HELPER]=edits.NEW_START
    spots=[]
    for before,after in edits.EDITS:
        hits=[i for i in range(len(word['program'])-len(before)+1)
              if word['program'][i:i+len(before)]==before]
        need(len(hits)==1,'One exact source run for each released edit')
        spots.append((hits[0],before,after))
    occupied=[i+j for i,before,after in spots for j in range(len(before))]
    need(len(set(occupied))==len(occupied),'Disjoint source runs')
    for i,before,after in sorted(spots,reverse=True):
        word['program'][i:i+len(before)]=deepcopy(after)
    word['schema']=schema
    return word

def trace_equivalent(base,word):
    public=base['public']
    def kinds(e):
        op=e[0]
        if op=='ADD':out=[(e[1],'T'),(e[2],'S')]
        elif op=='COPY':out=[(e[1],'X'),(e[2],'S')]
        elif op in ['INIT','CLEAR']:return [(e[1],'X'),('owned-lifetime','X')]
        else:out=[(e[1],'X')]
        return [(r,'X' if r>=public else t)for r,t in out]
    def strip(e):return tuple(e[:2]) if e[0]in ['INIT','CLEAR'] else tuple(e[:-1])
    old=[e for e in base['program']if e[0]!='MOVE']
    new=[e for e in word['program']if e[0]!='MOVE']
    need(Counter(map(strip,old))==Counter(map(strip,new)),'Same scalar and owned operations')
    slots=defaultdict(list)
    for j,e in enumerate(new):slots[strip(e)].append(j)
    for v in slots.values():v.reverse()
    match=[slots[strip(e)].pop()for e in old];groups=defaultdict(list)
    for i,e in enumerate(old):
        for r,t in kinds(e):
            g=groups[r]
            if g and g[-1][0]==t and t!='X':g[-1][1].append(i)
            else:g.append([t,[i]])
    for g in groups.values():
        for (_,a),(_,b)in zip(g,g[1:]):
            need(max(match[i]for i in a)<min(match[i]for i in b),'Dependent operations keep their order')
    changed=[[i,old[i][:-1],old[i][-1],new[match[i]][-1]]for i in range(len(old))if old[i]!=new[match[i]]]
    return dict(scalar_and_owned_events=len(old),changed_charts=changed)
