"""Check the old stream lower bound from the source events and a second depth calculation."""
from pathlib import Path
import json
import sys

work = Path("build/construction-scale-nov21")
package = Path("research/single-band-native-carriers").resolve()
sys.path.insert(0, str(package))
from upstream import install
install()
from compiler_cut import source_graph
from tagged_sources import BIT_CALLER, BIT_PIN, packed, pinned

graph = source_graph()
caller = packed(pinned(BIT_CALLER, BIT_PIN))
frames = {int(k): v["dim"] for k,v in graph["frames"].items()}
moves = [abs(frames[e[2]]-frames[e[3]]) for e in graph["program"] if e[0]=="MOVE"]
local = sum(e[0]=="LOCAL" for e in caller["events"])
outer = [(1 if e[0]=="QUERY" else 7*len(e[2])) for e in caller["events"] if e[0] in ("QUERY","IDLE_ONE")]
ports = caller["ports"]
total = local*len(moves)+ports*len(outer)
high = local*sum(r>=14 for r in moves)+ports*sum(r>=14 for r in outer)
assert total==107460 and high==3108
assert min(moves+outer)>0
ranks = {1+7*len(e["stages"]) for e in caller["diagonal"]}
ranks.update(abs(frames[graph["final"][j]]-frames[graph["start"][j]]) for j in range(2*ports,len(graph["start"])))
assert {8,64} <= ranks

# Use directed rational logarithms with 12 fractional bits. Discard at least
# one and up to 10,000 whole exit positions per level as a separate lower estimate.
scale=4096
def log_floor(x):
    return (x**scale).bit_length()-1
A=134**480
R=11**480
root=236800000000000
ambient_upper=log_floor(A)+1
block_lower=log_floor(R)
lo=root*scale
levels=[]
for i in range(7):
    # Uniformly allow for the round down to a multiple of four.
    m=1+(lo-3*scale-ambient_upper)//(4*480*scale)
    m=(m-1)//10000*10000
    assert m>=1
    levels.append(m)
    lo=block_lower*m
count=480*(sum(levels[:6])*log_floor(high)+levels[-1]*log_floor(total))//scale
upper=json.loads((package/"REQUIREMENTS.json").read_bytes())["price-match"]["log2_all_fixed_setup_upper"]
assert count>=3366293830151418 and count>2*upper
result=dict(raw_paid=total,high_paid=high,log2_expanded_lower_separate_calculation=count,rounded_down_levels=levels)
print(json.dumps(result))
print("PASS_INDEPENDENT_SETUP_LOWER")
