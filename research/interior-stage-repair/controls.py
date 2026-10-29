"""Compare the affine insertion formula with literal rational words."""
from fractions import Fraction as F
from importlib.util import module_from_spec,spec_from_file_location
from pathlib import Path
spec=spec_from_file_location('interior_repair',Path(__file__).with_name('discovery.py'))
m=module_from_spec(spec);spec.loader.exec_module(m)


def evaluate(word,z):
    n=len(word[0]);out=[[F(0)]*n for _ in range(n)]
    for degree,a in enumerate(word):out=m.add(out,[[F(x)*z**degree for x in row] for row in a])
    return out


def verify():
    cases=0
    for n in (4,6):
        channels=m.setup(n);ident=[[F(i==j) for j in range(n)] for i in range(n)]
        tests=[ident,[[F((i+2*j)%5-2,3) for j in range(n)] for i in range(n)]]
        for cut in range(n+1):
            for channel,(prefix,suffix) in enumerate(channels):
                linear=m.linear_product(suffix[cut+1],prefix[cut+1],n)
                for jmat in tests:
                    coeff=[[[sum(c*jmat[w[0]//n][w[0]%n] for w,c in p.items()) for p in row] for row in a] for a in linear]
                    for z in (F(-1),F(2),F(3,2)):
                        left,right=evaluate(suffix[cut+1],z),evaluate(prefix[cut+1],z)
                        if evaluate(coeff,z)!=m.mm(m.mm(left,jmat),right):raise ValueError('The affine insertion differs from the literal product')
                        for inverse in (False,True):
                            shift=(1/z if inverse else z)-1
                            factor=m.add(ident,[[shift*x for x in row] for row in jmat])
                            literal=m.mm(m.mm(left,factor),right)
                            affine=m.add(evaluate(prefix[-1],z),[[shift*x for x in row] for row in evaluate(coeff,z)])
                            if literal!=affine:raise ValueError('The inserted Laurent word differs')
                            cases+=1
    # A scalar projector preserves the query identity at every insertion cut.
    for n in (4,6):
        plus,minus=m.setup(n)
        for z in (F(2),F(3)):
            if any(x for row in m.add(evaluate(plus[0][-1],z),evaluate(minus[0][-1],z)) for x in row):raise ValueError('The query control fails')
    print('PASS_LITERAL_INTERNAL_STAGE_FORMULA',cases)


if __name__=='__main__':verify()
