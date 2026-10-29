"""Check the three-query phase model with exact integer arithmetic.

The all-width component proof is in COMPONENT.md. These finite checks do
not bind a complete native route or caller.
"""
import json


def require(ok,message):
    if not ok:raise ValueError(message)


def walsh(values):
    out=list(values);width=1
    while width<len(out):
        for first in range(0,len(out),2*width):
            for i in range(first,first+width):
                a,b=out[i],out[i+width];out[i],out[i+width]=a+b,a-b
        width*=2
    return out


def binary_rank(rows):
    pivots={}
    for r in rows:
        while r:
            p=r.bit_length()-1
            if p not in pivots:pivots[p]=r;break
            r^=pivots[p]
    return len(pivots)


def normalize_checks():
    count=0
    for m in range(1,7):
        images=set()
        for y in range(1<<(2*m)):
            bits=[(y>>j)&1 for j in range(2*m)]
            ab=[]
            for k in range(m):
                tail=sum(bits[2*k+2:])%2;ab.extend([bits[2*k]^tail,bits[2*k+1]^tail])
            hyper=sum(ab[2*k]*ab[2*k+1] for k in range(m))%2
            ell=sum((k%2)*(ab[2*k]+ab[2*k+1]) for k in range(m))%2
            q=(sum(bits)*(sum(bits)-1)//2)%2
            require(q==(hyper^ell),'The native phase normalization differs')
            # The native Q eigenvalue times G must have the same absolute sign.
            weight=sum(bits);power=(weight-(weight%2))%4
            require(power==2*q,'The absolute Q/G phase differs')
            images.add(tuple(ab));count+=1
        require(len(images)==1<<(2*m),'The address change is not invertible')
    return count


def component_checks():
    records=[]
    for m in range(1,4):
        size=1<<(3*m);mask=(1<<m)-1
        phase=[[],[],[]];counts={(1,1):0,(-1,1):0,(1,-1):0,(-1,-1):0}
        for y in range(size):
            a=y&mask;b=(y>>m)&mask;c=y>>(2*m)
            q=[(a&b).bit_count()%2,(a&c).bit_count()%2,(a&(b^c)).bit_count()%2]
            signs=[1-2*x for x in q]
            for i in range(3):phase[i].append(signs[i])
            require(signs[0]*signs[1]==signs[2],'The query product differs')
            counts[tuple(signs[:2])]+=1
        trace=1<<(2*m)
        require(counts[(1,1)]==(size+3*trace)//4,'The positive multiplicity differs')
        require(all(v==(size-trace)//4 for k,v in counts.items() if k!=(1,1)),'A nontrivial multiplicity differs')
        kernels=[walsh(p) for p in phase]
        for shift in range(size):
            product=sum(kernels[0][x]*kernels[1][x^shift] for x in range(size))
            require(product==size*kernels[2][shift],'The literal convolution endpoint differs')
        ranks=[]
        for p in phase:
            def q(y):return (1-p[y])//2
            rows=[]
            for i in range(3*m):
                rows.append(sum((q((1<<i)^(1<<j))^q(1<<i)^q(1<<j))<<j for j in range(3*m)))
            ranks.append(binary_rank(rows))
        require(ranks==[2*m]*3,'A query polar rank differs')
        changed=list(kernels[2]);changed[0]+=1
        require(sum(kernels[0][x]*kernels[1][x] for x in range(size))!=size*changed[0], 'A changed kernel passes')
        records.append({'m':m,'dimension':size,'trace':trace,'multiplicities':list(counts.values()),'polar_ranks':ranks})
    return records


if __name__=='__main__':
    result={'normalization_inputs':normalize_checks(),'components':component_checks(),'changed_kernel_rejected':True,'new_kappa':False}
    print(json.dumps(result,indent=2));print('PASS_KLEIN_QUERY_COMPONENT_ONLY')
