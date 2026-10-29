"""Check the literal blocks used in the all-length curvature proof."""
import json
from boundary_curvature import fixed, e


def multiply(a,b):
    return [[sum(x*y for x,y in zip(row,col)) for col in zip(*b)] for row in a]


def add(a,b,sign=1):
    return [[x+sign*y for x,y in zip(row,other)] for row,other in zip(a,b)]


def verify():
    rows=[]
    for n in (4,6,8,10,12,16,24,32):
        m=n//2;h=[[0]*n for _ in range(n)];c=[[0]*n for _ in range(n)]
        ident=[[int(i==j) for j in range(n)] for i in range(n)];b=ident
        for poly in fixed(n):
            p=[[int(x.get((),0)) for x in row] for row in poly]
            c=add(c,add(multiply(h,p),multiply(p,h),-1));h=add(h,p)
            b=multiply(add(ident,p,-1),b)
        frame=[[0]*n for _ in range(n)];inverse=[[0]*n for _ in range(n)]
        for k in range(m):
            frame[2*k][k]=1;frame[2*k+1][k]=-1;frame[2*k+1][m+k]=1
            inverse[k][2*k]=1;inverse[m+k][2*k]=inverse[m+k][2*k+1]=1
        transformed=multiply(multiply(inverse,c),frame)
        for i in range(m):
            for j in range(m):
                if transformed[i][j]!=(j-i if j>i else 0):raise ValueError('The upper curvature block differs')
                if transformed[m+i][m+j]!=(-(i-j) if i>j else 0):raise ValueError('The lower curvature block differs')
                if transformed[m+i][j]:raise ValueError('The lower left block is nonzero')
        if any(row[0] for row in b):raise ValueError('The negative extreme query kernel fails')
        last=e.canonical(n)[1][-1]
        if any(sum(last[i]*b[i][j] for i in range(n)) for j in range(n)):raise ValueError('The final negative extreme kernel fails')
        rows.append({'n':n,'minor_size':n-2,'minor_determinant':(-1)**(m-1),'negative_extreme_kernel':'PASS'})
    return rows


if __name__=='__main__':
    print(json.dumps(verify(),indent=2))
    print('PASS_ALL_LENGTH_CURVATURE_BLOCK_CONTROLS')
