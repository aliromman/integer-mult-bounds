"""Check the whole-width ordinary S3 phase and all three characters.

The proof uses Fourier eigenvalues and a quadratic Gauss sum. These
literal integer products are finite controls, not an all-size replay.
"""
import json


def mm(a, b):
    columns = list(zip(*b))
    return [[sum(x*y for x, y in zip(row, column)) for column in columns] for row in a]


def quadratic(x):
    weight = x.bit_count()
    return (weight*(weight-1)//2) % 2


def component(width, remove_radical=True):
    size = 1 << width
    even = width-width % 2
    denominator = 1 << (even//2)
    remove = int(width % 4 == 3 and remove_radical)
    phases = [(-1)**(quadratic(x)+remove*((x >> (width-1)) & 1)) for x in range(size)]
    walsh = [sum(value*(-1)**((x & y).bit_count() % 2) for y, value in enumerate(phases))
             for x in range(size)]
    divisor = size//denominator
    if any(value % divisor for value in walsh):
        raise ValueError('The kernel denominator differs from the quadratic rank')
    kernel = [value//divisor for value in walsh]
    mask = (1 << even)-1
    linear = (even//2+1) % 2
    ordinary = [(-1)**(quadratic(x & mask)+linear*((x & mask).bit_count() % 2))
                for x in range(size)]
    s = [[kernel[x ^ y] for y in range(size)] for x in range(size)]
    sj = [[value*ordinary[y] for y, value in enumerate(row)] for row in s]
    return denominator, s, ordinary, sj


def is_scaled_identity(matrix, scale):
    return all(value == (scale if i == j else 0)
               for i, row in enumerate(matrix) for j, value in enumerate(row))


def verify():
    records = []
    for width in range(1, 9):
        denominator, s, j, sj = component(width)
        if not is_scaled_identity(mm(s, s), denominator**2):
            raise ValueError('The query square differs from the identity')
        if not is_scaled_identity(mm(mm(sj, sj), sj), denominator**3):
            raise ValueError('The complete S3 relation fails')
        size = len(s)
        trace_s = sum(s[i][i] for i in range(size))//denominator
        trace_product = sum(sj[i][i] for i in range(size))//denominator
        if trace_s != sum(j):
            raise ValueError('The two transposition characters differ')
        standard = (size-trace_product)//3
        one_dimensional = (size+2*trace_product)//3
        trivial = (one_dimensional+trace_s)//2
        sign = (one_dimensional-trace_s)//2
        if trivial+sign+2*standard != size or min(trivial, sign, standard) < 0:
            raise ValueError('The character multiplicities fail')
        if width >= 4 and min(trivial, sign, standard) == 0:
            raise ValueError('The all-character control fails')
        # Q G has the same Fourier signs. Retain its phase as an integer
        # exponent of i, modulo four.
        for x in range(size):
            weight = x.bit_count()
            if (weight-weight % 2) % 4 != 2*quadratic(x):
                raise ValueError('The Q-to-involution phase differs')
        records.append({'width': width, 'dimension': size, 'denominator': denominator,
                        'trivial': trivial, 'sign': sign, 'standard': standard,
                        'square_and_cube': 'PASS'})
    denominator, _, _, sj = component(3, remove_radical=False)
    if is_scaled_identity(mm(mm(sj, sj), sj), denominator**3):
        raise ValueError('The changed odd-width normalization was not rejected')
    denominator, s, _, _ = component(4)
    if is_scaled_identity(mm(mm(s, s), s), denominator**3):
        raise ValueError('The omitted ordinary phase was not rejected')
    result = {'status': 'PASS_WHOLE_WIDTH_S3_COMPONENT', 'widths': records,
              'changed_normalization_rejected': True, 'new_kappa': False,
              'scope': 'Literal finite component checks; no full source route, caller or price'}
    print(json.dumps(result, indent=2))
    return result


if __name__ == '__main__':
    verify()
