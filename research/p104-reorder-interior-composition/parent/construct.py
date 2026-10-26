"""Apply four checked changes of the public register basis.

Ali Romman with substantial OpenAI Codex assistance. Apache-2.0.
"""
from copy import deepcopy
from fractions import Fraction as Q
from sources import need


def blocks(word):
    result, block, frame = [], [], None
    for event in word['program']:
        if event[0] == 'MOVE':
            continue
        if event[0] in ('INIT', 'CLEAR'):
            if block:
                result.append(block)
            result.append([event])
            block, frame = [], None
        elif frame != event[-1]:
            if block:
                result.append(block)
            block, frame = [event], event[-1]
        else:
            block.append(event)
    if block:
        result.append(block)
    return result


def add(dst, src, coefficient):
    for column, value in tuple(src.items()):
        value = dst.get(column, Q(0)) + coefficient * value
        if value:
            dst[column] = value
        else:
            dst.pop(column, None)


def matrix(events, roles, frame):
    rows = {r: {r: Q(1)} for r in roles}
    for event in events:
        op, r = event[:2]
        need(op in ('ADD', 'SCALE') and event[-1] == frame,
             'One ordinary scalar block at one chart')
        coefficient = Q(event[-2])
        need(r in roles and coefficient != 0, 'Declared scalar destination')
        if op == 'SCALE':
            rows[r] = {j: coefficient * value for j, value in rows[r].items()}
        else:
            need(event[2] in roles and event[2] != r, 'Declared scalar source')
            add(rows[r], rows[event[2]], coefficient)
    return rows


def protected(word):
    n = word['public']
    return [e for e in word['program'] if e[1] >= n
            or (e[0] in ('ADD', 'COPY') and e[2] >= n)]


def rebuild(word, scalar):
    result = deepcopy(word)
    frames = word['initial'][:] + [None] * word['owned_fields']
    program = []

    def move(r, frame):
        if frames[r] != frame:
            need(frames[r] is not None, 'Live transport source')
            program.append(['MOVE', r, frames[r], frame, r >= word['public']])
            frames[r] = frame

    for event in scalar:
        op, r = event[:2]
        if op == 'INIT':
            need(frames[r] is None, 'Unused owned field')
            frames[r] = event[2]
        elif op == 'CLEAR':
            move(r, event[2])
            frames[r] = None
        else:
            move(r, event[-1])
            if op in ('ADD', 'COPY'):
                move(event[2], event[-1])
        program.append(event[:])
    for r in range(word['public']):
        move(r, word['final'][r])
    result['program'] = program
    return result


def step(word, witness, local):
    source = blocks(word)
    a, b, q = witness['a'], witness['b'], Q(witness['coefficient'])
    start, end = witness['start_block'], witness['end_block']
    entry, exit_frame = witness['entry_frame'], witness['exit_frame']
    need(0 <= a < word['public'] and 0 <= b < word['public'] and a != b,
         'Two distinct public fields')
    need(q and q.denominator & (q.denominator - 1) == 0, 'Nonzero dyadic basis')
    need(0 < start < end < len(source), 'Internal interval of scalar blocks')
    need(0 <= entry < len(word['frames']) and 0 <= exit_frame < len(word['frames']),
         'Declared basis charts')
    replacements = {row['index']: row for row in witness['blocks']}
    need(len(replacements) == len(witness['blocks']), 'Distinct replacement blocks')
    expected = set()
    scalar = []
    for index, block in enumerate(source):
        if index == start:
            scalar.append(['ADD', a, b, str(q), entry])
        roles = {e[1] for e in block} | {
            e[2] for e in block if e[0] in ('ADD', 'COPY')}
        if start <= index < end and roles & {a, b}:
            expected.add(index)
            need(index in replacements, 'Every affected block has a certificate')
            row = replacements[index]
            need(block == row['before'], 'Exact preceding scalar block')
            frame = block[0][-1]
            roles |= {a, b}
            roles |= {e[1] for e in row['after']} | {
                e[2] for e in row['after'] if e[0] == 'ADD'}
            original = matrix(block, roles, frame)
            # R = I + q E_ab. The two operations form R M R^-1.
            add(original[a], original[b], q)
            for vector in original.values():
                add(vector, {b: Q(1)}, -q * vector.get(a, 0))
            need(matrix(row['after'], roles, frame) == original,
                 'Exact conjugation of the complete scalar block')
            scalar.extend(row['after'])
        else:
            scalar.extend(block)
        if index + 1 == end:
            scalar.append(['ADD', a, b, str(-q), exit_frame])
    need(set(replacements) == expected, 'No extra replacement certificate')
    result = rebuild(word, scalar)
    need(protected(result) == protected(word), 'Every owned instruction and transport stays fixed')
    need(result['initial'] == word['initial'] and result['final'] == word['final'],
         'All public boundaries stay fixed')
    return result, local.replay(result)


def construct(source, witnesses, local):
    word, receipts = source, []
    for witness in witnesses:
        word, receipt = step(word, witness, local)
        receipts.append(receipt)
    word['schema'] = 'p104-paid-interior-bases/1'
    return word, receipts
