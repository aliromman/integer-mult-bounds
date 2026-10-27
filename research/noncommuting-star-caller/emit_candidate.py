"""Emit and replay the small frame candidate on the selected source.

The signed gain is the known helper-307 construction from PR5601.
This replay is a positive control, not a new construction.
"""
from copy import deepcopy
import gzip
import hashlib
import importlib.util
import json
from pathlib import Path
import sys

import binary_subspaces as ms

HERE = Path(__file__).resolve().parent


def load(path, name):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


def emit(source, graph, candidate):
    word = deepcopy(source)
    tags = {ms.span(b): i for i, b in enumerate(word['frames'])}
    chosen = list(graph['nodes'])
    for edit in candidate['changed']:
        node = edit['node']
        if ms.span(edit['old']) != ms.span(word['frames'][chosen[node]]):
            raise ValueError('Candidate source frame differs')
        space = ms.span(edit['new'])
        if space not in tags:
            tags[space] = len(word['frames'])
            word['frames'].append(edit['new'])
        chosen[node] = tags[space]
    word['initial'] = [chosen[v] for v in graph['entry']]
    word['final'] = [chosen[v] for v in graph['exit']]
    current = word['initial']+[None]*8
    program = []
    for ei, original in enumerate(source['program']):
        op = original[:]
        name, r = op[:2]
        target = chosen[graph['events'][ei]]
        if name == 'INIT':
            op[2] = target
            current[r] = target
        elif name == 'CLEAR':
            op[2] = current[r]
            current[r] = None
        elif name == 'MOVE':
            op[2], op[3] = current[r], target
            current[r] = target
            if op[2] == op[3]:
                continue
        else:
            if current[r] != target or (name != 'SCALE' and current[op[2]] != target):
                raise ValueError('Emitted contact has different physical frames')
            op[-1] = target
        program.append(op)
    word['program'] = program
    strip = lambda w: [e[:-1] if e[0] not in ('INIT', 'CLEAR') else e[:2]
                       for e in w['program'] if e[0] != 'MOVE']
    if strip(source) != strip(word):
        raise ValueError('Logical scalar program changed')
    return word


def native_delta(source, word):
    compiler = load(HERE/'sources/clifford_exact.py', 'candidate_clifford')
    check = load(HERE/'sources/native_verify.py', 'candidate_gaussian')
    def edges(w):
        main = set(w['live']) | {x+112 for x in w['live']}
        pairs = {(tuple(w['frames'][e[2]]), tuple(w['frames'][e[3]]))
                 for e in w['program'] if e[0] == 'MOVE'}
        pairs |= {(tuple(w['frames'][w['initial'][r]]), tuple(w['frames'][w['final'][r]]))
                  for r in range(w['public']) if r not in main}
        return pairs
    added = sorted(edges(word)-edges(source))
    recipes = []
    for a, b in added:
        operator = compiler.canonical_frame(9, a).adjoint().then(compiler.canonical_frame(9, b))
        recipe, _ = compiler.lower_child(operator)
        check.equal(check.recipe(recipe), check.adjoint(check.frame(a, 9))+check.frame(b, 9),
                    'Changed native map and absolute phase')
        if recipe['rank'] != abs(len(a)-len(b)):
            raise ValueError('Changed native rank differs')
        recipes.append([a, b, recipe])
    return {'changed_maps_checked': len(added),
            'recipe_sha256': hashlib.sha256(json.dumps(recipes, sort_keys=True).encode()).hexdigest()}


def run():
    for pin in json.loads((HERE/'SOURCE-PINS.json').read_text()):
        if hashlib.sha256((HERE/pin['local_path']).read_bytes()).hexdigest() != pin['sha256']:
            raise ValueError('Source authentication failed')
    source = json.loads(gzip.decompress((HERE/'sources/selected-word.json.gz').read_bytes()))
    graph = json.loads((HERE/'p104-contact-graph.json').read_text())
    candidate = json.loads((HERE/'monotone-free-helpers.json').read_text())
    word = emit(source, graph, candidate)
    local = load(HERE/'sources/selected-local.py', 'candidate_source_local')
    report = local.replay(word)
    for key in ('paid', 'helpers'):
        if report[key] != {int(r): count for r, count in candidate[key].items()}:
            raise ValueError('Emitted full invoice differs from the model')
    report['transpose'] = local.reverse_columns(word)
    report['native_delta'] = native_delta(source, word)
    raw = (json.dumps(word, sort_keys=True, separators=(',', ':'))+'\n').encode()
    report['word_sha256'] = hashlib.sha256(raw).hexdigest()
    report['status'] = 'KNOWN_PR5601_HELPER307_POSITIVE_CONTROL; NO_NEW_KAPPA'
    (HERE/'FRAME-CANDIDATE.json.gz').write_bytes(gzip.compress(raw, mtime=0))
    (HERE/'FRAME-CANDIDATE-REPLAY.json').write_text(json.dumps(report, indent=2, sort_keys=True)+'\n')
    print(json.dumps(report, indent=2, sort_keys=True))


if __name__ == '__main__':
    run()
