"""Authenticate all predecessor files before loading retained code.

Ali Romman with substantial OpenAI Codex assistance. Apache-2.0.
"""
from contextlib import contextmanager
from hashlib import sha256
from pathlib import Path
import importlib.util,json,sys
HERE=Path(__file__).resolve().parent

def need(ok,message):
    if not ok:raise ValueError(message)

def load(path,name):
    spec=importlib.util.spec_from_file_location(name,path);module=importlib.util.module_from_spec(spec)
    sys.modules[name]=module;spec.loader.exec_module(module);return module

def authenticate(root=HERE):
    entries=json.loads((root/'INPUTS.json').read_bytes())
    for item in entries:
        need(sha256((root/item['local']).read_bytes()).hexdigest()==item['sha256'],'Frozen input '+item['local'])
    return len(entries)

@contextmanager
def sources():
    authenticate()
    parent=HERE/'parent';outer=load(parent/'sources.py','composition_outer_sources')
    with outer.sources() as(contact,outer_root):
        inner=load(contact/'sources.py','composition_inner_sources')
        with inner.sources() as(base,root):
            yield parent,contact,base,root
