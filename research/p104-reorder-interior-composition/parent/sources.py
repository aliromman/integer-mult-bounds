"""Authenticate the retained PR5660 export before any source import.

Ali Romman with substantial OpenAI Codex assistance. Apache-2.0.
The archive checks follow the repository exporter and flyingviper's PR5613.
"""
from contextlib import contextmanager
from hashlib import sha1, sha256
from pathlib import Path, PurePosixPath
import importlib.util
import io
import json
import sys
import tarfile
import tempfile

HERE = Path(__file__).resolve().parent


def need(ok, why):
    if not ok:
        raise ValueError(why)


def load(path, name):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


@contextmanager
def sources():
    pin = json.loads((HERE/'SOURCES.json').read_bytes())['source_export']
    raw = (HERE/pin['local_path']).read_bytes()
    need(sha256(raw).hexdigest() == pin['sha256'], 'Exact predecessor archive SHA256')
    need(sha1(b'blob '+str(len(raw)).encode()+b'\0'+raw).hexdigest() == pin['git_blob'],
         'Exact retained archive Git blob')
    with tempfile.TemporaryDirectory(prefix='p104-interior-source-') as directory:
        root, names = Path(directory), set()
        with tarfile.open(fileobj=io.BytesIO(raw), mode='r:') as bundle:
            for member in bundle.getmembers():
                path = PurePosixPath(member.name)
                need(member.isfile() and not path.is_absolute() and '..' not in path.parts
                     and '\\' not in member.name and ':' not in member.name,
                     'Regular relative source path')
                need(member.name not in names, 'Unique source path')
                names.add(member.name)
                target = root.joinpath(*path.parts)
                target.parent.mkdir(parents=True, exist_ok=True)
                target.write_bytes(bundle.extractfile(member).read())
        manifest = json.loads((root/'.research-export.json').read_bytes())
        need(manifest['commit'] == pin['commit'] and manifest['package'] == pin['package']
             and len(manifest['files']) == pin['files'], 'Exact source commit and inventory')
        need(names == {'.research-export.json'} | {r['path'] for r in manifest['files']},
             'No missing or additional source payload')
        for row in manifest['files']:
            payload = (root/row['path']).read_bytes()
            need(len(payload) == row['size'] and sha256(payload).hexdigest() == row['sha256'],
                 'Authenticated payload: '+row['path'])
        yield root/pin['package'], root
