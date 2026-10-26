"""Ambiguous coordination JSON must fail before any network mutation."""
import contextlib
import io
import json
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'scripts'))
import research_coordination as board
import post_update as post
sys.path.pop(0)


def value():
    return {'schema': board.SCHEMA, 'kind': 'HANDOFF', 'task': 'fixture',
            'scope': 'one missing check', 'base': 'a' * 40, 'artifact': None}


def fenced(raw):
    return 'A scoped handoff.\n\n```json\n' + raw + '\n```\n'


class CoordinationJSONTests(unittest.TestCase):
    def test_duplicate_fields_and_nonfinite_extensions_are_rejected(self):
        raw = json.dumps(value())
        for extra in ('"kind":"DONE"', '"base":"' + 'b' * 40 + '"',
                      '"extra":{"x":1,"x":2}', '"extra":NaN', '"extra":Infinity'):
            with self.subTest(extra=extra):
                rows, warnings = board.messages(fenced(raw[:-1] + ',' + extra + '}'))
                self.assertEqual(rows, [])
                self.assertEqual(len(warnings), 1)

    def test_prepared_body_rejects_duplicate_status(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / 'body.md'
            path.write_text(fenced(json.dumps(value())[:-1] + ',"kind":"DONE"}'))
            with self.assertRaisesRegex(ValueError, 'duplicate.*kind'):
                post.prepare_body(path)

    def test_metadata_file_rejects_before_post(self):
        with tempfile.TemporaryDirectory() as tmp:
            body, metadata = Path(tmp) / 'body.md', Path(tmp) / 'meta.json'
            body.write_text('A handoff.\n')
            metadata.write_text(json.dumps(value())[:-1] + ',"kind":"DONE"}')
            args = ['post_update.py', '--discussion', '1166', '--body-file', str(body), '--metadata-file', str(metadata)]
            error = io.StringIO()
            with patch.object(sys, 'argv', args), patch.object(post, 'post') as send, contextlib.redirect_stderr(error):
                self.assertEqual(post.main(), 1)
            send.assert_not_called()
            self.assertIn('duplicate coordination JSON key: kind', error.getvalue())

    def test_valid_unicode_metadata_and_prose_round_trip(self):
        data = dict(value(), scope='A κ check', extension={'count': 3})
        raw = json.dumps(data, ensure_ascii=False)
        self.assertEqual(board.parse_message(raw), data)
        self.assertEqual(board.messages(fenced(raw)), ([data], []))


if __name__ == '__main__':
    unittest.main()
