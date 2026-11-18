"""Consent, cache and model-selection contracts; no real weights in unit tests."""
import contextlib
import io
import json
import tempfile
import unittest
import zipfile
from types import SimpleNamespace
from pathlib import Path
from unittest.mock import patch

import bremap as br


class CheckpointTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.cache = Path(self.tmp.name) / 'models'
        self.stdout = contextlib.redirect_stdout(io.StringIO())
        self.stdout.__enter__()
        self.addCleanup(self.stdout.__exit__, None, None, None)

    @staticmethod
    def transfer(entry, filename, destination):
        if filename == 'args.json':
            destination.write_text(json.dumps(dict(embsize=8, nheads=2, d_hid=8,
                nlayers=2, n_layers_cls=1, dropout=0.1, pad_token='<pad>', pad_value=-2)))
        elif filename == 'vocab.json':
            destination.write_text('{"A": 0, "B": 1}')
        else:
            with zipfile.ZipFile(destination, 'w') as archive:
                archive.writestr('model/data.pkl', b'fixture; never deserialized')
                archive.writestr('model/version', b'3')

    def download(self, **kwargs):
        return br.download_checkpoint('scgpt-whole-human', cache_dir=self.cache, **kwargs)

    def test_catalog_is_offline_and_rejects_unknown_family(self):
        table = br.list_checkpoints('geneformer')
        self.assertIn('geneformer-v1-10m', table.checkpoint_id.tolist())
        self.assertEqual(set(table.model_version), {'V1', 'V2'})
        self.assertFalse(self.cache.exists())
        with self.assertRaisesRegex(ValueError, 'scgpt.*geneformer'):
            br.list_checkpoints('unknown')

    def test_refusal_and_invalid_options_have_no_side_effects(self):
        with patch('bremap.checkpoints._download_file', side_effect=AssertionError('network forbidden')):
            for answer in ['', 'n', 'maybe']:
                with patch('builtins.input', return_value=answer):
                    with self.assertRaisesRegex(RuntimeError, 'cancelled'):
                        self.download()
            with self.assertRaisesRegex(RuntimeError, 'cancelled'):
                self.download(confirm=False)
            with self.assertRaisesRegex(ValueError, 'confirm'):
                self.download(confirm='yes')
            with self.assertRaisesRegex(ValueError, 'Unknown checkpoint'):
                br.download_checkpoint('../weights', cache_dir=self.cache, confirm=True)
            with patch('builtins.input', side_effect=EOFError):
                with self.assertRaisesRegex(RuntimeError, 'confirm=True'):
                    self.download()
        self.assertFalse(self.cache.exists())

    def test_interactive_choice_and_confirmed_cache_reuse(self):
        with patch('builtins.input', side_effect=['1', 'yes']), \
             patch('bremap.checkpoints._download_file', side_effect=self.transfer):
            path = br.download_checkpoint(cache_dir=self.cache)
        self.assertEqual(path.name, 'scgpt-whole-human')
        self.assertEqual(json.loads((path/'vocab.json').read_text()), {'A': 0, 'B': 1})
        with patch('builtins.input', side_effect=AssertionError('cached; no prompt')), \
             patch('bremap.checkpoints._download_file', side_effect=AssertionError('cached; no network')):
            self.assertEqual(path, self.download())
        (path/'best_model.pt').write_bytes(b'changed')
        with self.assertRaisesRegex(ValueError, 'corrupt|changed'):
            self.download(confirm=True)

    def test_failed_transfer_resumes_completed_files_only(self):
        def fail(entry, name, dest):
            if name == 'best_model.pt':
                raise OSError('connection lost')
            self.transfer(entry, name, dest)
        with patch('bremap.checkpoints._download_file', side_effect=fail):
            with self.assertRaisesRegex(RuntimeError, 'retry|Retry'):
                self.download(confirm=True)
        self.assertFalse((self.cache/'scgpt-whole-human').exists())
        def resume(entry, name, dest):
            if name != 'best_model.pt':
                raise AssertionError('completed file downloaded twice')
            self.transfer(entry, name, dest)
        with patch('bremap.checkpoints._download_file', side_effect=resume):
            result = self.download(confirm=True)
        self.assertTrue((result/'best_model.pt').is_file())
        self.assertFalse((self.cache/'.scgpt-whole-human.partial').exists())

    def test_html_or_lfs_pointer_is_not_a_completed_checkpoint(self):
        for payload in [b'<!DOCTYPE html><html>quota exceeded</html>',
                        b'version https://git-lfs.github.com/spec/v1\noid sha256:123']:
            with self.subTest(payload=payload):
                def bad(entry, name, dest):
                    self.transfer(entry, name, dest)
                    if name == 'best_model.pt': dest.write_bytes(payload)
                with patch('bremap.checkpoints._download_file', side_effect=bad):
                    with self.assertRaisesRegex(RuntimeError, 'weights|HTML|pointer'):
                        self.download(confirm=True)
                self.assertFalse((self.cache/'scgpt-whole-human').exists())

    def test_unmanaged_directory_is_never_overwritten(self):
        path = self.cache/'scgpt-whole-human'
        path.mkdir(parents=True)
        (path/'note.txt').write_text('keep me')
        with self.assertRaisesRegex(ValueError, 'unmanaged|manifest'):
            self.download(confirm=True)
        self.assertEqual((path/'note.txt').read_text(), 'keep me')

    def test_truncated_weights_are_rejected_and_retry_replaces_them(self):
        def truncated(entry, name, dest):
            self.transfer(entry, name, dest)
            if name == 'best_model.pt': dest.write_bytes(b'PK\x03\x04truncated archive')
        with patch('bremap.checkpoints._download_file', side_effect=truncated):
            with self.assertRaisesRegex(RuntimeError, 'archive|weights'):
                self.download(confirm=True)
        with patch('bremap.checkpoints._download_file', side_effect=self.transfer):
            path = self.download(confirm=True)
        self.assertTrue(zipfile.is_zipfile(path/'best_model.pt'))

    def test_partial_cache_corruption_and_concurrent_writers_fail(self):
        lock = self.cache/'.scgpt-whole-human.lock'
        lock.mkdir(parents=True)
        with self.assertRaisesRegex(RuntimeError, 'already in progress'):
            self.download(confirm=True)
        lock.rmdir()
        def fail(entry, name, dest):
            if name == 'vocab.json': raise OSError('interrupted')
            self.transfer(entry, name, dest)
        with patch('bremap.checkpoints._download_file', side_effect=fail):
            with self.assertRaises(RuntimeError): self.download(confirm=True)
        partial = self.cache/'.scgpt-whole-human.partial'
        (partial/'args.json').write_text('{}')
        with self.assertRaisesRegex(ValueError, 'corrupt|changed'):
            self.download(confirm=True)
        self.assertFalse(lock.exists())

    def test_transports_use_only_selected_files_and_pinned_revision(self):
        from bremap.checkpoints import _download_file, _select
        entry = _select('geneformer-v1-10m')
        def hub(*, repo_id, filename, revision, local_dir, token):
            if (repo_id != 'ctheodoris/Geneformer' or
                    filename != 'Geneformer-V1-10M/config.json' or
                    revision != '1f7fbae4e469a5f4f1af8c111a529cfe1b3829f5' or token is not False):
                raise AssertionError('wrong checkpoint source')
            target = Path(local_dir)/filename
            target.parent.mkdir(parents=True)
            target.write_text('{"max_position_embeddings":2048}')
            return str(target)
        dest = Path(self.tmp.name)/'config.json'
        with patch.dict('sys.modules', {'huggingface_hub': SimpleNamespace(hf_hub_download=hub)}):
            _download_file(entry, 'config.json', dest)
        self.assertEqual(json.loads(dest.read_text())['max_position_embeddings'], 2048)
        def drive(*, id, output, quiet, resume, use_cookies):
            if id != 'wanted': raise AssertionError('unrelated file requested')
            Path(output).write_text('{"gene":0}')
            return output
        listing = [SimpleNamespace(id='wanted', path='vocab.json'),
                   SimpleNamespace(id='unrelated', path='extra_weights.pt')]
        dest = Path(self.tmp.name)/'vocab.json'
        with patch.dict('sys.modules', {'gdown': SimpleNamespace(
                download_folder=lambda **kw: listing, download=drive)}):
            _download_file(_select('scgpt-whole-human'), 'vocab.json', dest)
        self.assertEqual(json.loads(dest.read_text()), {'gene': 0})

    def test_factory_sets_version_and_validates_options_before_download(self):
        with patch('bremap.checkpoints.download_checkpoint', return_value=Path(self.tmp.name)):
            teacher = br.TeacherEncoder.from_pretrained('geneformer-v1-10m', confirm=True)
        self.assertEqual(teacher.model, 'geneformer')
        self.assertEqual(teacher.options['model_version'], 'V1')
        self.assertEqual(teacher.options['max_length'], 2048)
        with patch('bremap.checkpoints.download_checkpoint', side_effect=AssertionError('no download')):
            with self.assertRaisesRegex(ValueError, 'model_version'):
                br.TeacherEncoder.from_pretrained('geneformer-v1-10m', model_version='V2')
            with self.assertRaisesRegex(ValueError, 'device'):
                br.TeacherEncoder.from_pretrained('scgpt-whole-human', device='auto')


if __name__ == '__main__':
    unittest.main()
