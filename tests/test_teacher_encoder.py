import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch
import anndata as ad
import numpy as np
import pandas as pd
from bremap import TeacherEncoder


class TeacherTests(unittest.TestCase):
    def test_cache_resume_and_input_change(self):
        a=ad.AnnData(np.arange(24).reshape(6,4)+1,
            obs=pd.DataFrame(index=[f'c{i}' for i in range(6)]),
            var=pd.DataFrame(index=['A','B','C','D']))
        with tempfile.TemporaryDirectory() as d:
            p=Path(d); (p/'best_model.pt').write_bytes(b'test');(p/'args.json').write_text('{}');(p/'vocab.json').write_text('{}')
            e=TeacherEncoder('scgpt',p,device='cpu',chunk_size=2)
            calls=[]
            def encode(x, seed):
                calls.append(seed)
                if len(calls)==2: raise RuntimeError('interrupted')
                return np.asarray(x.X)[:,:2].astype(float)
            with patch.object(e,'_encode_chunk',side_effect=encode):
                with self.assertRaisesRegex(RuntimeError,'interrupted'): e.transform(a,cache_dir=p/'cache')
            with patch.object(e,'_encode_chunk',side_effect=lambda x,s: np.asarray(x.X)[:,:2].astype(float)) as f:
                result=e.transform(a,cache_dir=p/'cache');self.assertEqual(f.call_count,2)
                self.assertEqual(result.index.tolist(),a.obs_names.tolist())
                np.testing.assert_array_equal(result,e.transform(a,cache_dir=p/'cache'));self.assertEqual(f.call_count,2)
                b=a.copy();b.X[0,0]+=1
                with self.assertRaisesRegex(ValueError,'fingerprint'): e.transform(b,cache_dir=p/'cache')
            chunk=p/'cache'/'chunk_000000.npz'
            chunk.write_bytes(b'corrupted')
            with self.assertRaisesRegex(ValueError,'cache|Cache'): e.transform(a,cache_dir=p/'cache')

    def test_invalid_options(self):
        with self.assertRaisesRegex(ValueError,'model'): TeacherEncoder('other','.')
        with self.assertRaisesRegex(ValueError,'device'): TeacherEncoder('scgpt','.',device='auto')
        with self.assertRaisesRegex(ValueError,'model_version'): TeacherEncoder('geneformer','.',device='cuda')
        with self.assertRaisesRegex(ValueError,'2048'):
            TeacherEncoder('geneformer','.',model_version='V1',max_length=128)

    def test_attach_has_verified_shared_encoder_identity(self):
        a=ad.AnnData(np.ones((2,3)),obs=pd.DataFrame(index=['a','b']),var=pd.DataFrame(index=['A','B','C']))
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)
            for name in ['args.json','vocab.json','best_model.pt']: (p/name).write_text('{}')
            e=TeacherEncoder('scgpt',p,device='cpu')
            with patch.object(e,'_encode_chunk',return_value=np.ones((2,4))):
                result=e.attach(a,cache_dir=p/'cache',key='X_teacher')
                self.assertEqual(result.obsm['X_teacher'].shape,(2,4))
                self.assertNotIn('X_teacher',a.obsm)

    def test_changed_checkpoint_cannot_reuse_warm_model(self):
        a=ad.AnnData(np.ones((2,3)),obs=pd.DataFrame(index=['a','b']),var=pd.DataFrame(index=['A','B','C']))
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)
            for name in ['args.json','vocab.json','best_model.pt']: (p/name).write_text('{}')
            e=TeacherEncoder('scgpt',p,device='cpu')
            with patch.object(e,'_encode_chunk',return_value=np.ones((2,4))): e.transform(a,cache_dir=p/'first')
            e._bundle='stale weights';(p/'best_model.pt').write_text('changed weights')
            def check(data,seed):
                self.assertIsNone(e._bundle)
                return np.ones((2,4))
            with patch.object(e,'_encode_chunk',side_effect=check): e.transform(a,cache_dir=p/'second')


if __name__=='__main__': unittest.main()
