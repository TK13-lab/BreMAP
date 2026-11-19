"""Read actual native-R fixtures in Python, with no R runtime dependency."""
from pathlib import Path
import unittest
import numpy as np
import pandas as pd
from scipy import sparse
import bremap as br

FIXTURES = Path(__file__).parent/'fixtures'/'processed_rds'


class RDSInputTests(unittest.TestCase):
    def test_seurat_v4_v5_and_sce_exact_values_and_metadata(self):
        for filename in ['seurat_assay.rds','seurat_assay5.rds','sce.rds']:
            with self.subTest(filename=filename):
                data=br.read_data(FIXTURES/filename,assay='RNA' if filename!='sce.rds' else None,
                    matrix='counts',expression_type='counts')
                self.assertTrue(sparse.issparse(data.X))
                self.assertEqual(data.obs_names.tolist(),['cell2','cell1','cell4','cell3'])
                self.assertEqual(data.var_names.tolist(),['CD3D','MS4A1','LYZ'])
                np.testing.assert_array_equal(data.X.toarray(),[[1,0,3],[0,2,2],[4,1,0],[2,2,1]])
                self.assertEqual(data.obs.donor_id.tolist(),['donorB','donorA','donorB','donorC'])
                self.assertEqual(data.obs.state.astype(str).tolist(),['T','B','T','B'])
                self.assertTrue(pd.isna(data.obs.numeric_note.iloc[1]))
                key='PCA' if filename=='sce.rds' else 'pca'
                np.testing.assert_array_equal(data.obsm[key],[[1,5],[2,6],[3,7],[4,8]])

    def test_data_only_assay_and_inventory_need_no_counts(self):
        path=FIXTURES/'seurat_data_only.rds'
        report=br.inspect_data(path)
        self.assertEqual([(r['assay'],r['matrix']) for r in report['matrices']],[('RNA','data')])
        data=br.read_data(path,assay='RNA',matrix='data',expression_type='log1p_library_normalized')
        np.testing.assert_allclose(data.X.toarray(),np.log1p(np.array([[1,0,3],[0,2,2],[4,1,0],[2,2,1]])*.25))
        with self.assertRaisesRegex(ValueError,'counts|available'):
            br.read_data(path,assay='RNA',matrix='counts')

    def test_layer_membership_and_split_selection_are_explicit(self):
        data=br.read_data(FIXTURES/'seurat_subset_layer.rds',assay='RNA',matrix='small',
                         expression_type='log1p_library_normalized')
        self.assertEqual(data.obs_names.tolist(),['cell2','cell1'])
        self.assertEqual(data.var_names.tolist(),['CD3D','MS4A1'])
        np.testing.assert_allclose(data.X.toarray(),np.log1p([[.25,0],[0,.5]]))
        self.assertEqual(len(data.obsm['pca']),2)
        with self.assertRaisesRegex(ValueError,'layer|matrix|available'):
            br.read_data(FIXTURES/'seurat_split.rds',assay='RNA',matrix='counts')
        selected=br.read_data(FIXTURES/'seurat_split.rds',assay='RNA',matrix='counts.batchA')
        self.assertEqual(selected.obs_names.tolist(),['cell1','cell3'])

    def test_named_sparse_matrix_and_unsupported_class(self):
        data=br.read_data(FIXTURES/'sparse_matrix.rds')
        self.assertEqual(data.shape,(4,3));self.assertTrue(sparse.issparse(data.X))
        with self.assertRaisesRegex(ValueError,'class|supported'):
            br.read_data(FIXTURES/'unsupported.rds')

    def test_dense_duplicates_unnamed_and_optional_dependency(self):
        data=br.read_data(FIXTURES/'dense_matrix.rds')
        np.testing.assert_array_equal(data.X,[[1,0,3],[0,2,2],[4,1,0],[2,2,1]])
        for name in ['duplicate_cells.rds','unnamed_matrix.rds','missing_gene.rds']:
            with self.assertRaisesRegex(ValueError,'names|unique'):
                br.read_data(FIXTURES/name)
        with self.assertRaisesRegex(ValueError,'Delayed|external|Unsupported'):
            br.read_data(FIXTURES/'delayed_matrix.rds')
        from unittest.mock import patch
        with patch.dict('sys.modules', {'rdata': None}):
            with self.assertRaisesRegex(ImportError,'bremap\\[rds\\]'):
                br.read_data(FIXTURES/'dense_matrix.rds')

    def test_object_without_reductions_still_loads(self):
        data=br.read_data(FIXTURES/'seurat_no_reductions.rds',assay='RNA',matrix='data',
                         expression_type='log1p_library_normalized')
        self.assertEqual(data.shape,(4,3))
        self.assertEqual(list(data.obsm),[])


if __name__=='__main__':unittest.main()
