"""Processed-only inputs must reach useful analysis without fabricated counts."""
import tempfile
import unittest
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse
import bremap as br


def specification(data):
    return br.create_spec(spec_id='processed-test', tissue='immune',
        states=['synthetic_quiescent','synthetic_cytotoxic'],
        programs={'quiescent':['CCR7','TCF7','LEF1'],'cytotoxic':['NKG7','PRF1','GZMB']},
        reference_source='Synthetic test', eligibility='Synthetic only',
        minimum_genes_per_anchor=2, expression=data)


class ProcessedLoadingTests(unittest.TestCase):
    def example(self):
        return ad.AnnData(sparse.csr_matrix([[1., 3.], [2., 0.]]),
            obs=pd.DataFrame({'donor_id': ['a', 'b']}, index=['c1', 'c2']),
            var=pd.DataFrame(index=['CD3D', 'MS4A1']))

    def test_normalized_only_is_explicit_and_source_is_unchanged(self):
        source = self.example(); source.X.data *= .25
        with self.assertRaisesRegex(ValueError, 'integer'):
            br.read_data(source)
        data = br.read_data(source, expression_type='library_normalized')
        np.testing.assert_array_equal(data.X.toarray(), [[.25, .75], [.5, 0]])
        self.assertNotIn('bremap_expression', source.uns)
        self.assertEqual(data.uns['bremap_expression']['type'], 'library_normalized')

    def test_raw_selection_uses_raw_genes_and_preserves_cell_embeddings(self):
        source = self.example(); source.obsm['X_pca'] = np.array([[1., 2.], [3., 4.]])
        source.raw = source.copy(); source = source[:, ['CD3D']].copy()
        source.X.data = np.log1p(source.X.data)
        data = br.read_data(source, matrix='raw.X', expression_type='counts')
        self.assertEqual(data.var_names.tolist(), ['CD3D', 'MS4A1'])
        np.testing.assert_array_equal(data.X.toarray(), [[1, 3], [2, 0]])
        np.testing.assert_array_equal(data.obsm['X_pca'], source.obsm['X_pca'])
        self.assertEqual(source.n_vars, 1)

    def test_named_layer_and_old_counts_layer_agree(self):
        source = self.example(); source.layers['counts'] = source.X.copy()
        source.X.data = np.log1p(source.X.data)
        a = br.read_data(source, matrix='layers/counts', expression_type='counts')
        b = br.read_data(source, counts_layer='counts')
        np.testing.assert_array_equal(a.X.toarray(), b.X.toarray())
        with self.assertRaisesRegex(ValueError, 'matrix|counts_layer'):
            br.read_data(source, matrix='X', counts_layer='counts')
        with self.assertRaisesRegex(ValueError, 'missing|Missing'):
            br.read_data(source, matrix='layers/absent')

    def test_log_base_and_transform_errors_are_actionable(self):
        source = self.example(); source.X.data = np.log2(1+source.X.data)
        data = br.read_data(source, expression_type='log1p_library_normalized', log_base=2)
        self.assertEqual(data.uns['bremap_expression']['log_base'], 2.)
        for base in [0, 1, -2, True, float('inf'), 'unknown']:
            with self.assertRaisesRegex(ValueError, 'log_base'):
                br.read_data(source, expression_type='log1p_library_normalized', log_base=base)
        source.X.data[0] = 10000
        with self.assertRaisesRegex(ValueError, 'overflow|Overflow|range'):
            br.read_data(source, expression_type='log1p_library_normalized')

    def test_negative_processed_values_require_processing_identity(self):
        source = self.example(); source.X[0, 1] = -1
        with self.assertRaisesRegex(ValueError, 'processing_id'):
            br.read_data(source, expression_type='preprocessed')
        data = br.read_data(source, expression_type='preprocessed',
            processing_id='training-scaler-v1', units='gene z scores')
        np.testing.assert_array_equal(data.X.toarray(), [[1, -1], [2, 0]])
        self.assertEqual(data.uns['bremap_expression']['units'], 'gene z scores')
        with self.assertRaisesRegex(ValueError, 'nonnegative'):
            br.read_data(source, expression_type='library_normalized')

    def test_inventory_does_not_declare_integer_values_to_be_raw_counts(self):
        source = self.example(); source.layers['logcounts'] = source.X.copy()
        source.raw = source.copy(); source.obsm['X_umap'] = np.zeros((2, 2))
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp)/'processed.h5ad'; source.write_h5ad(path)
            report = br.inspect_data(path)
        self.assertEqual({r['matrix'] for r in report['matrices']}, {'X','raw.X','layers/logcounts'})
        self.assertIn('donor_id', report['obs_columns'])
        self.assertIn('X_umap', report['embeddings'])
        self.assertIn('declaration', report['notes'].lower())
        self.assertTrue(all('integer_like_sample' in row for row in report['matrices']))

    def test_inventory_exposes_ids_storage_and_conditional_operations(self):
        report = br.inspect_data(self.example())
        self.assertTrue(report['matrices'][0]['sparse'])
        self.assertEqual(report['matrices'][0]['cell_ids_preview'], ['c1','c2'])
        self.assertEqual(report['matrices'][0]['gene_ids_preview'], ['CD3D','MS4A1'])
        self.assertIn('donor_id',report['metadata_candidates']['donor'])
        self.assertFalse(report['operations_by_expression']['preprocessed']['teacher'])
        self.assertTrue(report['operations_by_expression']['preprocessed']['mapping'])

    def test_loader_keeps_strict_count_integrality_and_explicit_overrides(self):
        source = self.example(); source.X.data[0] += 1e-7
        with self.assertRaisesRegex(ValueError, 'integer'):
            br.read_data(source)
        data = br.read_data(source,expression_type='log1p_library_normalized',log_base=2)
        with self.assertRaisesRegex(ValueError,'declaration|expression_type'):
            br.read_data(data,log_base=10)
        with self.assertRaisesRegex(ValueError,'declaration|expression_type'):
            br.read_data(data,processing_id='ignored')

    def test_empty_names_and_empty_matrix_selection_are_errors(self):
        for selection in ['', '  ']:
            with self.assertRaisesRegex(ValueError,'matrix'):
                br.read_data(self.example(),matrix=selection)
        source=self.example();source.obs_names=['','c2']
        with self.assertRaisesRegex(ValueError,'cell|identifier'):
            br.read_data(source)
        source=self.example();source.var_names=['CD3D','']
        with self.assertRaisesRegex(ValueError,'gene|identifier'):
            br.read_data(source)


class ProcessedMappingTests(unittest.TestCase):
    def test_count_normalized_and_log_inputs_have_matching_mapping(self):
        source = br.load_example()
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp); reference = {}
            for index, kind in enumerate(['counts','library_normalized','log1p_library_normalized']):
                data = source.copy()
                if kind != 'counts':
                    data.X = data.X / data.X.sum(axis=1, keepdims=True) * .01
                    if kind.startswith('log1p'): data.X = np.log2(1+data.X)
                data = br.read_data(data, expression_type=kind, **({'log_base':2} if kind.startswith('log1p') else {}))
                config = specification(data)
                inputs = br.prepare_partitions(data, config, root/f'inputs{index}')
                for method in ['affinity','transport']:
                    result = br.run_analysis(inputs, config, root/f'{index}_{method}', method=method)
                    if index == 0: reference[method] = result
                    else:
                        pd.testing.assert_frame_equal(reference[method].table(),result.table(),rtol=1e-8,atol=1e-8)
                        pd.testing.assert_frame_equal(reference[method].table('programs'),result.table('programs'),rtol=1e-8,atol=1e-8)
                        replay = br.apply_reference(reference[method].path/'reference.zip', inputs['query'],root/f'replay{index}_{method}')
                        pd.testing.assert_frame_equal(result.table(),replay.table(),rtol=1e-8,atol=1e-8)

    def test_residual_workflow_scores_actual_units_and_rejects_wrong_processing(self):
        source = br.load_example(); source.X = source.X*.1-1
        source.X[0] = 0; source.X[0,0] = 2; source.X[0,1] = -2
        data = br.read_data(source, expression_type='preprocessed',
            processing_id='fixed-development-scaler', units='gene z score')
        config = specification(data)
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp); inputs = br.prepare_partitions(data, config, root/'inputs')
            for method in ['affinity','transport']:
                result = br.run_analysis(inputs,config,root/method,method=method)
                scored = result.table('programs').set_index('cell_id')
                indexes = data.obs_names.get_indexer(scored.index)
                genes = data.var_names.get_indexer(['CCR7','TCF7','LEF1'])
                np.testing.assert_allclose(scored.quiescent, source.X[indexes][:,genes].mean(axis=1),atol=1e-12)
                replay = br.apply_reference(result.path/'reference.zip', inputs['query'],root/(method+'_replay'))
                pd.testing.assert_frame_equal(result.table(),replay.table(),check_exact=True)
            import matplotlib
            matplotlib.use('Agg')
            import matplotlib.pyplot as plt
            fig = result.plot('programs', x='quiescent',y='cytotoxic')
            self.assertIn('gene z score',fig.axes[0].get_xlabel())
            self.assertNotIn('CP10k',fig.axes[0].get_xlabel()); plt.close(fig)
            query = ad.read_h5ad(inputs['query'])
            query.uns['bremap_expression']['processing_id'] = 'separately-fitted-scaler'
            query.write_h5ad(root/'wrong.h5ad')
            with self.assertRaisesRegex(ValueError,'processing|match'):
                br.apply_reference(result.path/'reference.zip',root/'wrong.h5ad',root/'wrong_output')
            self.assertFalse((root/'wrong_output').exists())

    def test_reference_gene_alignment_normalizes_the_selected_universe(self):
        from bremap.expression import features
        data = ad.AnnData(sparse.csr_matrix([[.001,.003,.02]]),var=pd.DataFrame(index=['A','B','C']))
        data = br.read_data(data,expression_type='library_normalized')
        np.testing.assert_allclose(features(data, [0,1]).toarray(),[[np.log1p(2500),np.log1p(7500)]])


class ProcessedStudentTests(unittest.TestCase):
    def inputs(self):
        rng = np.random.default_rng(31)
        data = ad.AnnData(sparse.csr_matrix(rng.poisson(4,(24,12))),
            obs=pd.DataFrame({'donor_id':np.repeat(['train1','train2','validation'],8)},index=[f'c{i}' for i in range(24)]),
            var=pd.DataFrame(index=[f'G{i}' for i in range(12)]))
        targets = pd.DataFrame(rng.normal(size=(24,4)),index=data.obs_names)
        return data, targets

    def fit(self,data,targets):
        return br.StudentEncoder.fit(data[:16],targets.iloc[:16],validation=data[16:],
            validation_teacher=targets.iloc[16:], n_genes=10,n_components=3,
            hidden_size=5,epochs=2,seed=7,teacher_provenance={'model':'custom','checkpoint':'test'})

    def test_student_accepts_equivalent_normalized_and_log_query(self):
        data, targets = self.inputs(); student = self.fit(data,targets)
        expected = student.transform(data)
        for kind in ['library_normalized','log1p_library_normalized']:
            processed = data.copy(); x = data.X.toarray().astype(float)
            processed.X = x / x.sum(axis=1,keepdims=True) * .01
            if kind.startswith('log1p'): processed.X = np.log1p(processed.X)
            processed = br.read_data(processed,expression_type=kind)
            np.testing.assert_allclose(student.transform(processed),expected,atol=2e-5,rtol=2e-5)
            trained = self.fit(processed,targets)
            self.assertTrue(np.isfinite(trained.transform(processed)).all().all())

    def test_processed_student_save_load_and_wrong_contract(self):
        data, targets = self.inputs(); data.X = data.X.toarray().astype(float)-4
        data = br.read_data(data,expression_type='preprocessed',processing_id='fixed-scaler',units='z score')
        student = self.fit(data,targets)
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp)/'student.npz'; student.save(path)
            loaded = br.StudentEncoder.load(path)
            pd.testing.assert_frame_equal(loaded.transform(data),student.transform(data),check_exact=True)
        wrong = data.copy(); wrong.uns['bremap_expression']['processing_id']='other-scaler'
        with self.assertRaisesRegex(ValueError,'processing|match'): student.transform(wrong)
        with self.assertRaisesRegex(ValueError,'processing|match'):
            br.StudentEncoder.fit(data[:16],targets.iloc[:16],validation=wrong[16:],
                validation_teacher=targets.iloc[16:],epochs=1,
                teacher_provenance={'model':'custom','checkpoint':'test'})

    def test_teacher_rejects_residuals_before_download_or_checkpoint_access(self):
        data, _ = self.inputs(); data.X = data.X.toarray().astype(float)-4
        data = br.read_data(data,expression_type='preprocessed',processing_id='residual-v1',units='residual')
        with tempfile.TemporaryDirectory() as tmp:
            cache = Path(tmp)/'unused'
            with self.assertRaisesRegex(ValueError,'preprocessed|residual'):
                br.TeacherEncoder.from_pretrained('geneformer-v1-10m',data=data,cache_dir=cache,confirm=True)
            self.assertFalse(cache.exists())
            teacher = br.TeacherEncoder('scgpt',Path(tmp)/'absent-weights',device='cpu')
            with self.assertRaisesRegex(ValueError,'preprocessed|residual'):
                teacher.transform(data,cache_dir=cache)
            self.assertFalse(cache.exists())


if __name__ == '__main__':
    unittest.main()
