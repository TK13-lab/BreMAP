import tempfile
import unittest
from pathlib import Path
import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse

from bremap import StudentEncoder, evaluate_student


class StudentTests(unittest.TestCase):
    def setUp(self):
        rng = np.random.default_rng(9)
        self.a = ad.AnnData(sparse.csr_matrix(rng.poisson(4, (24, 12))),
            obs=pd.DataFrame({'donor_id': np.repeat(['a','b','c'],8)},
                             index=[f'c{i}' for i in range(24)]),
            var=pd.DataFrame(index=[f'G{i}' for i in range(12)]))
        self.y = pd.DataFrame(rng.normal(size=(24, 4)), index=self.a.obs_names)
        self.kw = dict(n_genes=10, n_components=3, hidden_size=5, epochs=3,
                       teacher_provenance={'model':'custom','checkpoint':'test-v1'}, seed=3)

    def fit(self, **kwargs):
        return StudentEncoder.fit(self.a[:16], self.y.iloc[:16], validation=self.a[16:],
            validation_teacher=self.y.iloc[16:], **(self.kw | kwargs))

    def test_portability_exact_ids_and_gene_order(self):
        model = self.fit()
        result = model.transform(self.a[16:])
        self.assertEqual(result.index.tolist(), self.a.obs_names[16:].tolist())
        np.testing.assert_allclose(result, model.transform(self.a[16:, ::-1]), atol=1e-5)
        with tempfile.TemporaryDirectory() as d:
            p = Path(d)/'student.npz'; model.save(p)
            np.testing.assert_allclose(result, StudentEncoder.load(p).transform(self.a[16:]))
            with self.assertRaises(FileExistsError): model.save(p)
        shuffled = StudentEncoder.fit(self.a[:16], self.y.iloc[:16].iloc[::-1],
            validation=self.a[16:], validation_teacher=self.y.iloc[16:], **self.kw)
        np.testing.assert_allclose(result, shuffled.transform(self.a[16:]), atol=1e-6)

    def test_overlap_counts_coverage_and_wrong_targets_fail(self):
        b = self.a[16:].copy(); b.obs['donor_id'] = 'a'
        with self.assertRaisesRegex(ValueError, 'donor'):
            StudentEncoder.fit(self.a[:16], self.y.iloc[:16], validation=b,
                validation_teacher=self.y.iloc[16:], **self.kw)
        y = self.y.iloc[:16].copy(); y.index = ['wrong']*16
        with self.assertRaisesRegex(ValueError, 'cell|ID'):
            StudentEncoder.fit(self.a[:16], y, validation=self.a[16:],
                validation_teacher=self.y.iloc[16:], **self.kw)
        model = self.fit()
        with self.assertRaisesRegex(ValueError, 'coverage'): model.transform(self.a[:, :2])
        bad = self.a.copy(); bad.X.data[0] = -1
        with self.assertRaisesRegex(ValueError, 'negative|nonnegative'): model.transform(bad)
        with self.assertRaisesRegex(ValueError, 'loss'): self.fit(loss='unknown')

    def test_relational_loss_runs(self):
        model = self.fit(loss='mse_relational', relational_weight=.05)
        self.assertTrue(np.isfinite(model.transform(self.a)).all().all())

    def test_attach_and_tampered_artifact(self):
        model=self.fit();data=model.attach(self.a,key='X_student')
        self.assertEqual(data.obsm['X_student'].shape,(24,4))
        self.assertNotIn('X_student',self.a.obsm)
        with self.assertRaisesRegex(ValueError,'already exists'): model.attach(data,key='X_student')
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/'model.npz';model.save(p)
            with np.load(p,allow_pickle=False) as a: values={k:a[k] for k in a.files}
            values['w0'][0,0]+=1
            np.savez_compressed(p,**values)
            with self.assertRaisesRegex(ValueError,'integrity'): StudentEncoder.load(p)

    def test_incompatible_teacher_metadata_rejected(self):
        y=self.y.iloc[:16].copy();v=self.y.iloc[16:].copy()
        y.attrs['bremap_encoder_metadata']={'model':'custom','checkpoint':'test-v1','preprocessing':'A'}
        v.attrs['bremap_encoder_metadata']={'model':'custom','checkpoint':'other','preprocessing':'A'}
        with self.assertRaisesRegex(ValueError,'provenance|checkpoint'):
            StudentEncoder.fit(self.a[:16],y,validation=self.a[16:],validation_teacher=v,**self.kw)

    def test_target_memory_layout_does_not_change_training(self):
        results=[]
        for order in ['C','F']:
            y=pd.DataFrame(np.array(self.y.to_numpy(),order=order),index=self.y.index)
            model=StudentEncoder.fit(self.a[:16],y.iloc[:16],validation=self.a[16:],
                validation_teacher=y.iloc[16:],**self.kw)
            results.append(model.transform(self.a).to_numpy())
        np.testing.assert_array_equal(results[0],results[1])

    def test_policy_does_not_hide_missing_task_validation(self):
        t = pd.DataFrame({'candidate_state':['a']*20,'accepted':[True]*20,
                          'donor_id':['p']*20}, index=[f'x{i}' for i in range(20)])
        s = t.copy(); s.loc[s.index[:3], 'candidate_state'] = 'b'
        no_truth = evaluate_student(t, s, states=['a','b'], agreement_threshold=.85)
        self.assertTrue(no_truth['agreement_pass'])
        self.assertIsNone(no_truth['joint_pass'])
        truth = t.candidate_state
        r = evaluate_student(t, s, states=['a','b'], truth=truth,
                             agreement_threshold=.90)
        self.assertFalse(r['agreement_pass'])
        self.assertFalse(r['joint_pass'])
        with self.assertRaisesRegex(ValueError, 'threshold'):
            evaluate_student(t,s,states=['a','b'],agreement_threshold=85)


class StudentEvaluationTests(unittest.TestCase):
    def setUp(self):
        self.t = pd.DataFrame({'candidate_state':['a','a','a','a'],
            'accepted':[True]*4,'donor_id':['p','p','q','q']},index=['1','2','3','4'])
        self.truth = pd.Series(['a','b','a','b'],index=self.t.index)

    def test_saved_mapping_and_annotation_files_use_exact_ids_and_donors(self):
        with tempfile.TemporaryDirectory() as tmp:
            p=Path(tmp);self.t.to_csv(p/'teacher.csv',index_label='cell_id')
            self.t.iloc[::-1].to_csv(p/'student.csv',index_label='cell_id')
            labels=self.t[['donor_id']].assign(label=self.truth)
            labels.to_csv(p/'truth.csv',index_label='cell_id')
            r=evaluate_student(p/'teacher.csv',p/'student.csv',states=['a','b'],truth=p/'truth.csv')
            self.assertAlmostEqual(r['student_f1'],1/3)
            labels['donor_id']='wrong';labels.to_csv(p/'truth.csv',index_label='cell_id')
            with self.assertRaisesRegex(ValueError,'donor'):
                evaluate_student(p/'teacher.csv',p/'student.csv',states=['a','b'],truth=p/'truth.csv')

    def test_agreement_does_not_hide_poor_absolute_task_quality(self):
        r=evaluate_student(self.t,self.t.copy(),states=['a','b'],truth=self.truth,
                           minimum_task_f1=.75,minimum_coverage=.8)
        self.assertTrue(r['joint_pass'])  # Existing relative-fidelity contract.
        self.assertAlmostEqual(r['student_f1'],1/3)
        self.assertEqual(r['student_coverage'],1.)
        self.assertFalse(r['task_pass'])
        self.assertFalse(r['case_pass'])
        self.assertEqual(r['thresholds']['minimum_task_f1'],.75)

    def test_absolute_coverage_and_missing_truth_are_explicit(self):
        rejected=self.t.copy();rejected['accepted']=False
        r=evaluate_student(rejected,rejected,states=['a','b'],truth=self.truth,
                           minimum_task_f1=.3,minimum_coverage=.8)
        self.assertTrue(r['joint_pass'])
        self.assertFalse(r['case_pass'])
        r=evaluate_student(self.t,self.t,states=['a','b'],minimum_task_f1=.75)
        self.assertIsNone(r['case_pass'])
        self.assertIsNone(r['student_f1'])

    def test_user_input_errors_are_actionable(self):
        for value in [85,float('nan'),True,'0.85']:
            with self.subTest(value=value),self.assertRaisesRegex(ValueError,'minimum_task_f1'):
                evaluate_student(self.t,self.t,states=['a','b'],minimum_task_f1=value)
        blank=self.t.copy();blank['donor_id']=' '
        with self.assertRaisesRegex(ValueError,'donor'):
            evaluate_student(blank,blank,states=['a','b'])
        with self.assertRaisesRegex(TypeError,'truth.*Series'):
            evaluate_student(self.t,self.t,states=['a','b'],truth=self.truth.to_frame())


if __name__ == '__main__': unittest.main()
