import tempfile
import unittest
from pathlib import Path
import bremap as br


def specification():
    return br.create_spec(spec_id='tutorial', tissue='immune',
        states=['synthetic_quiescent','synthetic_cytotoxic'],
        programs={'quiescent':['CCR7','TCF7','LEF1'],'cytotoxic':['NKG7','PRF1','GZMB']},
        reference_source='Synthetic tutorial',eligibility='Synthetic only',minimum_genes_per_anchor=2)


class UsabilityTests(unittest.TestCase):
    def test_complete_workflow(self):
        with tempfile.TemporaryDirectory() as tmp:
            out=Path(tmp);data=br.load_example()
            inputs=br.prepare_partitions(data,specification(),out/'inputs')
            import anndata as ad
            self.assertNotIn('state',ad.read_h5ad(inputs['query']).obs)
            truth=br.example_annotations(out/'truth.csv')
            result=br.run_analysis(inputs,specification(),out/'run',annotations=truth)
            self.assertEqual(result.summary()['cells'],18)
            import matplotlib
            matplotlib.use('Agg')
            import matplotlib.pyplot as plt
            for kind, options in [('coverage',{}),('affinity',{}),('programs',dict(x='quiescent',y='cytotoxic')),('confusion',dict(annotations=truth))]:
                fig=result.plot(kind,save=out/f'{kind}.png',**options)
                self.assertGreater((out/f'{kind}.png').stat().st_size,1000)
                plt.close(fig)
            labels=br.example_annotations().iloc[:-1]
            with self.assertRaisesRegex(ValueError,'cover all'):
                result.plot('confusion',annotations=labels)
            with self.assertRaises(ValueError):result.plot('programs',x='missing',y='cytotoxic')
            self.assertEqual(br.WorkflowResult(out/'run').summary()['cells'],18)

    def test_reject_overlap_without_writing(self):
        data=br.load_example();data.obs['donor_id']='same'
        with tempfile.TemporaryDirectory() as tmp:
            dest=Path(tmp)/'inputs'
            with self.assertRaisesRegex(ValueError,'overlap'):br.prepare_partitions(data,specification(),dest)
            self.assertFalse(dest.exists())

    def test_reject_counts_and_roles(self):
        data=br.load_example();data.X[0,0]=.5
        with self.assertRaisesRegex(ValueError,'integer'):br.read_data(data)
        data=br.load_example();data.obs.loc[data.obs.role=='query','role']='test'
        with tempfile.TemporaryDirectory() as tmp:
            with self.assertRaisesRegex(ValueError,'exactly'):br.prepare_partitions(data,specification(),Path(tmp)/'inputs')

    def test_bad_calibration_atomic(self):
        data=br.load_example();data.obs.loc[data.obs.role=='calibration','state']='wrong'
        with tempfile.TemporaryDirectory() as tmp:
            dest=Path(tmp)/'inputs'
            with self.assertRaises(ValueError):br.prepare_partitions(data,specification(),dest)
            self.assertFalse(dest.exists())


if __name__=='__main__':unittest.main()
