"""Run one prepared biological task with installed BReMAP APIs only."""
import argparse
from pathlib import Path
import anndata as ad
import bremap as br
from bremap.workflow import run_workflow


def run(inputs, output, method):
    spec = br.ReferenceSpec.load(inputs / 'spec.json')
    donors = {}
    for role in ('reference', 'calibration', 'query'):
        path = inputs / f'{role}.h5ad'
        report = br.inspect_input(path, spec, role=role)
        assert report['valid'], (role, report['errors'])
        data = ad.read_h5ad(path, backed='r')
        donors[role] = set(data.obs[spec.data['donor_column']].astype(str))
        print(f'{role}: cells={data.n_obs} genes={data.n_vars} donors={len(donors[role])} valid=True')
        data.file.close()
    assert not donors['reference'] & donors['calibration'], 'Reference/calibration donor overlap'
    assert not donors['reference'] & donors['query'], 'Reference/query donor overlap'
    assert not donors['calibration'] & donors['query'], 'Calibration/query donor overlap'
    recipe = dict(schema='bremap-workflow-v1', name='prepared-case', spec=spec,
                  reference=inputs/'reference.h5ad', calibration=inputs/'calibration.h5ad',
                  query=inputs/'query.h5ad', method=method, policy='normal_donor_v1', output=output)
    if (inputs/'annotations.csv').is_file():
        recipe['annotations'] = inputs/'annotations.csv'
    result = run_workflow(recipe)
    assert (result/'reference.zip').is_file()
    assert (result/'programs/cells.csv').is_file()
    print('PASS: donor-disjoint inputs; mapping, programs and reference archive written')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--inputs', type=Path, required=True)
    parser.add_argument('--out', type=Path, required=True)
    parser.add_argument('--method', choices=['affinity', 'transport'], default='affinity')
    args = parser.parse_args()
    run(args.inputs, args.out, args.method)
