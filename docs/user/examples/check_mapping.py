"""Check a mapping directory; print aggregate diagnostics, never patient identifiers."""
import argparse
from pathlib import Path
import numpy as np
import pandas as pd


def check(path):
    frame = pd.read_csv(Path(path) / 'cells.csv')
    assert frame.cell_id.is_unique, 'Duplicate cell IDs'
    assert frame.donor_id.notna().all(), 'Missing donors'
    assert frame.accepted.isin([True, False]).all(), 'Invalid acceptance values'
    columns = [c for c in frame if c.startswith('conditional_affinity::')]
    assert len(columns) >= 2, 'At least two states required'
    scored = frame[columns].notna().all(axis=1)
    values = frame.loc[scored, columns].to_numpy()
    assert np.isfinite(values).all() and (values >= 0).all()
    assert np.allclose(values.sum(axis=1), 1, atol=1e-6), 'Affinities do not sum to one'
    assert scored.any(), 'No scored cells'
    print(f'cells={len(frame)} donors={frame.donor_id.nunique()} states={len(columns)}')
    print(f'scored={scored.sum()} unscored={(~scored).sum()}')
    print(f'acceptance_cell_weighted={frame.accepted.mean():.4f}')
    print(f'acceptance_donor_macro={frame.groupby("donor_id").accepted.mean().mean():.4f}')
    print('candidate_state_counts (all cells; not tissue proportions):')
    print(frame.candidate_state.value_counts().sort_index().to_string())
    print('PASS: unique cells, valid donors, acceptance flags and normalized affinities')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('mapping', type=Path)
    check(parser.parse_args().mapping)
