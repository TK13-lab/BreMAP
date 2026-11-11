"""Summarize existing study readouts; this does not fit downstream predictors."""
import argparse
import pandas as pd


def summarize(case, path):
    data = pd.read_csv(path)
    if case == 'cd4':
        required = ['donor', 'variant', 'n_query', 'acceptance', 'protein_rho']
        assert set(required) <= set(data)
        assert not data.duplicated(['donor', 'variant']).any()
        assert data.acceptance.between(0, 1).all()
        assert data.protein_rho.dropna().between(-1, 1).all()
        print(data.groupby('variant').agg(
            donors=('donor', 'nunique'), cells=('n_query', 'sum'),
            acceptance=('acceptance', 'mean'), rho=('protein_rho', 'mean')
        ).to_string(float_format=lambda x: f'{x:.4f}'))
    else:
        keys = ['contrast', 'endpoint', 'model', 'patient']
        assert set(keys + ['observed', 'prediction']) <= set(data)
        assert not data.duplicated(keys).any()
        assert data[['observed', 'prediction']].notna().all().all()
        data['squared_error'] = (data.observed - data.prediction)**2
        table = data.groupby(keys[:3]).agg(n=('patient', 'nunique'), mse=('squared_error', 'mean'))
        # Show the named primary contrast, without selecting on achieved accuracy.
        table = table.loc[('d0_d14', 'estrogen_mean')]
        print(table.to_string(float_format=lambda x: f'{x:.6f}'))
    print('PASS: readout keys unique and required values present')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('case', choices=['cd4', 'feline'])
    parser.add_argument('table')
    args = parser.parse_args()
    summarize(args.case, args.table)
