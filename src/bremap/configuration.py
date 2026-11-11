"""User-defined programs and input diagnostics; no biological validation implied."""
import json
from pathlib import Path

import anndata as ad
import numpy as np

from .core import _identifiers
from .legacy import _embedding
from .preprocessing import normalize_counts
from .specification import ReferenceSpec, validate_gene_symbols


def _programs(value):
    if not isinstance(value, dict) or not value:
        raise ValueError('Programs must be a nonempty mapping of names to gene lists')
    result = {}
    for name, genes in value.items():
        if not isinstance(name, str) or not name.strip() or not isinstance(genes, list) or not genes:
            raise ValueError('Programs require nonempty names and gene lists')
        validate_gene_symbols(genes)
        result[name] = sorted(g.upper() for g in genes)
    return dict(sorted(result.items()))


def _unique_object(pairs):
    result = {}
    for key, value in pairs:
        if key in result:
            raise ValueError(f'Duplicate program or JSON key: {key}')
        result[key] = value
    return result


def read_gene_programs(path):
    """Read a JSON name -> list mapping or GMT (name, description, genes)."""
    path = Path(path)
    if path.suffix.lower() == '.json':
        return _programs(json.loads(path.read_text(), object_pairs_hook=_unique_object))
    if path.suffix.lower() != '.gmt':
        raise ValueError('Gene programs require .json or .gmt input')
    rows = []
    for number, line in enumerate(path.read_text().splitlines(), 1):
        if not line.strip():
            continue
        fields = line.split('\t')
        if len(fields) < 3:
            raise ValueError(f'GMT line {number} requires name, description and genes')
        rows.append((fields[0], fields[2:]))
    return _programs(_unique_object(rows))


def create_spec(*, spec_id, tissue, states, programs, reference_source, eligibility,
                spec_version='1.0', donor_column='donor_id', label_column='state',
                study_column=None, embedding_key=None, cost_mode='biology',
                minimum_genes_per_anchor=1, heldout_validation=None):
    """Construct a custom human gene-symbol spec for the existing mapper."""
    if isinstance(programs, (str, Path)):
        programs = read_gene_programs(programs)
    return ReferenceSpec.from_dict({
        'schema': 'bremap-spec-v1', 'spec_id': spec_id, 'spec_version': spec_version,
        'tissue': tissue, 'species': 'human', 'gene_namespace': 'gene_symbol',
        'states': states, 'donor_column': donor_column, 'label_column': label_column,
        'study_column': study_column, 'reference': {'source': reference_source, 'eligibility': eligibility},
        'programs': {'anchors': _programs(programs), 'heldout_validation': {} if heldout_validation is None else heldout_validation,
                     'minimum_genes_per_anchor': minimum_genes_per_anchor},
        'representation': ({'kind': 'raw_svd'} if embedding_key is None else
                           {'kind': 'frozen_embedding', 'key': embedding_key}),
        'cost_mode': cost_mode, 'normalization': 'CP10k_log1p_on_reference_gene_universe'})


def inspect_input(path, spec, role='query'):
    """Report schema/QC errors separately from descriptive program warnings.

    Scores use the input universe. Mapping may use a smaller frozen reference
    universe; this report does not replace reference/query compatibility checks.
    """
    if role not in ('reference', 'calibration', 'query'):
        raise ValueError('role must be reference, calibration or query')
    spec = spec if isinstance(spec, ReferenceSpec) else ReferenceSpec.load(spec)
    a = ad.read_h5ad(path)
    errors, warnings = [], []
    result = {'schema': 'bremap-input-report-v1', 'role': role,
              'n_cells': a.n_obs, 'n_genes': a.n_vars, 'tissue': spec.data['tissue'],
              'programs': {}, 'overlapping_genes': {}, 'errors': errors, 'warnings': warnings,
              'biological_validation': False, 'score_universe': 'input genes; reference alignment checked during mapping'}
    try:
        spec.validate_identity(a)
        _identifiers(a.obs_names, 'cell IDs', a.n_obs)
        if not a.obs_names.is_unique or not a.n_obs or not a.n_vars:
            raise ValueError('Input requires nonempty, unique cell IDs and genes')
    except (ValueError, TypeError) as exc:
        errors.append(str(exc))
    for column in [spec.data['donor_column']] + ([spec.data['label_column']] if role != 'query' else []):
        if column not in a.obs:
            errors.append(f'Missing metadata column: {column}')
        else:
            try:
                values = _identifiers(a.obs[column], column, a.n_obs)
                result[column + '_counts'] = {str(k): int(v) for k, v in a.obs[column].value_counts().items()}
                if column == spec.data['label_column'] and set(values) != set(spec.data['states']):
                    errors.append('Reference/calibration labels must contain exactly the declared states')
            except ValueError as exc:
                errors.append(str(exc))
    logs = None
    try:
        if a.X is None:
            raise ValueError('Raw integer counts in X are required')
        logs = normalize_counts(a.X)
        zero = int((np.asarray(logs.sum(axis=1)).ravel() == 0).sum())
        result['zero_library_cells'] = zero
        if zero:
            (warnings if role == 'query' else errors).append(f'{zero} zero-library cells; query cells will abstain')
    except (ValueError, TypeError) as exc:
        errors.append(str(exc))
    genes = np.asarray([str(g).upper() for g in a.var_names])
    heldout = {g for values in spec.data['programs']['heldout_validation'].values() for g in values}
    result['validation_genes_present_in_input'] = sorted(heldout & set(genes))
    if result['validation_genes_present_in_input']:
        warnings.append('Validation genes remain in the input universe; program holdout alone does not exclude them from an encoder')
    membership = {}
    for name, symbols in spec.data['programs']['anchors'].items():
        matched = sorted(set(symbols) & set(genes))
        missing = sorted(set(symbols) - set(genes))
        row = {'matched_genes': matched, 'missing_genes': missing,
               'coverage_fraction': len(matched) / len(symbols)}
        result['programs'][name] = row
        for gene in symbols:
            membership.setdefault(gene, []).append(name)
        if len(matched) < spec.data['programs']['minimum_genes_per_anchor']:
            errors.append(f'Insufficient anchor coverage: {name}')
        if missing:
            warnings.append(f'{name}: {len(missing)} requested genes absent')
        if matched and logs is not None and a.n_obs:
            scores = np.asarray(logs[:, np.isin(genes, symbols)].mean(axis=1)).ravel()
            row.update(score_mean=float(scores.mean()), score_sd=float(scores.std()),
                       zero_score_fraction=float((scores == 0).mean()))
            if np.isclose(scores.std(), 0):
                warnings.append(f'{name}: constant program score in this input')
    result['overlapping_genes'] = {g: names for g, names in membership.items() if len(names) > 1}
    if result['overlapping_genes']:
        warnings.append('Programs share genes; their biological costs are not independent')
    rep = spec.data['representation']
    if rep['kind'] == 'frozen_embedding':
        try:
            z, provenance = _embedding(a, rep['key'])
            result['embedding'] = {'key': rep['key'], 'dimensions': z.shape[1], 'provenance': provenance}
        except (ValueError, KeyError) as exc:
            errors.append(str(exc))
    result['valid'] = not errors
    return result
