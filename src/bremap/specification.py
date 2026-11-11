"""Versioned tissue declarations; these validate identity, not biological origin."""
from dataclasses import dataclass
import json
from pathlib import Path
import re

from .core import _identifiers


@dataclass(frozen=True)
class ReferenceSpec:
    data: dict

    @classmethod
    def load(cls, path):
        return cls.from_dict(json.loads(Path(path).read_text()))

    @classmethod
    def from_dict(cls, data):
        if not isinstance(data, dict):
            raise ValueError('Reference specification must be a JSON object')
        data = json.loads(json.dumps(data))
        if data.get('schema') != 'bremap-spec-v1':
            raise ValueError('Unsupported reference specification schema')
        for key in ('spec_id', 'spec_version', 'tissue', 'species', 'gene_namespace', 'donor_column', 'label_column'):
            if not isinstance(data.get(key), str) or not data[key].strip():
                raise ValueError(f'Spec requires a nonempty {key}')
        if data['species'] != 'human' or data['gene_namespace'] != 'gene_symbol':
            raise ValueError('Only human gene_symbol specifications are supported')
        states = data.get('states')
        if not isinstance(states, list) or len(states) < 2:
            raise ValueError('Spec requires at least two states')
        states = _identifiers(states, 'states', len(states)).tolist()
        if len(set(states)) != len(states) or any(s.lower() == 'unknown' for s in states):
            raise ValueError('States must be unique; unknown is reserved for abstention')
        data['states'] = sorted(states)
        if data.get('study_column') is not None and (not isinstance(data['study_column'], str) or not data['study_column'].strip()):
            raise ValueError('study_column must be a nonempty name when supplied')
        reference = data.get('reference')
        if not isinstance(reference, dict) or any(not isinstance(reference.get(k), str) or not reference[k].strip()
                                                   for k in ('source', 'eligibility')):
            raise ValueError('Spec reference requires source and eligibility descriptions')
        if data.get('cost_mode') not in ('biology', 'embedding'):
            raise ValueError('Spec cost_mode must be biology or embedding')
        if data.get('normalization') != 'CP10k_log1p_on_reference_gene_universe':
            raise ValueError('Unsupported normalization/feature-universe policy')
        representation = data.get('representation', {})
        if not isinstance(representation, dict) or representation.get('kind') not in ('raw_svd', 'frozen_embedding'):
            raise ValueError('Unsupported representation kind')
        if representation['kind'] == 'raw_svd' and 'key' in representation:
            raise ValueError('raw_svd representation cannot declare an embedding key')
        if representation['kind'] == 'frozen_embedding' and (not isinstance(representation.get('key'), str) or not representation['key'].strip()):
            raise ValueError('Frozen embedding representation requires a key')
        programs = data.get('programs', {})
        anchors = programs.get('anchors') if isinstance(programs, dict) else None
        if not isinstance(anchors, dict) or not anchors:
            raise ValueError('Spec requires nonempty anchor programs')
        heldout = programs.get('heldout_validation', {})
        if not isinstance(heldout, dict):
            raise ValueError('heldout_validation must map program names to genes')
        normalized = []
        for collection in (anchors, heldout):
            result = {}
            for name, genes in collection.items():
                if not isinstance(name, str) or not name.strip() or not isinstance(genes, list) or not genes:
                    raise ValueError('Programs require nonempty names and gene lists')
                validate_gene_symbols(genes)
                result[name] = sorted(g.upper() for g in _identifiers(genes, 'program genes', len(genes)))
            normalized.append(dict(sorted(result.items())))
        if set(g for gs in normalized[0].values() for g in gs) & set(g for gs in normalized[1].values() for g in gs):
            raise ValueError('Anchor and held-out validation genes must be disjoint')
        minimum = programs.get('minimum_genes_per_anchor', 1)
        if isinstance(minimum, bool) or not isinstance(minimum, int) or minimum < 1:
            raise ValueError('minimum_genes_per_anchor must be a positive integer')
        data['programs'] = dict(anchors=normalized[0], heldout_validation=normalized[1], minimum_genes_per_anchor=minimum)
        return cls(data)

    def validate_identity(self, data):
        identity = data.uns.get('bremap')
        if not isinstance(identity, dict):
            raise ValueError("H5AD requires uns['bremap'] tissue/spec identity")
        for key in ('tissue', 'spec_id', 'spec_version'):
            if identity.get(key) != self.data[key]:
                raise ValueError(f'H5AD {key} does not match the reference specification')
        validate_gene_symbols(data.var_names)


def validate_gene_symbols(genes):
    """Validate human-symbol syntax and case collisions, without a gene database."""
    values = _identifiers(genes, 'gene symbols', len(genes))
    if len(set(g.upper() for g in values)) != len(values):
        raise ValueError('Case-insensitive gene symbols must be unique')
    if any(not re.fullmatch(r'[A-Za-z][A-Za-z0-9_.-]*', g) or re.fullmatch(r'ENSG\d+(?:\.\d+)?', g, re.I)
           for g in values):
        raise ValueError('Expected human gene symbols; unsupported gene identifiers')
