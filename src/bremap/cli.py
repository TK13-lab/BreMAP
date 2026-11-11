"""Public command dispatcher, sharing the path API and retaining legacy calls."""
import argparse
import json
from pathlib import Path
from zipfile import BadZipFile

from . import __version__, legacy
from .api import fit_reference, calibrate, map_query, assess_reliability, evaluate
from .specification import ReferenceSpec
from .configuration import create_spec, read_gene_programs, inspect_input
from .resources_io import export_reference, import_reference, export_embeddings


def _write_new_json(path, value):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open('x') as handle:
        json.dump(value, handle, indent=2)
        handle.write('\n')
    return path


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--version', action='version', version=__version__)
    sub = parser.add_subparsers(dest='command', required=True)
    fit = sub.add_parser('fit-reference', help='Fit frozen reference preprocessing and prototypes')
    fit.add_argument('--reference', required=True)
    fit.add_argument('--model', required=True)
    fit.add_argument('--spec', help='Versioned tissue specification; requires explicit --policy')
    fit.add_argument('--policy', choices=('normal_donor_v1', 'legacy_cli_v1'))
    fit.add_argument('--method', choices=('affinity', 'transport'))
    fit.add_argument('--programs', help='Legacy programs JSON (when no --spec is supplied)')
    fit.add_argument('--label-column', default='reference_state')
    fit.add_argument('--donor-column', default='global_donor_key')
    fit.add_argument('--embedding-key')
    cal = sub.add_parser('calibrate', help='Calibrate with normal donors; use --out for a separate artifact')
    cal.add_argument('--model', required=True)
    cal.add_argument('--calibration', required=True)
    cal.add_argument('--out', help='Required for immutable spec-based calibration')
    mapping = sub.add_parser('map', help='Map disjoint query donor groups')
    mapping.add_argument('--model', required=True)
    mapping.add_argument('--query', required=True)
    mapping.add_argument('--out', required=True)
    reliability = sub.add_parser('reliability', help='Assess frozen reference donor/study sensitivity')
    reliability.add_argument('--model', required=True)
    reliability.add_argument('--reference', required=True)
    reliability.add_argument('--query', required=True)
    reliability.add_argument('--out', required=True)
    reliability.add_argument('--n-bootstrap', type=int, default=20)
    reliability.add_argument('--seed', type=int, default=42)
    reliability.add_argument('--leave-study-out', action='store_true')
    evaluation = sub.add_parser('evaluate', help='Evaluate frozen mapping with exact cell/donor annotations')
    evaluation.add_argument('--mapping', required=True)
    evaluation.add_argument('--annotations', required=True)
    evaluation.add_argument('--out', required=True)
    validator = sub.add_parser('validate-spec')
    validator.add_argument('--spec', required=True)
    init = sub.add_parser('init-spec', help='Create a custom human tissue/state specification')
    init.add_argument('--spec-id', required=True)
    init.add_argument('--tissue', required=True)
    init.add_argument('--states', nargs='+', required=True)
    init.add_argument('--programs', required=True, help='Named gene lists in JSON or GMT')
    init.add_argument('--heldout-programs', help='Optional disjoint validation programs, JSON/GMT')
    init.add_argument('--reference-source', required=True)
    init.add_argument('--eligibility', required=True)
    init.add_argument('--spec-version', default='1.0')
    init.add_argument('--donor-column', default='donor_id')
    init.add_argument('--label-column', default='state')
    init.add_argument('--study-column')
    init.add_argument('--embedding-key')
    init.add_argument('--minimum-genes-per-anchor', type=int, default=1)
    init.add_argument('--cost-mode', choices=('biology', 'embedding'), default='biology')
    init.add_argument('--out', required=True)
    inspect = sub.add_parser('inspect-input', help='Report gene/program coverage, QC and embedding provenance')
    inspect.add_argument('--input', required=True)
    inspect.add_argument('--spec', required=True)
    inspect.add_argument('--role', choices=('reference', 'calibration', 'query'), default='query')
    inspect.add_argument('--out', required=True)
    export = sub.add_parser('export-reference', help='Export a portable reference archive without counts')
    export.add_argument('--model', required=True)
    export.add_argument('--out', required=True)
    restore = sub.add_parser('import-reference', help='Verify and import a portable reference archive')
    restore.add_argument('--archive', required=True)
    restore.add_argument('--out', required=True)
    vectors = sub.add_parser('export-embeddings', help='Export vectors and explicitly selected metadata, without counts')
    vectors.add_argument('--input', required=True)
    vectors.add_argument('--key', required=True)
    vectors.add_argument('--obs-columns', nargs='*', default=[])
    vectors.add_argument('--out', required=True)
    toy = sub.add_parser('toy')
    toy.add_argument('--out', required=True)
    args = parser.parse_args(argv)
    try:
        if args.command == 'init-spec':
            spec = create_spec(spec_id=args.spec_id, tissue=args.tissue, states=args.states,
                programs=args.programs, reference_source=args.reference_source, eligibility=args.eligibility,
                spec_version=args.spec_version, donor_column=args.donor_column, label_column=args.label_column,
                study_column=args.study_column, embedding_key=args.embedding_key, cost_mode=args.cost_mode,
                minimum_genes_per_anchor=args.minimum_genes_per_anchor,
                heldout_validation=read_gene_programs(args.heldout_programs) if args.heldout_programs else None)
            return _write_new_json(args.out, spec.data)
        if args.command == 'inspect-input':
            report = inspect_input(args.input, args.spec, args.role)
            _write_new_json(args.out, report)
            if not report['valid']:
                raise ValueError('Input inspection failed; see errors in ' + args.out)
            return report
        if args.command == 'export-reference':
            return export_reference(args.model, args.out)
        if args.command == 'import-reference':
            return import_reference(args.archive, args.out)
        if args.command == 'export-embeddings':
            return export_embeddings(args.input, args.out, key=args.key, obs_columns=args.obs_columns)
        if args.command == 'fit-reference':
            if args.spec:
                return fit_reference(args.reference, args.spec, args.model, args.method or 'transport', args.policy)
            if args.policy:
                raise ValueError('--policy requires --spec')
            if not args.programs:
                raise ValueError('Legacy fit requires --programs; tissue-aware fit requires --spec')
            args.method = args.method or 'affinity'
            return legacy.fit_reference(args)
        if args.command == 'calibrate':
            if args.out:
                return calibrate(args.model, args.calibration, args.out)
            if (Path(args.model) / 'artifact.json').exists():
                raise ValueError('Spec-based calibration requires --out for an immutable artifact')
            return legacy.calibrate(args)
        if args.command == 'map':
            return map_query(args.model, args.query, args.out)
        if args.command == 'reliability':
            return assess_reliability(args.model, args.reference, args.query, args.out,
                                      args.n_bootstrap, args.seed, args.leave_study_out)
        if args.command == 'evaluate':
            return evaluate(args.mapping, args.annotations, args.out)
        if args.command == 'toy':
            return legacy.toy(args)
        ReferenceSpec.load(args.spec)
    except (ValueError, OSError, KeyError, TypeError, BadZipFile) as exc:
        parser.error(str(exc))


def console_main():
    """Console launchers call sys.exit; successful API results are not exit codes."""
    main()


if __name__ == '__main__':
    console_main()
