"""Declarative orchestration of the existing immutable reference API."""
import json
import hashlib
import os
from pathlib import Path
import shutil
import tempfile
import time

import numpy as np
import pandas as pd

from .api import fit_reference, calibrate, map_query, evaluate, assess_reliability, load_reference, transform_input
from .configuration import inspect_input, _unique_object
from .legacy import _versions, _write_json
from .preprocessing import file_hash
from .resources_io import import_reference, export_reference
from .specification import ReferenceSpec


def _recipe(recipe, base_dir=None):
    if isinstance(recipe, dict):
        value = dict(recipe)
        if isinstance(value.get('spec'), ReferenceSpec):
            value['spec'] = value['spec'].data
        for key in ['reference', 'calibration', 'query', 'spec', 'archive', 'annotations', 'output']:
            if isinstance(value.get(key), Path):
                value[key] = str(value[key])
        value = json.loads(json.dumps(value))
        path, base = None, Path(base_dir or Path.cwd()).resolve()
    else:
        path = Path(recipe).resolve()
        if base_dir is not None:
            raise ValueError('base_dir is only for Python mapping recipes')
        base = path.parent
        value = json.loads(path.read_text(), object_pairs_hook=_unique_object)
    if not isinstance(value, dict) or value.get('schema') != 'bremap-workflow-v1':
        raise ValueError('Expected bremap-workflow-v1 recipe')
    allowed = {'schema', 'name', 'reference', 'calibration', 'query', 'spec', 'archive',
               'method', 'policy', 'annotations', 'reliability', 'output'}
    if set(value)-allowed:
        raise ValueError('Unknown recipe keys: '+', '.join(sorted(set(value)-allowed)))
    if not isinstance(value.get('name'), str) or not value['name'].strip():
        raise ValueError('Recipe requires a nonempty name')
    required = {'query', 'archive'} if 'archive' in value else {'query', 'reference', 'calibration', 'spec', 'method', 'policy'}
    if not required <= value.keys():
        raise ValueError('Missing recipe keys: '+', '.join(sorted(required-value.keys())))
    if 'archive' in value and set(value) & {'reference', 'calibration', 'spec', 'method', 'policy', 'reliability'}:
        raise ValueError('Archive recipes reuse frozen settings; do not supply fitting options')
    if 'archive' not in value and value['method'] not in ('affinity', 'transport'):
        raise ValueError('Method must be affinity or transport')
    reliability = value.get('reliability')
    if reliability is not None:
        if not isinstance(reliability, dict) or set(reliability)-{'n_bootstrap', 'seed', 'leave_study_out'}:
            raise ValueError('Invalid reliability options')
        for key, default, minimum in [('n_bootstrap', 20, 1), ('seed', 42, 0)]:
            n = reliability.get(key, default)
            if type(n) is not int or n < minimum:
                raise ValueError('Invalid reliability '+key)
        if type(reliability.get('leave_study_out', False)) is not bool:
            raise ValueError('leave_study_out must be boolean')
    paths = {}
    for key in ['reference', 'calibration', 'query', 'spec', 'archive', 'annotations']:
        if key in value:
            if key == 'spec' and isinstance(value[key], dict):
                ReferenceSpec.from_dict(value[key])
                continue
            if not isinstance(value[key], str) or not value[key].strip():
                raise ValueError('Invalid path: '+key)
            paths[key] = (base/value[key]).resolve()
            if not paths[key].is_file():
                raise ValueError('Input file does not exist: '+str(paths[key]))
    return path, base, value, paths


def _owned_output(out):
    """Never remove a caller's unrelated or modified files during replacement."""
    try:
        manifest = json.loads((out/'workflow.json').read_text())
        if manifest['schema'] != 'bremap-workflow-result-v1' or manifest['status'] != 'complete':
            raise ValueError('Not complete')
        checks = manifest['output_sha256']
        files = list(out.rglob('*'))
        if any(p.is_symlink() for p in files):
            raise ValueError('Links are not owned artifacts')
        actual = {str(p.relative_to(out)) for p in files if p.is_file()}
        if actual != set(checks) | {'workflow.json'}:
            raise ValueError('Unexpected files')
        for name, digest in checks.items():
            if (out/name).resolve().is_relative_to(out) is False or file_hash(out/name) != digest:
                raise ValueError('Modified output')
    except (OSError, ValueError, KeyError, TypeError) as exc:
        raise ValueError('Overwrite requires a recognized, unmodified workflow output') from exc


def run_workflow(recipe, *, out_dir=None, overwrite=False, base_dir=None):
    """Run a Python mapping or JSON recipe; return its output Path.

    Python recipes accept a ReferenceSpec object and Path inputs. Their relative
    paths resolve against base_dir (default cwd); JSON paths against the file.
    An explicit
    out_dir is relative to the caller. Existing results survive failed replacement.
    Frozen embeddings must already be supplied with matched provenance in the H5AD.
    """
    started = time.perf_counter()
    path, base, config, inputs = _recipe(recipe, base_dir)
    if out_dir is None:
        if not isinstance(config.get('output'), str) or not config['output'].strip():
            raise ValueError('Specify output in the recipe or --out')
        out_dir = base/config['output']
    supplied_out = Path(out_dir).absolute()
    if supplied_out.is_symlink():
        raise ValueError('Output must not be a symbolic link')
    out = supplied_out.resolve()
    if any(p.is_relative_to(out) for p in ([path] if path is not None else [])+list(inputs.values())):
        raise ValueError('Inputs/recipe must not be inside the output directory')
    if out.exists():
        if not overwrite:
            raise ValueError('Output already exists; use --overwrite for owned workflow results')
        _owned_output(out)
    out.parent.mkdir(parents=True, exist_ok=True)
    lock = out.parent/('.'+out.name+'.bremap-lock')
    try:
        lock.mkdir()
    except FileExistsError as exc:
        raise ValueError('Another workflow owns the output lock: '+str(lock)) from exc
    stage = None
    try:
        # Repeat under the lock in case another run finished during validation.
        if out.exists():
            if not overwrite:
                raise ValueError('Output already exists')
            _owned_output(out)
        stage = Path(tempfile.mkdtemp(prefix='.'+out.name+'-', dir=out.parent))
        stages = []
        def step(name, function, *args, **kwargs):
            begin = time.perf_counter()
            result = function(*args, **kwargs)
            stages.append({'stage': name, 'wall_seconds': time.perf_counter()-begin})
            return result
        _write_json(stage/'recipe.json', config)
        input_hashes = {key: {'path': str(p), 'sha256': file_hash(p)} for key, p in inputs.items()}
        if 'archive' in inputs:
            calibrated = step('import_reference', import_reference, inputs['archive'], stage/'calibrated')
            meta, _, _, cal = load_reference(calibrated)
            if cal is None:
                raise ValueError('Archive must contain a calibrated reference')
            spec = ReferenceSpec.from_dict(meta['spec'])
        else:
            spec = (ReferenceSpec.from_dict(config['spec']) if isinstance(config['spec'], dict)
                    else ReferenceSpec.load(inputs['spec']))
        inspection = {}
        for role in ['reference', 'calibration', 'query']:
            if role in inputs:
                check = step('inspect_'+role, inspect_input, inputs[role], spec, role)
                inspection[role] = check
                if not check['valid']:
                    raise ValueError(role+': '+'; '.join(check['errors']))
        _write_json(stage/'input_inspection.json', inspection)
        if 'archive' not in inputs:
            model = step('fit_reference', fit_reference, inputs['reference'], spec, stage/'reference',
                         method=config['method'], policy=config['policy'])
            calibrated = step('calibrate', calibrate, model, inputs['calibration'], stage/'calibrated')
        mapped = step('map', map_query, calibrated, inputs['query'], stage/'mapping')
        meta, arrays, _, _ = load_reference(calibrated)
        _, programs, ids, donors, _, zero = transform_input(inputs['query'], meta, arrays)
        scores = pd.DataFrame(programs, columns=list(meta['programs']['anchors']))
        scores.insert(0, 'donor_id', donors); scores.insert(0, 'cell_id', ids)
        scores['qc_zero_library'] = zero
        (stage/'programs').mkdir()
        scores.to_csv(stage/'programs/cells.csv', index=False)
        summaries = []
        for donor, frame in scores.groupby('donor_id', sort=True):
            scored = frame[~frame.qc_zero_library]
            for name in meta['programs']['anchors']:
                summaries.append(dict(donor_id=donor, program=name, n_cells=len(frame), n_scored=len(scored),
                                      mean=scored[name].mean(), sd=scored[name].std(ddof=0)))
        pd.DataFrame(summaries).to_csv(stage/'programs/donors.csv', index=False)
        if 'annotations' in inputs:
            step('evaluate', evaluate, mapped, inputs['annotations'], stage/'evaluation')
        if config.get('reliability') is not None:
            step('reliability', assess_reliability, calibrated, inputs['reference'], inputs['query'],
                 stage/'reliability', **config['reliability'])
        step('export_reference', export_reference, calibrated, stage/'reference.zip')
        # Detect changed inputs rather than claiming provenance for a moving file.
        for key, p in inputs.items():
            if file_hash(p) != input_hashes[key]['sha256']:
                raise ValueError('Input changed during workflow: '+key)
        cells = pd.read_csv(stage/'mapping/cells.csv')
        try:
            import resource
            import sys
            rss = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/(1024 if sys.platform == 'darwin' else 1)
        except ImportError:
            rss = None
        record = dict(schema='bremap-workflow-result-v1', status='complete', name=config['name'],
                      recipe_sha256=(file_hash(path) if path is not None else hashlib.sha256(
                          json.dumps(config, sort_keys=True).encode()).hexdigest()),
                      recipe_origin='json_file' if path is not None else 'python_mapping',
                      inputs=input_hashes, output=str(out),
                      representation=meta['representation'], method=meta['method'], states=meta['classes'],
                      n_cells=len(cells), n_donors=cells.donor_id.nunique(),
                      acceptance=float(cells.accepted.mean()), stages=stages,
                      wall_seconds=time.perf_counter()-started, process_peak_rss_kib=rss,
                      memory_definition='Process-lifetime peak RSS; isolate one run per process for comparison.',
                      versions=_versions(), code_sha256={p.name: file_hash(p) for p in Path(__file__).parent.glob('*.py')},
                      limits=['Reference resemblance is not lineage ancestry or functional ground truth.',
                              'Programs are mean log-expression on the frozen reference gene universe.',
                              'No automatic encoder download/training or outcome-driven recipe selection.',
                              'Underlying build-time temporary paths are provenance, not artifact lookup paths.'])
        (stage/'report.md').write_text(
            f"# BReMAP workflow: {config['name']}\n\nCompleted {len(cells)} cells across {record['n_donors']} donor groups.\n\n"
            f"Representation: {meta['representation']['kind']}; method: {meta['method']}. "
            f"Cell acceptance: {record['acceptance']:.3f}.\n\n"
            "Outputs: mapping/cells.csv and donors.csv; programs/cells.csv and donors.csv; "
            "reference.zip for reuse; optional evaluation/ and reliability/.\n\n"
            "Conditional affinities describe resemblance to declared states, not validated cell fractions. "
            "Acceptance is a diagnostic, not proof of biological correctness. "
            "Consult input_inspection.json for missing genes and workflow.json for hashes/settings/runtime.\n")
        record['output_sha256'] = {str(p.relative_to(stage)): file_hash(p) for p in sorted(stage.rglob('*')) if p.is_file()}
        _write_json(stage/'workflow.json', record)
        if out.exists():
            with tempfile.TemporaryDirectory(prefix='.'+out.name+'-previous-', dir=out.parent) as tmp:
                backup = Path(tmp)/'result'
                out.rename(backup)
                try:
                    stage.rename(out)
                except BaseException:
                    backup.rename(out)
                    raise
        else:
            stage.rename(out)
        return out
    finally:
        if stage is not None and stage.exists():
            shutil.rmtree(stage)
        lock.rmdir()
