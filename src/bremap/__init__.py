"""BReMAP provisional research software for normal-reference analog mapping."""

from .core import ReferenceMapper
from .transport import transport_plan, transport_predict

__version__ = "0.4.1"
__all__ = ["ReferenceMapper", "transport_plan", "transport_predict"]

from .specification import ReferenceSpec
from .api import (Artifact, AnalysisResult, ReferenceModel, CalibratedReference, MappingResult,
                  ReliabilityResult, EvaluationResult, fit_reference, calibrate, map_query,
                  assess_reliability, evaluate)

__all__ += ['ReferenceSpec', 'Artifact', 'ReferenceModel', 'CalibratedReference', 'MappingResult',
            'AnalysisResult', 'ReliabilityResult', 'EvaluationResult',
            'fit_reference', 'calibrate', 'map_query', 'assess_reliability', 'evaluate']

from .configuration import read_gene_programs, create_spec, inspect_input
from .resources_io import export_reference, import_reference, export_embeddings

__all__ += ['read_gene_programs', 'create_spec', 'inspect_input',
            'export_reference', 'import_reference', 'export_embeddings']
