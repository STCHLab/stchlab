"""STCHLab: tools for solar thermochemical redox-cycle analysis.

The package exposes reusable modules for material thermodynamics, reactor-level
models, and system-level analysis in both counter-current and fixed-bed
configurations.
"""

from .materials import BCZM, CEF, CeO2
from .core import STCHLab

__all__ = ['STCHLab', 'CEF', 'CeO2', 'BCZM']
