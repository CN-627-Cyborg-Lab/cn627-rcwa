"""
CN-627: Lightweight 1D-RCWA validator.
CN-627: 轻量级一维RCWA验证器。
"""

__version__ = "0.1.0"

from .core import rcwa_te, rcwa_spectrum
from .tmm import tmm_slab
from .materials import get_si_eps, MATERIALS