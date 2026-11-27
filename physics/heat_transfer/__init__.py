"""Shared heat-transfer API; correlations are supplied by Roles 3 and 4."""
# Load the kernel submodules first. Their names equal the public function names, and
# importing a submodule rebinds the package attribute to the module object. Binding the
# base functions afterwards keeps `heat_transfer.conduction(...)` callable on every call.
from . import conduction as _conduction_kernel, convection as _convection_kernel  # noqa: F401
from . import radiation as _radiation_kernel  # noqa: F401
from .base import (conduction, convection, radiation,
                   heat_transfer_coefficient)

__all__ = ["conduction", "convection", "radiation", "heat_transfer_coefficient"]
