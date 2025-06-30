######################################################################
# Ghostly: Ghost atom bonded modifications for alchemical free energy
# simulations.
#
# Copyright: 2024-2025
#
# Authors: The OpenBioSim Team <team@openbiosim.org>
#
# Ghsotly is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# Ghostly is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with Ghostly. If not, see <http://www.gnu.org/licenses/>.
#####################################################################

# Make sure we used the mixed API so we can use BioSimSpace.
try:
    import sire as _sr

    _sr.use_mixed_api(support_old_module_names=False)
    _sr.convert.supported_formats()

    del _sr
except ImportError:
    pass

from ._ghostly import *
