# Copyright (c) 2009-2025 The Regents of the University of Michigan.
# Part of HOOMD-blue, released under the BSD 3-Clause License.

"""ExpandedYukawa pair potential."""

# Import the C++ module.
from hoomd.data.parameterdicts import TypeParameterDict
from hoomd.data.typeparam import TypeParameter
from hoomd.expanded_Yukawa import _expanded_Yukawa

# Import the hoomd Python package and other necessary components.
from hoomd.md import pair


class ExpandedYukawa(pair.Pair):
    """ExpandedYukawa pair potential."""

    # set static class data
    _ext_module = _expanded_Yukawa
    _cpp_class_name = 'PotentialPairExpandedYukawa'
    _accepted_modes = ('none', 'shift', 'xplor')

    def __init__(self, nlist, default_r_cut=None, default_r_on=0.0, mode='none'):
        super().__init__(nlist, default_r_cut, default_r_on, mode)
        params = TypeParameter(
            'params',
            'particle_types',
            TypeParameterDict(epsilon=float, kappa=float, delta=float, len_keys=2),
        )
        self._add_typeparam(params)
