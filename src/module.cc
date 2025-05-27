// Copyright (c) 2009-2025 The Regents of the University of Michigan.
// Part of HOOMD-blue, released under the BSD 3-Clause License.

// TODO: Include the header files of classes that will be exported to Python.
#include "EvaluatorPairExpandedYukawa.h"

#include "hoomd/md/PotentialPair.h"
#include <pybind11/pybind11.h>

#ifdef ENABLE_HIP
#include "hoomd/md/PotentialPairGPU.h"
#endif

namespace hoomd
    {
namespace md
    {

// TODO: Set the name of the python module to match ${COMPONENT_NAME} (set in
// CMakeLists.txt), prefixed with an underscore.
PYBIND11_MODULE(_expanded_Yukawa, m)
    {
    // TODO: Call export_Class(m) for each C++ class to be exported to Python.

    detail::export_PotentialPair<EvaluatorPairExpandedYukawa>(m, "PotentialPairExpandedYukawa");
#ifdef ENABLE_HIP
    // TODO: Call export_ClassGPU(m) for each GPU enabled C++ class to be exported
    // to Python.
    detail::export_PotentialPairGPU<EvaluatorPairExample>(m, "PotentialPairExpandedYukawaGPU");
#endif
    }

    } // end namespace md
    } // end namespace hoomd
