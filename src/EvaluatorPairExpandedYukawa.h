// Copyright (c) 2009-2025 The Regents of the University of Michigan.
// Part of HOOMD-blue, released under the BSD 3-Clause License.

#ifndef __PAIR_EVALUATOR_EXPANDEDYUKAWA_H__
#define __PAIR_EVALUATOR_EXPANDEDYUKAWA_H__

#ifndef __HIPCC__
#include <string>
#endif

#include "hoomd/HOOMDMath.h"

/*! \file EvaluatorPairExpandedYukawa.h
    \brief Defines the pair evaluator class for the expanded Yukawa potential
*/

// need to declare these class methods with __device__ qualifiers when building in nvcc
// DEVICE is __host__ __device__ when included in nvcc and blank when included into the host
// compiler
#ifdef __HIPCC__
#define DEVICE __device__
#define HOSTDEVICE __host__ __device__
#else
#define DEVICE
#define HOSTDEVICE
#endif

namespace hoomd
    {
namespace md
    {

class EvaluatorPairExpandedYukawa
    {
    public:
    //! Define the parameter type used by this pair potential evaluator
    struct param_type
        {
        Scalar epsilon;
        Scalar kappa;
        Scalar delta;

        DEVICE void load_shared(char*& ptr, unsigned int& available_bytes) { }

        HOSTDEVICE void allocate_shared(char*& ptr, unsigned int& available_bytes) const { }

#ifdef ENABLE_HIP
        //! Set CUDA memory hints
        void set_memory_hint() const
            {
            // default implementation does nothing
            }
#endif

#ifndef __HIPCC__
        param_type() : epsilon(0), kappa(0), delta(0) { }

        param_type(pybind11::dict v, bool managed = false)
            {
            epsilon = v["epsilon"].cast<Scalar>();
            kappa = v["kappa"].cast<Scalar>();
            delta = v["delta"].cast<Scalar>();
            }

        pybind11::dict asDict()
            {
            pybind11::dict v;
            v["epsilon"] = epsilon;
            v["kappa"] = kappa;
            v["delta"] = delta;
            return v;
            }
#endif
        }
#if HOOMD_LONGREAL_SIZE == 32
        __attribute__((aligned(8)));
#else
        __attribute__((aligned(16)));
#endif

    //! Constructs the pair potential evaluator
    /*! \param _rsq Squared distance between the particles
        \param _rcutsq Squared distance at which the potential goes to 0
        \param _params Per type pair parameters of this potential
    */
    DEVICE EvaluatorPairExpandedYukawa(Scalar _rsq, Scalar _rcutsq, const param_type& _params)
        : rsq(_rsq), rcutsq(_rcutsq), epsilon(_params.epsilon), kappa(_params.kappa), delta(_params.delta)
        {
        }

    //! Expanded Yukawa doesn't use charge
    DEVICE static bool needsCharge()
        {
        return false;
        }
    //! Accept the optional charge value
    /*! \param qi Charge of particle i
        \param qj Charge of particle j
    */
    DEVICE void setCharge(Scalar qi, Scalar qj) { }

    //! Evaluate the force and energy
    /*! \param force_divr Output parameter to write the computed force divided by r.
        \param pair_eng Output parameter to write the computed pair energy
        \param energy_shift If true, the potential must be shifted so that
        V(r) is continuous at the cutoff
        \note There is no need to check if rsq < rcutsq in this method.
        Cutoff tests are performed in PotentialPair.

        \return True if they are evaluated or false if they are not because
        we are beyond the cutoff
    */
    DEVICE bool evalForceAndEnergy(Scalar& force_divr, Scalar& pair_eng, bool energy_shift)
        {
        // compute the force divided by r in force_divr
        if (rsq < rcutsq && epsilon != 0)
            {
            Scalar rinv = fast::rsqrt(rsq);
            Scalar r = Scalar(1.0) / rinv;
            Scalar rmd = r - delta;
            Scalar rmdinv = Scalar(1.0) / rmd;
            Scalar rmd_sq = (r - delta) * (r - delta);
            Scalar rmd2inv = Scalar(1.0) / rmd_sq;
            Scalar exp_val = fast::exp(-kappa * rmd);

            force_divr = epsilon * exp_val * rmd2inv * (rinv + kappa * rmd / r);
            pair_eng = epsilon * exp_val * rmdinv;

            if (energy_shift)
                {
                Scalar rcutinv = fast::rsqrt(rcutsq);
                Scalar rcut = Scalar(1.0) / rcutinv;
                Scalar rcutmd = rcut - delta;
                Scalar rcutmdinv = Scalar(1.0) / rcutmd;
                pair_eng -= epsilon * fast::exp(-kappa * rcutmd) * rcutmdinv;
                }
            return true;
            }
        else
            return false;
        }

    //! Expanded Yukawa doesn't eval LRC integrals
    DEVICE Scalar evalPressureLRCIntegral()
        {
        return 0;
        }

    //! Expanded Yukawa doesn't eval LRC integrals
    DEVICE Scalar evalEnergyLRCIntegral()
        {
        return 0;
        }

#ifndef __HIPCC__
    //! Get the name of this potential
    /*! \returns The potential name.
     */
    static std::string getName()
        {
        return std::string("expanded_yukawa");
        }

    std::string getShapeSpec() const
        {
        throw std::runtime_error("Shape definition not supported for this pair potential.");
        }
#endif

    protected:
    Scalar rsq;     //!< Stored rsq from the constructor
    Scalar rcutsq;  //!< Stored rcutsq from the constructor
    Scalar epsilon; //!< Stored epsilon from the constructor
    Scalar kappa;   //!< Stored kappa from the constructor
    Scalar delta;   //!< Stored delta from the constructor
    };

    } // end namespace md
    } // end namespace hoomd

#endif // __PAIR_EVALUATOR_EXPANDEDYUKAWA_H__
