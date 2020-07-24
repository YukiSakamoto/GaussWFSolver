#pragma once

#include "common.hpp"
#include "MolecularSystem.hpp"
#include "SCFBase.hpp"
#include "libint_adapter.hpp"
#include "numerical_operation.hpp"

namespace wf_solver {
class HartreeFock: public SCFBase
{
public:
    HartreeFock(const MolecularSystem &mol) : system_(mol)
    {
    }
    ~HartreeFock() {
        
    }
    virtual bool compute();
    virtual bool gradient() {
        return true;
    }
private:
    MatrixXReal form_D(const MatrixXReal& C, size_t n_occ_orbitals) const;
    REAL check_scf_convergence(const MatrixXReal& D, const MatrixXReal& D_prev, REAL *maxdp) const;
    REAL calculate_E0(const MatrixXReal &D, const MatrixXReal &Hcore, const MatrixXReal &F) const;
    const MolecularSystem &system_;
};

};
