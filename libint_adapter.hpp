#pragma once
#include <vector>
#include <libint2.hpp>
#include <Eigen/Core>
#include <libint2/diis.h>

#include "common.hpp"
#include "MolecularSystem.hpp"
#include "BasisFunctions.hpp"


namespace wf_solver {

typedef libint2::DIIS<MatrixXReal> DIIS;

class BasisFunctionsImpl: public BasisFunctions
{
public:
    BasisFunctionsImpl(const MolecularSystem &system, const std::string &basis_set_type): 
        system_(system), basis_set_type_(basis_set_type)
    { this->update_molecules_();}
    virtual ~BasisFunctionsImpl() {}
    virtual MatrixXReal compute_overlap_matrix() const;
    virtual MatrixXReal compute_kinetic_matrix() const;
    virtual MatrixXReal compute_nuclear_attraction_matrix() const;
    virtual MatrixXReal compute_fock_2body_matrix(const MatrixXReal &D) const;
    virtual size_t nbasis() const;

    virtual MatrixXReal compute_fock_2body_matrix_parallel(const MatrixXReal &D, size_t threads = 1) const;
private:
    const MolecularSystem &system_;
    const std::string &basis_set_type_;

    std::vector<libint2::Atom> atom_;       // this is used for libint2;
    libint2::BasisSet shells_;
    std::vector<size_t> shell2bf_;          // first index of the basis function 

    void update_molecules_();
    void map_shell_to_basis_function_();
    MatrixXReal compute_1body_ints(const libint2::Operator obtype) const;
};

void integral_engine_init(bool diagnostic = false);
void integral_engine_finalize();

inline size_t max_nprim(const std::vector<libint2::Shell> &shells) {
    size_t n = 0;
    for(size_t i = 0; i < shells.size(); i++) {
        n = std::max(shells[i].nprim(), n);
    }
    return n;
}

inline 
int max_l(const std::vector<libint2::Shell> &shells) {
    int l = 0;
    for(size_t i = 0; i < shells.size(); i++) {
        for(size_t j = 0; j < shells[i].ncontr(); j++) {
            l = std::max(shells[i].contr[j].l, l);
        }
    }
    return l;
}

};
