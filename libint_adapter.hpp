#pragma once
#include <vector>
#include <libint2.hpp>
#include <Eigen/Core>

#include "common.hpp"
#include "MolecularSystem.hpp"

namespace wf_solver {

void integral_engine_init(bool diagnostic = false);
void integral_engine_finalize();

inline size_t max_nprim(const std::vector<libint2::Shell> &shells) {
    size_t n = 0;
    for(size_t i = 0; i < shells.size(); i++) {
        n = std::max(shells[i].nprim(), n);
    }
    return n;
}

inline int max_l(const std::vector<libint2::Shell> &shells) {
    int l = 0;
    for(size_t i = 0; i < shells.size(); i++) {
        for(size_t j = 0; j < shells[i].size(); j++) {
            l = std::max(shells[i].contr[j].l, l);
        }
    }
    return l;
}

std::vector<libint2::Atom> convert_molecules(const MolecularSystem &mol_sys);
std::vector<size_t> map_shell_to_basis_function(const std::vector<libint2::Shell> &shells);

MatrixXReal compute_1body_ints(const std::vector<libint2::Shell> &shells, 
        const libint2::Operator obtype, const std::vector<libint2::Atom> &atoms);
MatrixXReal compute_fock_2body_matrix(const std::vector<libint2::Shell> &shells, const MatrixXReal &D);

MatrixXReal compute_overlap_matrix(const std::vector<libint2::Shell> &shells);
MatrixXReal compute_kinetic_matrix(const std::vector<libint2::Shell> &shells);
MatrixXReal compute_nuclear_attraction_matrix(
        const std::vector<libint2::Shell> &shells, const std::vector<libint2::Atom> &atoms);
};
