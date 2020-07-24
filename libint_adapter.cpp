#include "libint_adapter.hpp"
#include <cstdio>

namespace wf_solver {

void integral_engine_init(bool diagnostic)
{
    libint2::initialize(diagnostic);
}
 
void integral_engine_finalize()
{
    libint2::finalize();
}

std::vector<libint2::Atom> convert_molecules(const MolecularSystem &mol_sys)
{
    std::vector<libint2::Atom> ret;
    for(int i = 0; i < mol_sys.size(); i++) {
        libint2::Atom atom;
        std::array<REAL,3> pos = mol_sys.atom_position(i);
        atom.atomic_number = mol_sys.atomic_number(i);
        atom.x = pos[0];
        atom.y = pos[1];
        atom.z = pos[2];
        ret.push_back(atom);
    }
    return ret;
}

std::vector<size_t> map_shell_to_basis_function(const std::vector<libint2::Shell> &shells) {
    // map the index of the shell
    std::vector<size_t> result;
    size_t n = 0; 
    for(size_t i = 0; i < shells.size(); i++) {
        result.push_back(n);
        n += shells[i].size();
    }
    return result;
}

MatrixXReal compute_1body_ints(const std::vector<libint2::Shell> &shells, 
        const libint2::Operator obtype, const std::vector<libint2::Atom> &atoms)
{
    size_t nbasis = 0;
    for(size_t i = 0; i < shells.size(); i++) {
        nbasis += shells[i].size();
    }
    MatrixXReal result = MatrixXReal::Zero(nbasis, nbasis);
    std::vector<size_t> shell2bf = map_shell_to_basis_function(shells);

    libint2::Engine engine(obtype, max_nprim(shells), max_l(shells));

    if (obtype == libint2::Operator::nuclear) {
        std::vector<std::pair<REAL, std::array<REAL, 3>>> q;
        for(size_t i = 0; i < atoms.size(); i++) {
            q.push_back({static_cast<REAL>(atoms[i].atomic_number), {{atoms[i].x, atoms[i].y, atoms[i].z}} });
        }
        engine.set_params(q);
    }
    const auto& buf = engine.results();

    for(int i = 0; i < shells.size(); i++) {
        size_t bf1 = shell2bf[i];  // index of the first basis function in this shell
        size_t n1 = shells[i].size();

        for(int j = 0; j <= i; j++) {
            size_t bf2 = shell2bf[j];
            size_t n2  = shells[j].size();

            engine.compute(shells[i], shells[j]);
            Eigen::Map<const MatrixXReal> buf_mat(buf[0], n1, n2);
            result.block(bf1, bf2, n1, n2) = buf_mat;
            if (i != j)  {
                result.block(bf2, bf1, n2, n1) = buf_mat.transpose();
            }
        }
    }
    return result;
}

MatrixXReal compute_overlap_matrix(const std::vector<libint2::Shell> &shells)
{
    std::vector<libint2::Atom> __dummy_atoms;
    return compute_1body_ints(shells, libint2::Operator::overlap, __dummy_atoms);
}

MatrixXReal compute_kinetic_matrix(const std::vector<libint2::Shell> &shells)
{
    std::vector<libint2::Atom> __dummy_atoms;
    return compute_1body_ints(shells, libint2::Operator::kinetic, __dummy_atoms);
}

MatrixXReal compute_nuclear_attraction_matrix(
        const std::vector<libint2::Shell> &shells, const std::vector<libint2::Atom> &atoms)
{
    return compute_1body_ints(shells, libint2::Operator::nuclear, atoms);
}

MatrixXReal compute_fock_2body_matrix(const std::vector<libint2::Shell> &shells, const MatrixXReal &D)
{
    size_t nbasis = 0;
    for(size_t i = 0; i < shells.size(); i++) {
        nbasis += shells[i].size();
    }
    MatrixXReal G = MatrixXReal::Zero(nbasis,nbasis);
    libint2::Engine engine(libint2::Operator::coulomb, max_nprim(shells), max_l(shells));
    const auto& buf = engine.results();
    std::vector<size_t> shell2bf = map_shell_to_basis_function(shells);
    for(size_t s1 = 0; s1 < shells.size(); s1++) {
        size_t bf1_first = shell2bf[s1]; // the index of the basis function corresponds to shell i
        size_t n1  = shells[s1].size();
        for(size_t s2 = 0; s2 < shells.size(); s2++) {
            size_t bf2_first = shell2bf[s2];
            size_t n2  = shells[s2].size();
            for(size_t s3 = 0; s3 < shells.size(); s3++) {
                size_t bf3_first = shell2bf[s3];
                size_t n3  = shells[s3].size();
                for(size_t s4 = 0; s4 < shells.size(); s4++) {
                    size_t bf4_first = shell2bf[s4];
                    size_t n4  = shells[s4].size();

                    //========================================
                    //  Coulomb Integral
                    //========================================
                    engine.compute(shells[s1], shells[s2], shells[s3], shells[s4]);
                    const auto* buf_1234 = buf[0];
                    if (buf_1234 == nullptr) { continue; }

                    size_t f1234 = 0;
                    for(size_t f1=0; f1 != n1; f1++) {
                        const size_t bf1 = f1 + bf1_first;
                        for(size_t f2 = 0; f2 != n2; f2++) {
                            const size_t bf2 = f2 + bf2_first;
                            for(size_t f3 = 0; f3 != n3; f3++) {
                                const size_t bf3 = f3 + bf3_first;
                                for(size_t f4 = 0; f4 != n4; f4++, f1234++) {
                                    const size_t bf4 = f4 + bf4_first;
                                    G(bf1, bf2) += D(bf3, bf4) * buf_1234[f1234];
                                }
                            }
                        }
                    }
                    //========================================
                    //  Exchange Integral
                    //========================================
                    engine.compute(shells[s1], shells[s3], shells[s2], shells[s4]);
                    const auto* buf_1324 = buf[0];
                    size_t f1324 = 0;
                    for(size_t f1 = 0; f1 != n1; f1++) {
                        const size_t bf1 = f1 + bf1_first;
                        for(size_t f3 = 0; f3 != n3; f3++) {
                            const size_t bf3 = f3 + bf3_first;
                            for(size_t f2 = 0; f2 != n2; f2++) {
                                const size_t bf2 = f2 + bf2_first;
                                for(size_t f4 = 0; f4 != n4; f4++, f1324++) {
                                    const size_t bf4 = f4 + bf4_first;
                                    G(bf1, bf2) -= D(bf3, bf4) * 0.5 * buf_1324[f1324];
                                }
                            }
                        }
                    }

                }
            }
        }
    }
    return G;
}

};
