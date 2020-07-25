#include "HartreeFock.hpp"
#include <cstdio>

namespace wf_solver {

MatrixXReal HartreeFock::form_D(const MatrixXReal& C, size_t n_occ_orbitals) const
{
    // Szabo. pp. 139 (3.145)
    size_t row = C.rows();
    size_t col = C.cols();
    MatrixXReal D = MatrixXReal::Zero(row, col);
    for(size_t u = 0 ; u < row; u++) {
        for(size_t v = 0; v < col; v++) {
            for(size_t a = 0 ; a < n_occ_orbitals; a++) {
                D(u,v) += 2.0 * C(u,a) * C(v,a);
            }
        }
    }
    return D;
}

REAL HartreeFock::check_scf_convergence(
        const MatrixXReal& D, const MatrixXReal& D_prev, REAL *maxdp) const
{
    if (D.rows() != D_prev.rows() || D.cols() != D_prev.cols()) {
        throw;
    }
    MatrixXReal D_diff = D - D_prev;
    REAL maxdp_acc = 0.;
    REAL rmsdp_acc = 0.;

    size_t row = D_diff.rows();
    size_t col = D_diff.cols();
    for(size_t i = 0; i < row; i++) {
        for(size_t j = 0; j < col; j++) {
            rmsdp_acc += std::pow(D_diff(i,j), 2);
            maxdp_acc = std::max(maxdp_acc, std::abs(D_diff(i,j)) );
        }
    }
    if (maxdp != NULL) {
        *maxdp = maxdp_acc;
    }
    REAL rmsdp = std::sqrt( rmsdp_acc/(row*col) );
    return rmsdp;
}

REAL HartreeFock::calculate_E0(
        const MatrixXReal &D, const MatrixXReal &Hcore, const MatrixXReal &F) const
{
    // Szabo. pp.150 (3.184): 
    REAL E0 = 0.;
    MatrixXReal H_F = Hcore + F;
    size_t row = H_F.rows();
    size_t col = H_F.cols();
    for(size_t u = 0; u < row; u++) {
        for(size_t v = 0; v < col; v++) {
            E0 += D(v,u) * H_F(u,v);
        }
    }
    E0 *= 0.5;
    return E0;
}

bool HartreeFock::compute() {
    integral_engine_init(true);

    BasisFunctionsImpl bf(system_, "sto-3g");
    size_t nao = bf.nbasis(); // number of basis_set. 
    int num_occ_orbitals = system_.num_electrons() / 2;
    std::cout << num_occ_orbitals << std::endl;

    MatrixXReal S = bf.compute_overlap_matrix();
    MatrixXReal T = bf.compute_kinetic_matrix();
    MatrixXReal V = bf.compute_nuclear_attraction_matrix();
    MatrixXReal Hcore = T + V;
    std::cout << "1 electron integrals done" << std::endl;
    MatrixXReal X = canonical_orthogonalization(S);
    REAL nei = this->system_.nuclear_repulsion();
    std::cout << "NEI: " << nei << std::endl;

    std::printf("%2s\t%10s\t%10s\t\n", "i", "Energy", "MaxDP");
    // Initial guess is 0, at present.
    MatrixXReal D = MatrixXReal::Zero(nao, nao);
    for(int i = 0; i < 100; i++) {
        MatrixXReal G = bf.compute_fock_2body_matrix(D);
        MatrixXReal F = Hcore + G;
        MatrixXReal F_prim = X.adjoint() * F * X;
        Eigen::SelfAdjointEigenSolver<MatrixXReal> es(F_prim);
        if (es.info() != Eigen::Success) {  throw;  }

        // Energies and New coefficients 
        VectorXReal e = es.eigenvalues();
        MatrixXReal C_new_prime = es.eigenvectors();
        MatrixXReal C_new = X * C_new_prime;
        MatrixXReal D_new = this->form_D(C_new, num_occ_orbitals);
        REAL maxdp = 0.;
        REAL rmsdp = 0.;
        REAL E0  = this->calculate_E0(D, Hcore, F);
        rmsdp = this->check_scf_convergence(D_new, D, &maxdp);

        std::printf("%02d\t%10.8f\t%10.8f\t \n", i, E0, maxdp);
        if (rmsdp < 1.0e-8) {
            this->scf_convergence_ = true;
            break;
        } else {
            D = D_new;
        }
    }

    integral_engine_finalize();
    return scf_convergence_;
};

};
