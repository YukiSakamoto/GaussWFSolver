#pragma once
#include "common.hpp"
#include <Eigen/Eigen> 
#include <Eigen/Core>

namespace wf_solver {
inline
REAL squared_norm( const REAL x1, const REAL y1, const REAL z1, 
        const REAL x2, const REAL y2, const REAL z2) {
    return std::pow(x2-x1, 2) + std::pow(y2-y1, 2) + std::pow(z2-z1, 2);
}

inline
REAL norm( const REAL x1, const REAL y1, const REAL z1, 
        const REAL x2, const REAL y2, const REAL z2) {
    return std::sqrt( squared_norm(x1,y1,z1,x2,y2,z2) );
}

inline 
MatrixXReal canonical_orthogonalization(const MatrixXReal& S) 
{
    // Szabo. pp.144 (3.169 - 3.172)
    // Overlap Integral Matrix (S) should be the 
    //  self-adjoint matrix by the definition.
    Eigen::SelfAdjointEigenSolver<MatrixXReal> es(S);
    if (es.info() != Eigen::Success) {  throw;  }

    VectorXReal l = es.eigenvalues();
    MatrixXReal U = es.eigenvectors();

    size_t row = S.rows();
    size_t col = S.cols();
    MatrixXReal X = MatrixXReal::Zero(row, col);
    for (size_t i = 0; i < row; i++) {
        for(size_t j = 0; j < col; j++) {
            X(i,j) = U(i,j) / std::sqrt(l[j]);
        }
    }
    return X;
}

}
