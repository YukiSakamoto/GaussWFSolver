
namespace wf_solver {

class BasisFunctions
{
public:
    BasisFunctions() {}
    virtual ~BasisFunctions(){};
    virtual MatrixXReal compute_overlap_matrix() const = 0;
    virtual MatrixXReal compute_kinetic_matrix() const = 0;
    virtual MatrixXReal compute_nuclear_attraction_matrix() const = 0;
    virtual MatrixXReal compute_fock_2body_matrix(const MatrixXReal &D) const = 0;
    virtual size_t nbasis() const = 0;
};


};
