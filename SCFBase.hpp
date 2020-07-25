#pragma once
#include "common.hpp"

namespace wf_solver
{

class SCFBase
{
public:
    SCFBase() {
    }
    virtual 
    ~SCFBase() {
    }
    bool convergence() const {
        return scf_convergence_;
    }
    REAL energy() const {
        return this->energy_;
    }
    virtual bool compute() = 0;
    virtual bool gradient() = 0; 
    
protected:
    bool scf_convergence_ = false;
    REAL energy_ = 0.;
};
};
