#pragma once

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
    virtual bool compute() = 0;
    virtual bool gradient() = 0; 
    bool scf_convergence_ = false;
};
};
