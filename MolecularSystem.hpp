#pragma once

#include "common.hpp"
#include "numerical_operation.hpp"
#include <vector>
#include <array>
#include <numeric>
#include <algorithm>

#include <libint2.hpp>



namespace wf_solver {

class MolecularSystem
{
public:
    MolecularSystem(const int charge = 0, const int spin_multiplicity = 1) : 
        charge_(charge), spin_multiplicity_(spin_multiplicity)
    { }

    MolecularSystem( const std::vector<int> &atom_list, 
            const std::vector<REAL> &x_pos, const std::vector<REAL> &y_pos, const std::vector<REAL> &z_pos, 
            const int charge = 0, const int spin_multiplicity = 1) :
        charge_(charge), spin_multiplicity_(spin_multiplicity)
    {
        if (atom_list.size() != x_pos.size() 
                || atom_list.size() != y_pos.size() 
                || atom_list.size() != z_pos.size()) {
            throw;
        }
        std::copy( atom_list.begin(), atom_list.end(), std::back_inserter(this->atom_list_) );
        std::copy( x_pos.begin(), x_pos.end(), std::back_inserter(this->x_positions_) );
        std::copy( y_pos.begin(), y_pos.end(), std::back_inserter(this->y_positions_) );
        std::copy( z_pos.begin(), z_pos.end(), std::back_inserter(this->z_positions_) );
    }

    void add_atom(const int atom_number, const REAL x, const REAL y, const REAL z)
    {
        atom_list_.push_back(atom_number);
        x_positions_.push_back(x);
        y_positions_.push_back(y);
        z_positions_.push_back(z);

        update_internal_();
    }

    void set_charge(const int charge) 
    {   
        this->charge_ = charge; 
        update_internal_();
    }

    int charge() const
    {   return this->charge_;   }

    std::array<REAL,3> atom_position(const size_t index) const
    {
        if (index < atom_list_.size() ) {
            std::array<REAL,3> ret = { x_positions_[index], y_positions_[index], z_positions_[index] };
            return ret;
        } 
        throw;
    }

    int atomic_number(const size_t index) const
    {
        if (index < atom_list_.size() ) {
            return atom_list_[index];
        } 
        throw;
    }

    size_t size(void) const
    {
        return this->atom_list_.size();
    }

    //std::array<REAL,3> operator[](const size_t index) const
    //{   return this->atom_position(index);  }
    int num_electrons() const
    {
        return this->num_electrons_;
    }

    REAL nuclear_repulsion() const
    {
        REAL ret = 0.;
        for(size_t i = 0; i < this->atom_list_.size(); i++) {
            for(size_t j = i + 1; j < this->atom_list_.size();  j++) {
                if (i != j) {
                    REAL r_ij = norm(x_positions_[i],y_positions_[i],z_positions_[i], x_positions_[j],y_positions_[j],z_positions_[j]);
                    ret +=  (atom_list_[i] * atom_list_[j]) / r_ij;
                }
            }
        }
        return ret;
    }
    
private:
    void update_internal_()
    {
        if (spin_multiplicity_ <= 0) {  throw;  }
        int n_spin = (spin_multiplicity_-1) / 2;

        num_electrons_ = std::accumulate( atom_list_.begin(), atom_list_.end(), 0 );
        num_alpha_electrons_ = (num_electrons_ + n_spin) / 2;
        num_beta_electrons_  = (num_electrons_ - n_spin) / 2;
    }

private:
    // SOA(Structure of Array) is  used rather than AOS(Array of Structure) to improve cache hit ratio
    std::vector<REAL> x_positions_;
    std::vector<REAL> y_positions_;
    std::vector<REAL> z_positions_;
    std::vector<int>  atom_list_;
    
    int num_electrons_;
    int num_alpha_electrons_;
    int num_beta_electrons_;
    int charge_;
    int spin_multiplicity_;
};


};
