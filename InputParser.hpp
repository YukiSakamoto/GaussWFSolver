#pragma once
#include <string>
#include <map>
#include "MolecularSystem.hpp"

namespace wf_solver {

static const double factor_Angstrom2Bohr = 1.8897259885789;

const std::string Elements[] = {
    "0", 
    "H", "He", 
    "Li", "Be", "B", "C", "N", "O", "F", "Ne", 
    "Na", "Mg", "Al", "Si", "P" , "S", "Cl", "Ar",
    "K", "Ca", "Sc", "Ti", "V", "Cr", 
    "Mn", "Fe", "Co", "Ni", "Cu", "Zn", "Ga", "Ge", "As", "Se", "Br", "Kr",
    "",
};

class InputParser
{
public:
    typedef std::map<std::string, std::string> container_type;
    // Constructor
    InputParser(std::string &input_filename): input_filename_(input_filename) 
    {   this->parse();    }
    bool parse();
    void clear();
    const container_type &data() const
    {   return this->keyword_map_; };
private:
    std::string input_filename_;
    container_type keyword_map_;
};

void
read_xyz_file(const std::string xyz_filename, MolecularSystem &system, std::string &title);

};
