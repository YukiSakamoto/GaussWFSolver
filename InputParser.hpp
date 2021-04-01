#pragma once
#include <string>
#include <map>
#include "MolecularSystem.hpp"

namespace wf_solver {

static const double factor_Angstrom2Bohr = 1.8897259885789;

class InputParser
{
public:
    typedef std::map<std::string, std::string> container_type;
    // Constructor
    InputParser(std::string &input_filename): input_filename_(input_filename) 
    {   this->parse();    }
    bool parse();
    void clear();
    const container_type &data(){   return this->keyword_map_; };
private:
    std::string input_filename_;
    container_type keyword_map_;
};

void
read_xyz_file(const std::string xyz_filename, MolecularSystem &system, std::string &title);

};
