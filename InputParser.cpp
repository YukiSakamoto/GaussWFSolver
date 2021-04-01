#include "InputParser.hpp"

#include <boost/tokenizer.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/format.hpp>
#include <iostream>
#include <string>
#include <fstream>
#include <map>

namespace wf_solver {

std::vector<std::string> split_line_by_space(const std::string &str)
{
    std::vector<std::string> retval;
    boost::char_separator<char> sep(" \t\n");
    typedef boost::tokenizer<boost::char_separator<char> > tokenizer;
    tokenizer tokens(str,sep);
    for(tokenizer::iterator it = tokens.begin(); it != tokens.end(); it++) {
        retval.push_back(*it);
    }
    return retval;
}

std::vector<std::string> split_line_by_eql(const std::string &str)
{
    std::vector<std::string> retval;
    boost::char_separator<char> sep("= \t\n");
    typedef boost::tokenizer<boost::char_separator<char> > tokenizer;
    tokenizer tokens(str,sep);
    for(tokenizer::iterator it = tokens.begin(); it != tokens.end(); it++) {
        retval.push_back(*it);
    }
    return retval;
}

std::string remove_comment(const std::string &line)
{
    std::string ret;
    for(std::string::const_iterator it = line.begin(); it != line.end(); it++) {
        if (*it == '#') { break; }
        ret.push_back(*it);
    }
    return ret;
}

void InputParser::clear()
{
    this->keyword_map_.empty();
}

bool InputParser::parse()
{
    this->clear();
    std::ifstream ifs(input_filename_.c_str(), std::ios::in);
    if (!ifs.is_open()) {
        std::cerr << "Error!: maybe the file doesn't exist: " << input_filename_ << std::endl;
        throw;
    }
    size_t line_count = 0;
    while ( !ifs.eof() ) {
        std::string line_buf;
        std::getline(ifs, line_buf); 
        line_count++;

        line_buf =  remove_comment(line_buf);
        std::vector<std::string> tokens = split_line_by_eql(line_buf);
        if (tokens.size() == 0) {
            // just comment
        } else if (tokens.size() == 2) {
            std::string name, val;
            std::transform(tokens[0].begin(), tokens[0].end(), std::back_inserter(name), ::tolower);
            std::transform(tokens[1].begin(), tokens[1].end(), std::back_inserter(val) , ::tolower);
            this->keyword_map_.insert(std::make_pair(name, val));
        } else {
            throw;
        }
    }
    return true;
}

const std::string Elements[] = {
    "0", 
    "H", "He", 
    "Li", "Be", "B", "C", "N", "O", "F", "Ne", 
    "Na", "Mg", "Al", "Si", "P" , "S", "Cl", "Ar",
    "K", "Ca", "Sc", "Ti", "V", "Cr", 
    "Mn", "Fe", "Co", "Ni", "Cu", "Zn", "Ga", "Ge", "As", "Se", "Br", "Kr",
    "",
};

int get_atomic_number(const std::string &element)
{
    for(int i = 0; Elements[i] != ""; i++) {
        if (Elements[i] == element) {
            return i;
        }
    }
    return -1;
}

void
read_xyz_file(const std::string xyz_filename, MolecularSystem &system, std::string &title)
{
    std::ifstream ifs(xyz_filename.c_str(), std::ios::in);
    if (!ifs.is_open()) {
        std::cerr << "Error! read_xyz_file(): maybe ths file doesn't exist: " << xyz_filename << std::endl;
        throw;
    }

    size_t num_atoms = 0;
    size_t line_count = 0;
    while ( !ifs.eof() ) {
        std::string line_buf;
        std::vector<std::string> tokens;
        std::getline(ifs, line_buf); 
        line_count++;

        if (line_count == 1) {
            tokens = split_line_by_space(line_buf);
            // the first line is the number of atoms
            if (tokens.size()  == 1) { 
                num_atoms = boost::lexical_cast<int>(tokens[0]);
            } else {
                throw; 
            }
        } else if (line_count == 2) {
            title = line_buf;
        } else {
            tokens = split_line_by_space(line_buf);
            if (tokens.size() == 4) {
                int atomic_number = get_atomic_number(tokens[0]);
                REAL x = factor_Angstrom2Bohr * ( boost::lexical_cast<REAL>(tokens[1]) );
                REAL y = factor_Angstrom2Bohr * ( boost::lexical_cast<REAL>(tokens[2]) );
                REAL z = factor_Angstrom2Bohr * ( boost::lexical_cast<REAL>(tokens[3]) );
                system.add_atom(atomic_number, x, y, z);
            } else {
                continue;
            }
        }
    }
    if (system.size() != num_atoms) {
        std::cerr << "Number of atoms are not matched" << std::endl;
        throw;
    }

    return;
}

}
