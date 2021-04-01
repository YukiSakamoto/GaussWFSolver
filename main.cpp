#include <boost/tokenizer.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/format.hpp>

#include "SCFBase.hpp"
#include "HartreeFock.hpp"
#include "InputParser.hpp"
#include <iostream>
#include <cstdio>


static const double factor_Angstrom2Bohr = 1.8897259885789;

int main(int argc, char **argv)
{
    wf_solver::MolecularSystem mol;
    if (argc < 2) {
        std::cerr << "Input file is not specifies" << std::endl;
        throw;
    }
    std::string input_filename(argv[1]);
    wf_solver::InputParser input(input_filename);

    std::printf("#==================================================\n");
    for(wf_solver::InputParser::container_type::const_iterator it = input.data().begin(); it != input.data().end(); it++) {
        std::printf("%-10s =  %s\n", it->first.c_str(), it->second.c_str() );
    }
    std::printf("#==================================================\n");

    mol.add_atom(6,  0., 0., 0.);
    mol.add_atom(1,  0.000000*factor_Angstrom2Bohr,    0.000000*factor_Angstrom2Bohr,    1.083010*factor_Angstrom2Bohr) ;
    mol.add_atom(1,  0.000000*factor_Angstrom2Bohr,    1.021071*factor_Angstrom2Bohr,   -0.361003*factor_Angstrom2Bohr) ;
    mol.add_atom(1,  0.884274*factor_Angstrom2Bohr,   -0.510536*factor_Angstrom2Bohr,   -0.361003*factor_Angstrom2Bohr) ;
    mol.add_atom(1, -0.884274*factor_Angstrom2Bohr,   -0.510536*factor_Angstrom2Bohr,   -0.361003*factor_Angstrom2Bohr) ;
    wf_solver::HartreeFock hf_scf(mol);
    hf_scf.compute();
    if (hf_scf.convergence() == true) {
        std::printf("Total Energy: %17.10f\n", hf_scf.energy() );
    }
    return 0;
}
