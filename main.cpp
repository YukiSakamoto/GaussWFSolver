#include <boost/tokenizer.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/format.hpp>

#include "SCFBase.hpp"
#include "HartreeFock.hpp"
#include "InputParser.hpp"
#include <iostream>
#include <cstdio>


static const double factor_Angstrom2Bohr = 1.8897259885789;
static const double factor_Bohr2Angstrom = 1./1.8897259885789;

int main(int argc, char **argv)
{
    if (argc < 2) {
        std::cerr << "Input file is not specifies" << std::endl;
        throw;
    }
    std::string input_filename(argv[1]);

    wf_solver::InputParser input(input_filename);
    wf_solver::MolecularSystem mol;
    std::string title;
    std::string xyz_file =  input.data().at("system");
    read_xyz_file(xyz_file, mol, title);

    std::printf("%s\n", title.c_str());
    std::printf("#============================================================\n");
    for(wf_solver::InputParser::container_type::const_iterator it = input.data().begin(); it != input.data().end(); it++) {
        std::printf("%-10s =  %s\n", it->first.c_str(), it->second.c_str() );
    }
    std::printf("#============================================================\n");
    for(size_t i = 0; i < mol.size(); i++) {
        std::array<wf_solver::REAL,3> pos( mol.atom_position(i));
        std::printf("%2s \t%12.8f\t%12.8f\t%12.8f\n", 
                wf_solver::Elements[mol.atomic_number(i)].c_str() , 
                pos[0] * factor_Bohr2Angstrom, pos[1] * factor_Bohr2Angstrom, pos[2] * factor_Bohr2Angstrom);
    }
    std::printf("#============================================================\n");

    std::printf("  Entering Hartree Fock  \n");
    wf_solver::HartreeFock hf_scf(mol);
    hf_scf.compute();
    if (hf_scf.convergence() == true) {
        std::printf("Total Energy: %17.10f\n", hf_scf.energy() );
    }
    return 0;
}
