#include "SCFBase.hpp"
#include "HartreeFock.hpp"
#include <iostream>

static const double factor_Angstrom2Bohr = 1.8897259885789;

int main()
{
    wf_solver::MolecularSystem mol;
    //mol.add_atom(1, 0., 0., 0.);
    //mol.add_atom(1, 1.4, 0., 0.);

    //mol.add_atom(6,Angstrom2Bohr(0., 0., 0.) );
    //mol.add_atom(1,Angstrom2Bohr(  0.000000,    0.000000,    1.083010) );
    //mol.add_atom(1,Angstrom2Bohr(  0.000000,    1.021071,   -0.361003) );
    //mol.add_atom(1,Angstrom2Bohr(  0.884274,   -0.510536,   -0.361003) );
    //mol.add_atom(1,Angstrom2Bohr( -0.884274,   -0.510536,   -0.361003) );

    mol.add_atom(6,  0., 0., 0.);
    mol.add_atom(1,  0.000000*factor_Angstrom2Bohr,    0.000000*factor_Angstrom2Bohr,    1.083010*factor_Angstrom2Bohr) ;
    mol.add_atom(1,  0.000000*factor_Angstrom2Bohr,    1.021071*factor_Angstrom2Bohr,   -0.361003*factor_Angstrom2Bohr) ;
    mol.add_atom(1,  0.884274*factor_Angstrom2Bohr,   -0.510536*factor_Angstrom2Bohr,   -0.361003*factor_Angstrom2Bohr) ;
    mol.add_atom(1, -0.884274*factor_Angstrom2Bohr,   -0.510536*factor_Angstrom2Bohr,   -0.361003*factor_Angstrom2Bohr) ;
    wf_solver::HartreeFock hf_scf(mol);
    hf_scf.compute();
    return 0;
}
