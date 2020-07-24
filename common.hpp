#pragma once
#include <array>
#include <cmath>
#include <string>
#include "Eigen/Core"

namespace wf_solver {

// types
typedef double REAL;
typedef std::array<REAL,3> Position;

typedef Eigen::VectorXd VectorXReal;
typedef Eigen::Vector3d Vector3Real;
typedef Eigen::Vector2d Vector2Real;
typedef Eigen::MatrixXd MatrixXReal;

//const std::string atom_table[] = {0,
//    "H",  "He", 
//    "Li", "Be", "B",  "C",  "N", "O", "F", "Ne",
//    "Na", "Mg", "Al", "Si", "P", "S", "Cl", "Ar", 
//    "K",  "Ca"
//};

}
