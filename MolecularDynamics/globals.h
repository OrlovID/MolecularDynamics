#pragma once

namespace constant
{
	extern const double Kb = 1.380649e-23;	//boltzmann, J/K
	extern const double aem = 1.66053906660e-27;	//kg
	extern const double Na = 6.0221413e23;	//avogadro
}

const int fcc_neighbours_pos[12][3] = { { 0, 1, 1 }, { 0, 1, -1 }, { 0, -1, 1 }, { 0, -1, -1 }, 
										{ 1, 0, 1 }, { 1, 0, -1 }, { 1, 1, 0 }, { 1, -1, 0 },
										{ -1, 0, 1 }, { -1, 0, -1 },{ -1, 1, 0 }, { -1, -1, 0 }, };
const int fcc_npos[27] = { -1, -1, -1, -1, 0, 1, -1, 2, 3,
						-1, 4, 5, 6, -1, -1, 7, -1, -1, 
						-1, 8, 9, 10, -1, -1, 11, -1, -1 };
