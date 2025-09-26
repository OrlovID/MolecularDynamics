#include "stdafx.h"

Molecule::Molecule()
{
	reset_force();
	set_init_coord(0, 0, 0);
	set_init_vel(0, 0, 0);
	m = 1;
}


Molecule::~Molecule()
{
}

void Molecule::set_init_coord(double x1, double x2, double x3)
{
	x[0] = coord[0] = init_coord[0] = x1;
	x[1] = coord[1] = init_coord[1] = x2;
	x[2] = coord[2] = init_coord[2] = x3;
	flag[0] = flag[1] = flag[2] = 0;
}

void Molecule::set_init_vel(double v1, double v2, double v3)
{
	init_vel[0] = v[0] = v1;
	init_vel[1] = v[1] = v2; 
	init_vel[2] = v[2] = v3;
}

void Molecule::set_mass(double mass)
{
	m = mass;
}

double Molecule::velocity()
{
	double buf = 0;
	int k;
	for (k = 0; k < 3; k++) {
		buf += (v[k] * v[k]);
	}
	buf = sqrt(buf);
	v_abs = buf;
	return buf;
}

double Molecule::energy()
{
	double buf;
	velocity();
	buf = 0.5 * m * v_abs * v_abs;
	E = buf;
	return buf;
}

void Molecule::reset_force()
{
	int k;
	for (k = 0; k < 3; k++) {
		f[k] = 0;
	}
}