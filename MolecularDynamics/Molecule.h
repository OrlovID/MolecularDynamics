#pragma once
class LJgas;

class Molecule
{
public:
	Molecule(); 
	~Molecule();
	void set_init_coord(double x1, double x2, double x3);
	void set_init_vel(double v1, double v2, double v3);
	void set_mass(double mass);
	double velocity();
	double energy();
	void reset_force();
	friend class LJgas;

	double x[3];	//coordinates with periodic border conditions //time: t, t + dt, t + 2dt, ...
	double coord[3]; //real coordinates
	double init_coord[3]; //initial coordinates
	double init_vel[3]; //initial velocities
	double v[3];	//velocity //time: t - dt/2, t + dt/2, t + 3dt/2, ...
	double f[3]; //force //time: t, t + dt, t + 2dt, ...
	double m; //mass
	double v_abs; //absolut velocity
	double E; //Energy
	int flag[3]; //number of cell were the real molecule flag[k] = (int) (coord[k]/length)
	int lattice_number;
private:
};

