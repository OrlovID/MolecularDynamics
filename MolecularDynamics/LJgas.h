#pragma once
#include "stdafx.h"

using namespace std;

//вычисления в еденицах Леннарда Джонса
//длина в sigma, sigma в А = 1e-10 [м] 
//энергия в epsilon * Kb, epsilon в К
//масса в нуклоны * аем ( == 1)
//тогда время: в еденицах t_LJ, где
//t_LJ = time_param секунд
//time_param = sigma*1e-10 * sqrt(mass*аем / epsilon*Kb)
struct Param {
	double epsilon;		//in K				//real energy = E * epsilon * Kb
	double sigma;		//in A = 1e-10 m	//real length = x * sigma * 1e-10
	double mass;		//in aem			//real mass = m * aem
	double time_param;	//in s				//real time = t * time_param
};
struct Lattice {
	double coord[3];
	int type; // 0 = unknown, 1 = cube, 2 = fcc, ...
	int *neighbours;
	int pos[3];
	int mol_number;
	int neighbours_amount;
	//double a; //rib length
};

class LJgas
{
public:
	LJgas();
	~LJgas();
	void copy(LJgas *new_gase);
	void load_conditions(string file_in_name);
	void set_std_conditions(int number, double length);
	void set_dt(double dt);
	void save_conditions(string file_out_name);
	void save_step_xyz(string path, int flag_vel, double radius);
	void set_param(Param par);
	void save_stat(string file_out_name);
	void force();
	void verlet_vel();
	double get_mol_coord(int index, int k);
	void generate_coord_node(int number_per_rib);
	double generate_fcc(double l, int number_per_rib, int vacancy);
	void generate_coord_rand(double min_dr);
	void generetae_rand_vel(double temperature);
	double kinetic_energy();
	void velocity_mass_center();
	void velocity_mass_center(double *out);
	void null_vmc();
	double temperature_energy();
	double potential_energy();
	double temperature(); //in K
	double pressure();
	double maxvell_velocity(string file_out_name);
	double nevazka_r(LJgas *pgase);
	double nevazka_v(LJgas *pgase);
	double dr_square();
	double dr_square(int count, int d_step); //not ready
	double drdv_square();
	double auto_corr_vel();
	void reset_init_xv();
	//with time averaging
	void auto_corr_visc(string path);
	double thermostat_andersen(double nu, double T);
	double thermostat_andersen(double nu, double T, double *outer_array);
	double thermostat_berendsen(double tau, double T0);
	double ** get_vac_coord();
	int check_node_transition(); //check that molecule went to another lattice node
	int get_trans_count();
	bool m_flag_grid;
private:
	Molecule *Arr;
	//double *m_force_ij;
	Param m_par;
	double m_length;	//rib of cell; in sigma (len_LJ)
	double m_dt;		//step time of calc in time_LJ
	int m_number;		//number of molecules
	int m_step;			//step number; m_time = m_step*m_dt
	double m_time;
	int m_flags;
	double m_init_density;
	double m_init_T;
	//верны ли значения для газа для данных координат:
	//0x01 = coordinates set;
	//0x02 = velocity set;
	//0x04 = force() calculated
	//0x08 = VelocityMassCenter[3] calculated
	//0x10 = Kinetic Energy calculated
	//0x20 = Temperature Kinetic Energy calculated
	//0x40 = Temperature calculated
	//0x80 = Potential Energy calculated
	double m_EnergyKinetic;
	double m_EnergyTemperature;
	double m_EnergyPotential;
	double m_VelocityMassCenter[3];
	double m_Temperature;	
	double m_Pressure;		
	double m_pressureVirial;
	Molecule *m_vacancies;	
	int *m_vac_number;		//numbers of deleted molecules
	int m_vn;				//amount of --//--
	
	int m_lat_per_rib;		//number of llattice cells per calc cell rib
	Lattice *m_lattice_nodes; 
	int m_lat_number;		//number of nodes, including vacansies, if not ideal: != m_number
	double m_lattice_rib;	//length of lattice rib
	int m_transitions_count;
};


