#include "stdafx.h"

using namespace std;

LJgas::LJgas()
{
	Arr = new Molecule[1000];
	m_vacancies = NULL;
	m_vac_number = NULL;
	//m_force_ij = new double[1000 * 1000];
	m_number = 1000;
	m_length = 10;
	m_init_density = 1.0;
	m_dt = 1e-3;
	m_step = 0;
	m_pressureVirial = 0;
	m_flags = 0;
	m_vacancies = nullptr;
	m_vac_number = nullptr;
	m_lattice_nodes = nullptr;
	m_transitions_count = 0;
}

LJgas::~LJgas()
{
	int i;
	delete[] Arr;
	if(m_vacancies != nullptr) delete[] m_vacancies;
	if (m_vac_number != nullptr) delete[] m_vac_number;
	if (m_lattice_nodes != nullptr) {
		for (i = 0; i < m_lat_number; i++) {
			delete[] m_lattice_nodes[i].neighbours;
		}
		delete[] m_lattice_nodes;
	}
	//delete[] m_force_ij;
}

void LJgas::copy(LJgas *new_gase) {
	int i, k;
	delete[] new_gase->Arr;
	*new_gase = *this;
	Molecule *new_Arr = new Molecule[m_number];
	for (i = 0; i < m_number; i++) {
		new_Arr[i] = Arr[i];
	}
	new_gase->Arr = new_Arr;
}

void LJgas::load_conditions(string file_in_name)
{
	int i, k;
	double buf[3];
	ifstream file_in(file_in_name);
	file_in >> m_number >> m_length;
	delete[] Arr;
	Arr = new Molecule[m_number];
	for (i = 0; i < m_number; i++) {
		for (k = 0; k < 3; k++) {
			file_in >> buf[k];
		}
		Arr[i].set_init_coord(buf[0], buf[1], buf[2]);
		for (k = 0; k < 3; k++) {
			file_in >> buf[k];
		}
		Arr[i].set_init_vel(buf[0], buf[1], buf[2]);
	}
	m_flags = 0x01 | 0x02;
	file_in.close();
	m_init_density = m_number / (m_length * m_length * m_length);
	m_init_T = temperature();
	m_step = 0;
	m_pressureVirial = 0;
}

void LJgas::set_std_conditions(int number, double length)
{
	delete[] Arr;
	Arr = new Molecule[number];
	//m_force_ij = new double[number * number];
	m_number = number;
	m_length = length;
	m_init_density = number / (length * length * length);
	m_flags = 0;
	m_step = 0;
	m_pressureVirial = 0;
}

void LJgas::set_dt(double dt) {
	m_dt = dt;
}

void LJgas::set_param(Param par)
{
	m_par.epsilon = par.epsilon;
	m_par.sigma = par.sigma;
	m_par.mass = par.mass;
	m_par.time_param = par.time_param;
}

void LJgas::save_stat(string file_out_name)
{
	std::ofstream file_stat(file_out_name);
	file_stat << "dt step time: " << m_dt << endl;
	file_stat << "number: " << m_number << endl;
	file_stat << "initial temperature: " << m_init_T << endl;
	file_stat << "initial density: " << m_init_density << endl;
	file_stat << "length of cell: " << m_length << endl;
	file_stat.close();
}

void LJgas::save_conditions(string file_out_name) {
	ofstream file_out(file_out_name);
	//file_out << "    initialization \n";
	//file_out << "number \t rib length \n";
	file_out << m_number << " \t" << m_length << endl << endl;
	//file_out << "    coordinates and velocity\n";
	int i, k;
	for (i = 0; i < m_number; i++) {
		for (k = 0; k < 3; k++) {
			file_out << Arr[i].x[k] << " \t";
		}
		for (k = 0; k < 3; k++) {
			file_out << Arr[i].v[k] << " \t";
		}
		file_out << endl;
	}
	file_out.close();
}

//flag_vel != 0  <=> save velocities
//radius == 0 <=> do not print
void LJgas::save_step_xyz(string path, int flag_vel, double radius)
{
	int len = 0;
	int buf = m_step;
	int i, k;
	char *num_step;
	do {
		buf /= 10;
		len++;
	} while (buf);
	num_step = new char[len + 1];
	buf = m_step;
	for (i = len - 1; i >= 0; i--) {
		num_step[i] = 48 + buf % 10;
		buf /= 10;
	}
	num_step[len] = 0;
	string file_step(num_step);
	file_step += ".XYZ";
	//this->save_conditions(path + file_step);
	ofstream file_out(path + file_step);
	file_out << m_number << endl << endl;
	for (i = 0; i < m_number; i++) {
		for (k = 0; k < 3; k++) {
			file_out << Arr[i].x[k] << " \t";
		}
		if (flag_vel) {
			for (k = 0; k < 3; k++) {
				file_out << Arr[i].v[k] << " \t";
			}
		}
		if (radius) {
			file_out << radius << " \t";
		}
		file_out << endl;
	}
	file_out.close();
	delete[] num_step;
}

void LJgas::force()
{
	//also estimate part of pressure virial;
	double r2; // min (ri - rj)**2
	double buf[3]; //components of dr
	double r6; // min (ri - rj)**-6
	double force_buf;
	int i, j, k;
	m_pressureVirial = 0;
	for (i = 0; i < m_number; i++) {
		for (k = 0; k < 3; k++) {
			Arr[i].f[k] = 0;
		}
	}
	for (i = 0; i < m_number; i++) {
		for (j = i + 1; j < m_number; j++) {
			r2 = 0;
			for (k = 0; k < 3; k++) {
				buf[k] = Arr[i].x[k] - Arr[j].x[k];
				buf[k] -= m_length * rint(buf[k] / m_length);
				/*
				if (buf[k] > m_length / 2)
					buf[k] -= m_length;
				else if (buf[k] < -m_length / 2)
					buf[k] += m_length;
				*/
				r2 += (buf[k] * buf[k]);
			}
			r6 = pow(r2, -3); //r^-6
			force_buf = 24 * (2 * r6 * r6 - r6);
			m_pressureVirial += force_buf;
			force_buf /= r2;
			for (k = 0; k < 3; k++) {
				Arr[i].f[k] += (force_buf * buf[k]);
				Arr[j].f[k] -= (force_buf * buf[k]);	//3 law Newton
			}
		}
	}
	m_flags |= 0x04;
}

void LJgas::verlet_vel()
{
	//velocity verlet
	int i, k;
	int flbuf;
	double bufx;
	if ((m_flags & 0x04) == 0) {
		force();
	}
	for (i = 0; i < m_number; i++) {
		for (k = 0; k < 3; k++) {
			//v(t + dt/2)
			Arr[i].v[k] += 0.5 * Arr[i].f[k] * m_dt;
			//x(t+dt)
			bufx = Arr[i].v[k] * m_dt;
			Arr[i].x[k] += bufx;
			Arr[i].coord[k] += bufx;
			//check that molecule in cell, PBC
			flbuf = (int)(Arr[i].x[k] / m_length);
			Arr[i].flag[k] += flbuf;
			Arr[i].x[k] -= m_length * flbuf;
			if (Arr[i].x[k] > m_length) {
				Arr[i].x[k] -= m_length;
				Arr[i].flag[k]++;
			}
			else if (Arr[i].x[k] < 0) {
				Arr[i].x[k] += m_length;
				Arr[i].flag[k]--;
			}
		}
	}
	m_flags = 0x01 | 0x02;
	force();
	for (i = 0; i < m_number; i++) {
		for (k = 0; k < 3; k++) {
			//v(t+dt)
			Arr[i].v[k] += 0.5 * Arr[i].f[k] * m_dt;
		}
	}
	++m_step;
}

double LJgas::get_mol_coord(int index, int dim)
{
	return Arr[index].x[dim];
}

void LJgas::generate_coord_node(int number_per_rib)
{
	int i, j, k;
	double buf = m_length / number_per_rib;
	for (i = 0; i < number_per_rib; i++) {
		for (j = 0; j < number_per_rib; j++) {
			for (k = 0; k < number_per_rib; k++) {
				Arr[k + number_per_rib * j + number_per_rib * number_per_rib * i].set_init_coord(buf * i, buf * j, buf* k);
			}
		}
	}
	m_flags |= 0x01;
}

double LJgas::generate_fcc(double l, int number_per_rib, int vacancy)
{
	//vacancy == 0 no vacancy
	//vacancy == 1 with vacancy 
	int i, j, k, q, buf_number;
	int N0 = m_number;
	int N_fcc = 14; //number of atoms in fcc cell
	int fl_fccx[3] = { 1, 1, 0 };
	int fl_fccy[3] = { 1, 0, 1 };
	int fl_fccz[3] = { 0, 1, 1 };
	int pos[3] = { 0, 0, 0 };
	int ***grid;
	int buf;
	double Ei, Ef, Efv;
	m_lat_per_rib = number_per_rib;
	m_lat_number = m_number;
	m_flag_grid = true;
	m_lattice_nodes = new Lattice[m_number];
	m_lattice_rib = l;
	//set coords of atoms and lattice nodes
	i = 0;
	grid = new int**[2 * number_per_rib];
	for (pos[0] = 0; pos[0] < 2 * number_per_rib; pos[0]++) {
		grid[pos[0]] = new int*[2 * number_per_rib];
		for (pos[1] = 0; pos[1] < 2 * number_per_rib; pos[1]++) {
			grid[pos[0]][pos[1]] = new int[2 * number_per_rib];
			for (pos[2] = 0; pos[2] < 2 * number_per_rib; pos[2]++) {
				grid[pos[0]][pos[1]][pos[2]] = -2;
				if ((pos[0] + pos[1] + pos[2]) % 2 == 0) {
					Arr[i].set_init_coord(pos[0] * l / 2, pos[1] * l / 2, pos[2] * l / 2);
					for (k = 0; k < 3; k++) {
						m_lattice_nodes[i].coord[k] = Arr[i].x[k];
						m_lattice_nodes[i].pos[k] = pos[k];
					}
					m_lattice_nodes[i].mol_number = i;
					Arr[i].lattice_number = i;
					m_lattice_nodes[i].type = 2;
					m_lattice_nodes[i].neighbours_amount = 12;
					m_lattice_nodes[i].neighbours = new int[12];
					for (j = 0; j < 12; j++) {
						m_lattice_nodes[i].neighbours[j] = -1;
					}
					grid[pos[0]][pos[1]][pos[2]] = i;
					//coord of node #i
					i++;
				}
			}
		}
	}
	//find neighbours
	for (i = 0; i < m_lat_number; i++) {
		for (j = 0; j < 12; j++) {
			for (k = 0; k < 3; k++) {
				pos[k] = m_lattice_nodes[i].pos[k] + fcc_neighbours_pos[j][k];
				if (pos[k] >= 2 * number_per_rib) pos[k] -= 2 * number_per_rib;
				else if (pos[k] < 0) pos[k] += 2 * number_per_rib;
			}
			m_lattice_nodes[i].neighbours[j] = grid[pos[0]][pos[1]][pos[2]];
		}
	}
	for (i = 0; i < 2 * number_per_rib; i++) {
		for (j = 0; j < 2 * number_per_rib; j++) {
			delete[] grid[i][j];
		}
		delete[] grid[i];
	}
	delete[] grid;
	random_device rd;
	mt19937 mersenne(rd());
	/*
	for (i = 0; i < number_per_rib; i++) {
		for (j = 0; j < number_per_rib; j++) {
			for (k = 0; k < number_per_rib; k++) {
				Arr[k + number_per_rib * j + number_per_rib * number_per_rib * i].set_init_coord(l * i, l * j, l* k);
			}
		}
	}
	buf_number = pow(number_per_rib, 3);
	for (q = 0; q < 3; q++) {
		for (i = 0; i < number_per_rib; i++) {
			for (j = 0; j < number_per_rib; j++) {
				for (k = 0; k < number_per_rib; k++) {
					Arr[buf_number * (1 + q) + k + number_per_rib * j + number_per_rib * number_per_rib * i].set_init_coord(l * i + (l / 2) * fl_fccx[q], l * j + (l / 2) * fl_fccy[q], l * k + (l / 2) * fl_fccz[q]);
					//m_grid[0][buf_number * (1 + q) + k + number_per_rib * j + number_per_rib * number_per_rib * i] = l * i + (l / 2) * fl_fccx[q];
					//m_grid[1][buf_number * (1 + q) + k + number_per_rib * j + number_per_rib * number_per_rib * i] = l * j + (l / 2) * fl_fccy[q];
					//m_grid[2][buf_number * (1 + q) + k + number_per_rib * j + number_per_rib * number_per_rib * i] = l * k + (l / 2) * fl_fccz[q];
				}
			}
		}
	}
	*/
	m_flags |= 0x01;

	m_vn = vacancy;
	if (vacancy == 0) {
		m_vacancies = NULL;
		m_vac_number = NULL;
		return 0;
	}
	m_vac_number = new int[m_vn];
	m_vacancies = new Molecule[m_vn];
	Ei = this->potential_energy();
	Molecule* buf_arr = new Molecule[m_number];
	for (i = 0; i < m_number; i++) {
		buf_arr[i] = Arr[i];
	}
	delete[] Arr;
	m_number -= vacancy;
	m_init_density = m_number / (m_length * m_length * m_length);
	Arr = new Molecule[m_number];
	int vac_atom, flag_rand;
	for (j = 0; j < vacancy; j++) {
		do {
			vac_atom = mersenne() % m_number;
			//# of deleted atom
			flag_rand = 1;
			for (k = 0; k < j; k++) {
				if (vac_atom == m_vac_number[k]) flag_rand = 0;
			}
		} while (!flag_rand);
		m_vac_number[j] = vac_atom;
		m_vacancies[j] = buf_arr[vac_atom];
		m_lattice_nodes[vac_atom].mol_number = -1; //vacant (void) node
	}
	sort(m_vac_number, m_vac_number + vacancy);
	j = 0;
	for (i = 0; i < m_number; i++) {
		if (j != vacancy) {
			if (i + j == m_vac_number[j]) {
				j++;
			}
		}
		Arr[i] = buf_arr[i + j];
	}
	delete[] buf_arr;
	m_flags &= ~0x80;
	Ef = this->potential_energy();
	Efv = Ef - (N0 - vacancy) / (double)N0 * Ei;
	//Efv = Ef - Ei;
	Efv /= vacancy;
	return Efv;
}

void LJgas::generate_coord_rand(double min_dr)
{
	int i, j, k, num;
	double buf;
	double drij[3], drij0;
	random_device rd;
	mt19937 mersenne(rd());
	for (i = 0; i < m_number; i++) {
		for (k = 0; k < 3; k++) {
			buf = (double)mersenne() / pow(2.0, 32) * m_length;
			Arr[i].x[k] = buf;
		}
	}
	do {
		num = 0;
		for (i = 0; i < m_number; i++) {
			for (j = i + 1; j < m_number; j++) {
				drij0 = 0;
				for (k = 0; k < 3; k++) {
					drij[k] = (Arr[i].x[k] - Arr[j].x[k]);
					drij[k] -= m_length * rint(drij[k] / m_length);
					drij0 += drij[k] * drij[k];
				}
				drij0 = sqrt(drij0);
				if (drij0 < min_dr) {
					num++;
					buf = (min_dr - drij0) / drij0 / 2 + 0.05;
					for (k = 0; k < 3; k++) {
						Arr[i].x[k] += drij[k] * buf;
						Arr[j].x[k] -= drij[k] * buf;
						if (Arr[i].x[k] < 0) Arr[i].x[k] += m_length;
						else if (Arr[i].x[k] > m_length) Arr[i].x[k] -= m_length;
						if (Arr[j].x[k] < 0) Arr[j].x[k] += m_length;
						else if (Arr[j].x[k] > m_length) Arr[j].x[k] -= m_length;
					}
				}
			}
		}
	} while (num);
	for (i = 0; i < m_number; i++) {
		for (k = 0; k < 3; k++) {
			Arr[i].coord[k] = Arr[i].init_coord[k] = Arr[i].x[k];
		}
	}
	m_flags |= 0x01;
}

//maxvell distribution
void LJgas::generetae_rand_vel(double temperature)
{
	int i, k;
	for (k = 0; k < 3; k++) {
		std::random_device rd{};
		std::mt19937 gen{ rd() };
		std::normal_distribution<double> distribution(0.0, sqrt(temperature));
		for (i = 0; i < m_number; i++) {
			Arr[i].init_vel[k] = Arr[i].v[k] = distribution(gen);
		}
	}
	//- mass center velocity
	this->null_vmc();
	for (i = 0; i < m_number; i++) {
		Arr[i].energy();
	}
	m_init_T = temperature;
	m_flags |= 0x02;
}

double LJgas::kinetic_energy()
{
	if ((m_flags & 0x10) == 0) {
		int i;
		m_EnergyKinetic = 0;
		for (i = 0; i < m_number; i++) {
			m_EnergyKinetic += Arr[i].energy();
		}
		m_flags |= 0x10;
	}
	return m_EnergyKinetic;
}

void LJgas::velocity_mass_center()
{
	int i, k;
	for (k = 0; k < 3; k++) {
		m_VelocityMassCenter[k] = 0;
	}
	for (i = 0; i < m_number; i++) {
		for (k = 0; k < 3; k++) {
			m_VelocityMassCenter[k] += Arr[i].v[k];
		}
	}
	for (k = 0; k < 3; k++) {
		m_VelocityMassCenter[k] /= m_number;
	}
}

void LJgas::velocity_mass_center(double * out)
{
	int i, k;
	for (k = 0; k < 3; k++) {
		m_VelocityMassCenter[k] = 0;
	}
	for (i = 0; i < m_number; i++) {
		for (k = 0; k < 3; k++) {
			m_VelocityMassCenter[k] += Arr[i].v[k];
		}
	}
	for (k = 0; k < 3; k++) {
		m_VelocityMassCenter[k] /= m_number;
	}
	for (k = 0; k < 3; k++) {
		out[k] = m_VelocityMassCenter[k];
	}
}

void LJgas::null_vmc()
{
	int i, k;
	this->velocity_mass_center();
	for (i = 0; i < m_number; i++) {
		for (k = 0; k < 3; k++) {
			Arr[i].v[k] -= m_VelocityMassCenter[k];
			Arr[i].init_vel[k] = Arr[i].v[k];
		}
	}
	m_flags = 0x01 | 0x02;
}

double LJgas::temperature_energy()
{
	if ((m_flags & 0x20) == 0) {
		int i, k;
		m_EnergyTemperature = 0;
		velocity_mass_center();
		for (i = 0; i < m_number; i++) {
			for (k = 0; k < 3; k++) {
				m_EnergyTemperature += 0.5 * Arr[i].m * (Arr[i].v[k] - m_VelocityMassCenter[k]) * (Arr[i].v[k] - m_VelocityMassCenter[k]);
			}
		}
		m_flags |= 0x20;
	}
	return m_EnergyTemperature;
}

double LJgas::potential_energy()
{
	if ((m_flags & 0x80) == 0) {
		int i, j, k;
		double r2, buf_r2;
		m_EnergyPotential = 0;
		for (i = 0; i < m_number; i++) {
			for (j = i + 1; j < m_number; j++) {
				r2 = 0;
				for (k = 0; k < 3; k++) {
					buf_r2 = Arr[i].x[k] - Arr[j].x[k];
					if (buf_r2 > m_length / 2)
						buf_r2 -= m_length;
					else if (buf_r2 < -m_length / 2)
						buf_r2 += m_length;
					r2 += (buf_r2 * buf_r2);
				}
				buf_r2 = pow(r2, -3); //r^-6
				m_EnergyPotential += 4 * (buf_r2 * buf_r2 - buf_r2); // LJ Energy
			}
		}
		m_flags |= 0x80;
	}
	return m_EnergyPotential;
}

double LJgas::temperature()
{
	if ((m_flags & 0x40) == 0) {
		m_Temperature = 2.0 / 3.0 * kinetic_energy() / m_number;
		m_flags |= 0x40;
	}
	return m_Temperature;
}

double LJgas::pressure()
{
	double V = pow(m_length, 3);
	if ((m_flags & 0x04) == 0) {
		force();
	}
	m_Pressure = (m_number * temperature() + m_pressureVirial / 3) / V;
	return m_Pressure;
}

//todo
double LJgas::maxvell_velocity(string file_out_name)
{
	int i;
	for (i = 0; i < m_number; i++) {

	}
	return 0.0;
}

double LJgas::nevazka_r(LJgas *pgase)
{
	int i, k;
	double res = 0;
	for (i = 0; i < m_number; i++) {
		for (k = 0; k < 3; k++) {
			res += (Arr[i].x[k] - pgase->Arr[i].x[k]) * (Arr[i].x[k] - pgase->Arr[i].x[k]);
		}
	}
	res /= m_number;
	return res;
}

double LJgas::nevazka_v(LJgas *pgase)
{
	int i, k;
	double res = 0;
	for (i = 0; i < m_number; i++) {
		for (k = 0; k < 3; k++) {
			res += (Arr[i].v[k] - pgase->Arr[i].v[k]) * (Arr[i].v[k] - pgase->Arr[i].v[k]);
		}
	}
	res /= m_number;
	return res;
}

double LJgas::dr_square()
{
	int i, k;
	double res = 0;
	for (i = 0; i < m_number; i++) {
		for (k = 0; k < 3; k++) {
			res += (Arr[i].coord[k] - Arr[i].init_coord[k]) * (Arr[i].coord[k] - Arr[i].init_coord[k]);
		}
	}
	return res / m_number;
}

double LJgas::dr_square(int count, int d_step)
{
	int i, k;
	double res = 0;
	for (i = 0; i < m_number; i++) {
		for (k = 0; k < 3; k++) {
			res += (Arr[i].coord[k]  - Arr[i].init_coord[k]) * (Arr[i].coord[k] - Arr[i].init_coord[k]);
		}
	}
	return res / m_number;
}

double LJgas::drdv_square()
{
	int i, k;
	double res = 0;
	double V = m_length * m_length * m_length;
	for (i = 0; i < m_number; i++) {
		res += (Arr[i].coord[0] * Arr[i].v[1] - Arr[i].init_coord[0] * Arr[i].init_vel[1]) * Arr[i].m;
	}
	return (res * res / this->temperature() / V / 2);
}

//useful for diffusion
double LJgas::auto_corr_vel()
{
	int i, k;
	double res = 0;
	for (i = 0; i < m_number; i++) {
		for (k = 0; k < 3; k++) {
			res += Arr[i].v[k] * Arr[i].init_vel[k];
		}
	}
	return res / m_number;
}

void LJgas::reset_init_xv()
{
	int i, k;
	for (i = 0; i < m_number; i++) {
		for (k = 0; k < 3; k++) {
			Arr[i].init_coord[k] = Arr[i].coord[k] = Arr[i].x[k];
			Arr[i].init_vel[k] = Arr[i].v[k];
			Arr[i].flag[k] = 3;
		}
	}
}

void LJgas::auto_corr_visc(string path)
{
	double sigma = 0; 
	int i, j;
	//j = : 0 - xy; 1 - xz; 2 - yz
	int k1[3] = { 0, 0, 1 };
	int k2[3] = { 1, 2, 2 };
	double V = m_length * m_length * m_length;
	ofstream file_out(path, ios::app);
	file_out << m_step << "\t";
	for (j = 0; j < 3; j++) {
		sigma = 0;
		for (i = 0; i < m_number; i++) {
			sigma += Arr[i].m * Arr[i].v[k1[j]] * Arr[i].v[k2[j]] + Arr[i].x[k1[j]] * Arr[i].f[k2[j]];
		}
		sigma /= V;
		file_out << sigma << "\t";
	}
	file_out << endl;
	file_out.close();
}

double LJgas::thermostat_andersen(double nu, double T)
{
	double prob = nu * m_dt;
	int i, k, num = 0;
	std::random_device rd{};
	std::mt19937 mersenne{ rd() };
	for (i = 0; i < m_number; i++) {
		if ((double)mersenne() / pow(2.0, 32) < prob) {
			for (k = 0; k < 3; k++) {
				std::random_device rd2{};
				std::mt19937 gen{ rd2() };
				std::normal_distribution<double> distribution(0.0, sqrt(T));
				Arr[i].v[k] = distribution(gen);
			}
			num++;
		}
	}
	m_flags = 0x01 | 0x02;
	//- mass center velocity
	this->null_vmc();
	//
	for (i = 0; i < m_number; i++) {
		Arr[i].energy();
	}
	//
	return num;
}

double LJgas::thermostat_andersen(double nu, double T, double *outer_array)
{
	double prob = nu * m_dt;
	int i, k, num = 0;
	std::random_device rd{};
	std::mt19937 mersenne{ rd() };
	for (i = 0; i < m_number; i++) {
		if ((double)mersenne() / pow(2.0, 32) < prob) {
			for (k = 0; k < 3; k++) {
				std::random_device rd2{};
				std::mt19937 gen{ rd2() };
				std::normal_distribution<double> distribution(0.0, sqrt(T));
				Arr[i].v[k] = distribution(gen);
			}
			outer_array[num] = i;
			num++;
		}
	}
	m_flags = 0x01 | 0x02;
	//- mass center velocity
	this->null_vmc();
	//
	for (i = 0; i < m_number; i++) {
		Arr[i].energy();
	}
	//
	return num;
}

double LJgas::thermostat_berendsen(double tau, double T0) {
	int i, k;
	double lambda;
	this->temperature();
	lambda = sqrt(1 + m_dt / tau * (T0 / m_Temperature - 1));
	for (i = 0; i < m_number; i++) {
		for (k = 0; k < 3; k++) {
			Arr[i].v[k] *= lambda;
		}
	}
	this->null_vmc();
	return lambda;
}

double ** LJgas::get_vac_coord()
{
	int i, j, k;
	double **vac_coord = new double*[m_vn];
	for (i = 0; i < m_vn; i++) {
		vac_coord[i] = new double[3];
		for (k = 0; k < 3; k++) {
			vac_coord[i][k] = m_vacancies[i].x[k];
		}
	}
	return vac_coord;
}

//fcc lattice only
int LJgas::check_node_transition()
{
	int i, j, k;
	int count = 0;
	int shift[3] = { 0, 0, 0 };
	int a[3] = { 9, 3, 1 };
	int buf, buf2;
	double dr[3];
	//int *new_lattice = new int[m_lat_number];
	for (i = 0; i < m_number; i++) {
		buf = 0;
		for (k = 0; k < 3; k++) {
			dr[k] = Arr[i].x[k] - m_lattice_nodes[Arr[i].lattice_number].coord[k];
			if (dr[k] > m_length / 2) { 
				//m_lattice_rib < m_length /2 !
				dr[k] -= m_length;
			}
			else if (dr[k] < -m_length / 2) {
				dr[k] += m_length;
			}
			if (dr[k] > m_lattice_rib / 2) shift[k] = 1;
			else if (dr[k] < - m_lattice_rib / 2) shift[k] = 2;
			else shift[k] = 0;
			buf += shift[k] * a[k];
		}
		buf2 = fcc_npos[buf];
		if (buf2 != -1) {
			Arr[i].lattice_number = m_lattice_nodes[Arr[i].lattice_number].neighbours[buf2];
			count++;
		}
	}
	return count;
}

int LJgas::get_trans_count()
{
	return m_transitions_count;
}


