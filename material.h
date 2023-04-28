#pragma once

class material {
public:
	/*elastic parameters*/
	double c11, c12, c44;
	double density;
	double elast_mass_damping, elast_stiff_damping;
	double cijkl_glb[6][6];

	/*magnetic parameters*/
	bool if_FM;
	double Ms, gyro, FM_damping;
	unsigned int anisotropy_type;
	double K1, K2;
	double Aex;
	double lamda100, lamda111, B1, B2;
	double iDMI;

	bool if_spin_pump;
	double spin_hall_angle, spin_mix_cond, spin_diffus_length;

	/*antiferromagnetic parameters (2 sublattices)*/
	bool if_AFM;
	double Ms_AFM1, Ms_AFM2;
	double gyro_AFM1, gyro_AFM2;
	double damping_AFM1, damping_AFM2;
	double K1_AFM1, K1_AFM2;
	double uniaxis_AFM[3];
	double Aex_AFM1, Aex_AFM2;
	double J_AFM;
	double lamda100_AFM1, lamda100_AFM2;
	double lamda111_AFM1, lamda111_AFM2;
	double B1_AFM1, B1_AFM2;
	double B2_AFM1, B2_AFM2;

	//bool if_spin_pump_AFM;
	//double spin_hall_angle_AFM, spin_mix_cond_AFM, spin_diffus_length_AFM;

	/*Ferroelectric */
	bool if_FE;
	double FE_mass, FE_damping;
	double a1, a11, a12, a111, a112, \
		a123, a1111, a1112, a1122, a1123;
	double G11;
	double Q11, Q12, Q44;
	double tv1, tv2, tv3, tv4, tv5, tv6, tv7, tv8;
	double f11, f12, f44;
	double F11, F12, F44;

	/*electromagnetic parameters*/
	bool if_PEC;
	double r_permittivity;
	double conductivity;
	double comp_n1, omega_plasma_n1, tao_e_n1;
	double comp_n2, omega_plasma_n2, tao_e_n2;
	double comp_n3, omega_plasma_n3, tao_e_n3;
};
