#pragma once
//#include <vector>
//#include <stddef.h>
#include <complex>
#include <fstream>
#include <string>
#include "material.h"
#include "matrix.h"
#include "geometry.h"

class global_parameters {
private:
	double PI = 3.14159265358979323846;
	class geometry_parameters* pt_geo = &(geometry_parameters::geo);

	unsigned int* pt_material_layer_unit;
	unsigned int* pt_material_layer;
	long int nx, ny, nz, n;
	double dx, dy, dz;
private:
	void set_global();
	void read_materials();
	void read_struct();
	void check_material(material& mat);
	void get_cijkl(material& mat);
	void read_prescribe_Eext();
	void read_prescribe_Hext();
	void read_Hext_nonunif();
public:
	unsigned long long int step;
	double dt;
	double time_device = 0.;
	long prescribe_index_E = 0;
	long prescribe_index_H = 0;
	unsigned int num_materials;
	unsigned long int NFE, NFM, NAFM;
	unsigned long long* magcell_index;

	double theta_mag, phi_mag;
	double px_init, py_init, pz_init;

	//External magnetic field
	double Hext[3];
	double Hext_stat[3];
	double Hext_altr[3];
	double H_altr_freq[3];
	bool if_read_Hext_nonunif;
	matrix3d<double> Hextx_nonunif;
	matrix3d<double> Hexty_nonunif;
	matrix3d<double> Hextz_nonunif;

	bool if_prescribe_Hext;
	long num_prescribe_Hext;
	double* prescribe_Hext_time;
	double* prescribe_Hx;
	double* prescribe_Hy;
	double* prescribe_Hz;

	//External electric field
	double Eext[3];
	double Eext_stat[3];
	double Eext_altr[3];
	double E_altr_freq[3];

	bool if_prescribe_Eext;
	long num_prescribe_Eext;
	double* prescribe_Eext_time;
	double* prescribe_Ex;
	double* prescribe_Ey;
	double* prescribe_Ez;

	bool if_magnetostatics, if_mag_1Dmodel, if_elec_1Dmodel;
	bool if_demag_factor;
	bool if_elec_1D_compensate; // , if_elec_1D_dynamical_compensate;
	bool if_EM_backaction_on_P;
	bool if_EM_backaction_on_M;
	bool if_EMdynamic;
	bool if_input_planeEM_E;
	double planeEM_E[3];
	double planeEM_E_freq[3];
	bool if_EM_fromP;
	bool if_EM_fromM;
	bool if_elastodynamics;
	bool if_elasto_backaction_from_pORm;
	bool if_elasto_on_pORm;

	bool if_OrientRotate;
	double Tg2c[3][3], Tc2g[3][3];
	double demag_fac_x, demag_fac_y, demag_fac_z;

	material* material_parameters;
	matrix3d<unsigned int> material_cell;

	bool if_elastostatic, if_elastostatic_1D, if_elastostatic_film;

	long free_charge_nzi, free_charge_nzf;
	double free_charge;//, total_free_charge;
	double* free_charge_surfaces;
	double* free_charge_cells;

	//bool if_open_circuit_BC; //OC
	//long OC_nzi, OC_nzf;

	bool if_read_strain_ext;
	long strain_ext_nzi, strain_ext_nzf;
	double exx_ext_glb, eyy_ext_glb, exy_ext_glb;

	double scale_elec, scale_elasto;
	unsigned int elec_solver_limit, elasto_solver_limit;
	double elec_solver_tol, elasto_solver_tol;

	bool if_gaussian_strain_pulse;
	unsigned int input_strain_type;
	unsigned int input_strain_sourcez;
	char input_strain_component;
	double sigma_gauss, amplitude_gauss;
	unsigned long num_strain_cycles;

	bool if_FM_all, if_FE_all, if_AFM_all;
	bool if_PEC_all;
	bool if_periodic_allsurface;

	bool if_1D_ABC, if_1D_ABC_ELAST_onlytop, if_1D_ABC_EM_onlytop, if_PEC_XY, if_PEC_Z;
	
	//bool if_elastic_Liaoabsorb;

	double weighting_factor_fourth_Liao;

	bool if_flexo;

	bool if_spin_pumping;
	bool* if_J_ISHE;
	long int* source_z;
	double* thickness_pumping_layer;
	
	bool if_output_ave;
	unsigned long long int output_step_ave;

	bool if_output_uandv;
	unsigned long long int output_step_uandv;
	
	bool if_output_m, if_output_AFMm, if_output_p, if_output_q, if_output_em, if_output_strain;
	unsigned long long int output_step_m, output_step_AFMm, output_step_p, output_step_q, output_step_em, output_step_strain;
	bool if_output_only_magcell;
	unsigned long long int output_step_magcell;
	bool if_output_em_onecell;
	unsigned long long int em_onecell_index, output_step_em_onecell;

	bool if_output_Hstat, if_output_Estat, if_output_eigenstraint0_crt;
	unsigned long long int output_step_Hstat, output_step_Estat, output_step_eigenstraint0_crt;

	bool if_output_elastoforce;
	unsigned long long int output_step_elastoforce;

	bool if_output_emYee, if_output_Jp, if_output_Jishe;
	unsigned long long int output_step_emYee, output_step_Jp, output_step_Jishe;

	bool if_input_m, if_input_AFMm, if_input_pandq, if_input_em, if_input_uandv;
	bool if_input_Hstat, if_input_Estat, if_input_straint0, if_input_eigenstraint0_crt;
	bool if_input_Jp, if_input_elastoforce;

	bool if_Jf_input;
	unsigned int Jf_input_type;
	unsigned long Jfin_xi, Jfin_yi, Jfin_zi;
	unsigned long Jfin_xf, Jfin_yf, Jfin_zf;
	char Jf_input_component;
	double Jf_input_amp, Jf_input_freq, Jf_input_sigma;
	unsigned long num_Jf_cycles;

	bool if_prescribed_m;
	double precess_angle, precess_frequency;

	double dkx, dky, dkz;
	double* kx, * ky, * kz;

	matrix3d<double> k_norm;

	//matrix3d<std::complex<double>> test;

	static global_parameters glb;

public:
	//global_parameters(geometry_parameters*);
	void read_global();
	void copy_to_device();
	void update_to_device();

#pragma acc routine seq nohost
	void time_marching_device_serial();
};

////---------------CONSTANTS----------------------//
//constexpr double PI = 3.14159265358979323846;
////----------------------------------------------//
