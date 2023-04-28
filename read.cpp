#include "global.h"
#include "geometry.h"

void global_parameters::read_global() {
	long int num_layer_unit = geometry_parameters::geo.num_layer_unit;
	std::string aline;

	std::ifstream file("system_setting.in");
	if (file.is_open()) {
		for (long int i = 0; i < 8;) {
			std::getline(file, aline);
			if (aline.length() > 0) {
				i = i + 1;
			}
		}
		file >> step;			file.ignore(1000, '\n');
		file >> dt;				file.ignore(1000, '\n');
		file >> num_materials;	file.ignore(1000, '\n');

		pt_material_layer_unit = new unsigned int[num_layer_unit];
		for (long int i = 0; i < num_layer_unit; i++) {
			file >> pt_material_layer_unit[i];
		}
		file.ignore(1000, '\n');

		file >> std::boolalpha >> if_OrientRotate;		file.ignore(1000, '\n');
		file >> Tg2c[0][0] >> Tg2c[0][1] >> Tg2c[0][2]; file.ignore(1000, '\n');
		file >> Tg2c[1][0] >> Tg2c[1][1] >> Tg2c[1][2]; file.ignore(1000, '\n');
		file >> Tg2c[2][0] >> Tg2c[2][1] >> Tg2c[2][2]; file.ignore(1000, '\n');

		file >> theta_mag >> phi_mag; file.ignore(1000, '\n');
		file >> px_init >> py_init >> pz_init; file.ignore(1000, '\n');

		file >> std::boolalpha >> if_read_Hext_nonunif;											file.ignore(1000, '\n');
		file >> Hext_stat[0] >> Hext_stat[1] >> Hext_stat[2];									file.ignore(1000, '\n');
		file >> Hext_altr[0] >> Hext_altr[1] >> Hext_altr[2];									file.ignore(1000, '\n');
		file >> H_altr_freq[0] >> H_altr_freq[1] >> H_altr_freq[2];								file.ignore(1000, '\n');
		file >> std::boolalpha >> if_prescribe_Hext >> num_prescribe_Hext;						file.ignore(1000, '\n');

		file >> Eext_stat[0] >> Eext_stat[1] >> Eext_stat[2];									file.ignore(1000, '\n');
		file >> Eext_altr[0] >> Eext_altr[1] >> Eext_altr[2];									file.ignore(1000, '\n');
		file >> E_altr_freq[0] >> E_altr_freq[1] >> E_altr_freq[2];								file.ignore(1000, '\n');
		file >> std::boolalpha >> if_prescribe_Eext >> num_prescribe_Eext;						file.ignore(1000, '\n');

		file >> std::boolalpha >> if_magnetostatics >> if_mag_1Dmodel;							file.ignore(1000, '\n');
		file >> std::boolalpha >> if_demag_factor >> demag_fac_x >> demag_fac_y >> demag_fac_z; file.ignore(1000, '\n');
		file >> std::boolalpha >> if_elec_1Dmodel >> if_elec_1D_compensate;						file.ignore(1000, '\n');
		file >> std::boolalpha >> if_EMdynamic >> if_input_planeEM_E;;							file.ignore(1000, '\n');
		file >> planeEM_E[0] >> planeEM_E[1] >> planeEM_E[2];									file.ignore(1000, '\n');
		file >> planeEM_E_freq[0] >> planeEM_E_freq[1] >> planeEM_E_freq[2];					file.ignore(1000, '\n');
		//file >> std::boolalpha >> if_elec_1D_compensate >> if_elec_1D_dynamical_compensate		file.ignore(1000, '\n');
		file >> std::boolalpha >> if_EM_fromP >> if_EM_fromM;									file.ignore(1000, '\n');
		file >> std::boolalpha >> if_EM_backaction_on_P >> if_EM_backaction_on_M;				file.ignore(1000, '\n');
		file >> free_charge_nzi >> free_charge_nzf >> free_charge;								file.ignore(1000, '\n');
		//file >> std::boolalpha >> if_open_circuit_BC >> OC_nzi >> OC_nzf;						file.ignore(1000, '\n');
		file >> scale_elec >> elec_solver_limit >> elec_solver_tol;								file.ignore(1000, '\n');

		file >> std::boolalpha >> if_elastodynamics >> if_elasto_backaction_from_pORm >> if_elasto_on_pORm;		file.ignore(1000, '\n');
		file >> std::boolalpha >> if_gaussian_strain_pulse >> input_strain_type >> input_strain_sourcez;		file.ignore(1000, '\n');
		file >> input_strain_component >> sigma_gauss >> amplitude_gauss >> num_strain_cycles;					file.ignore(1000, '\n');
		file >> std::boolalpha >> if_elastostatic >> if_elastostatic_1D;										file.ignore(1000, '\n');
		file >> std::boolalpha >> if_read_strain_ext >> if_elastostatic_film;									file.ignore(1000, '\n');
		file >> strain_ext_nzi >> strain_ext_nzf >> exx_ext_glb >> eyy_ext_glb >> exy_ext_glb;					file.ignore(1000, '\n');
		file >> scale_elasto >> elasto_solver_limit >> elasto_solver_tol;										file.ignore(1000, '\n');

		file >> std::boolalpha >> if_1D_ABC >> if_1D_ABC_ELAST_onlytop >> if_1D_ABC_EM_onlytop >> weighting_factor_fourth_Liao; file.ignore(1000, '\n');
		file >> std::boolalpha >> if_PEC_XY >> if_PEC_Z;										  file.ignore(1000, '\n');

		file >> std::boolalpha >> if_flexo; file.ignore(1000, '\n');
		file >> std::boolalpha >> if_spin_pumping; file.ignore(1000, '\n');

		file >> std::boolalpha >> if_input_m >> if_input_AFMm >> if_input_Hstat;						file.ignore(1000, '\n');
		file >> std::boolalpha >> if_input_pandq >> if_input_Estat;										file.ignore(1000, '\n');
		file >> std::boolalpha >> if_input_em >> if_input_Jp;											file.ignore(1000, '\n');
		file >> std::boolalpha >> if_input_uandv >> if_input_straint0 >> if_input_eigenstraint0_crt;	file.ignore(1000, '\n');
		file >> std::boolalpha >> if_input_elastoforce;													file.ignore(1000, '\n');

		file >> std::boolalpha >> if_output_ave >> output_step_ave;			file.ignore(1000, '\n');

		file >> std::boolalpha >> if_output_m >> output_step_m >> if_output_only_magcell >> output_step_magcell;	file.ignore(1000, '\n');
		file >> std::boolalpha >> if_output_AFMm >> output_step_AFMm;		file.ignore(1000, '\n');
		file >> std::boolalpha >> if_output_p >> output_step_p;				file.ignore(1000, '\n');
		file >> std::boolalpha >> if_output_em >> output_step_em >> if_output_em_onecell >> output_step_em_onecell >> em_onecell_index;	file.ignore(1000, '\n');
		file >> std::boolalpha >> if_output_strain >> output_step_strain;	file.ignore(1000, '\n');

		file >> std::boolalpha >> if_output_Hstat >> output_step_Hstat;							file.ignore(1000, '\n');
		file >> std::boolalpha >> if_output_Estat >> output_step_Estat;							file.ignore(1000, '\n');
		//file >> std::boolalpha >> if_output_straint0 >> output_step_straint0;					file.ignore(1000, '\n');
		file >> std::boolalpha >> if_output_eigenstraint0_crt >> output_step_eigenstraint0_crt;	file.ignore(1000, '\n');
		file >> std::boolalpha >> if_output_q >> output_step_q;									file.ignore(1000, '\n');
		file >> std::boolalpha >> if_output_uandv >> output_step_uandv;							file.ignore(1000, '\n');
		file >> std::boolalpha >> if_output_elastoforce >> output_step_elastoforce;				file.ignore(1000, '\n');
		file >> std::boolalpha >> if_output_emYee >> output_step_emYee;							file.ignore(1000, '\n');
		file >> std::boolalpha >> if_output_Jp >> output_step_Jp;								file.ignore(1000, '\n');
		file >> std::boolalpha >> if_output_Jishe >> output_step_Jishe;							file.ignore(1000, '\n');

		file >> std::boolalpha >> if_Jf_input;													file.ignore(1000, '\n');
		file >> Jfin_xi >> Jfin_yi >> Jfin_zi;													file.ignore(1000, '\n');
		file >> Jfin_xf >> Jfin_yf >> Jfin_zf;													file.ignore(1000, '\n');
		file >> Jf_input_type >> Jf_input_component >> Jf_input_amp >> Jf_input_freq >> Jf_input_sigma >> num_Jf_cycles;			file.ignore(1000, '\n');
		file >> std::boolalpha >> if_prescribed_m >> precess_angle >> precess_frequency;		file.ignore(1000, '\n');
	}
	file.close();
	//----------Build transformation matrix-------------//
	Tc2g[0][0] = Tg2c[0][0]; Tc2g[0][1] = Tg2c[1][0]; Tc2g[0][2] = Tg2c[2][0];
	Tc2g[1][0] = Tg2c[0][1]; Tc2g[1][1] = Tg2c[1][1]; Tc2g[1][2] = Tg2c[2][1];
	Tc2g[2][0] = Tg2c[0][2]; Tc2g[2][1] = Tg2c[1][2]; Tc2g[2][2] = Tg2c[2][2];

	read_materials();
	set_global();
}

void global_parameters::read_materials() {
	material mat;
	std::string aline;
	double norm;

	material_parameters = new material[num_materials];

	std::ifstream file("materials.in");
	if (file.is_open()) {
		for (long int i = 0; i < num_materials; i++) {
			std::getline(file, aline);
			//while (aline.length() == 0)
			//	std::getline(file, aline);

			//Elastic parameters
			file >> mat.c11 >> mat.c12 >> mat.c44; file.ignore(1000, '\n');
			file >> mat.density; file.ignore(1000, '\n');
			file >> mat.elast_mass_damping >> mat.elast_stiff_damping; file.ignore(1000, '\n');

			//Ferromagnetic parameters
			file >> std::boolalpha>> mat.if_FM; file.ignore(1000, '\n');
			file >> mat.Ms >> mat.gyro >> mat.FM_damping; file.ignore(1000, '\n');
			file >> mat.anisotropy_type; file.ignore(1000, '\n');
			file >> mat.K1 >> mat.K2; file.ignore(1000, '\n');
			file >> mat.Aex; file.ignore(1000, '\n');
			file >> mat.iDMI; file.ignore(1000, '\n');
			file >> mat.lamda100 >> mat.lamda111; file.ignore(1000, '\n');

			mat.B1 = -3. / 2. * mat.lamda100 * (mat.c11 - mat.c12);
			mat.B2 = -3. * mat.lamda111 * mat.c44;

			//Anti-Ferromagnetic parameters
			file >> std::boolalpha >> mat.if_AFM;													file.ignore(1000, '\n');
			file >> mat.Ms_AFM1 >> mat.Ms_AFM2;														file.ignore(1000, '\n');
			file >> mat.gyro_AFM1 >> mat.gyro_AFM2;													file.ignore(1000, '\n');
			file >> mat.damping_AFM1 >> mat.damping_AFM2;											file.ignore(1000, '\n');
			file >> mat.K1_AFM1 >> mat.K1_AFM2;														file.ignore(1000, '\n');

			file >> mat.uniaxis_AFM[0] >> mat.uniaxis_AFM[1] >> mat.uniaxis_AFM[2];					file.ignore(1000, '\n');
			norm = sqrt(mat.uniaxis_AFM[0] * mat.uniaxis_AFM[0] + mat.uniaxis_AFM[1] * mat.uniaxis_AFM[1] + mat.uniaxis_AFM[2] * mat.uniaxis_AFM[2]);
			mat.uniaxis_AFM[0] = mat.uniaxis_AFM[0] / norm;
			mat.uniaxis_AFM[1] = mat.uniaxis_AFM[1] / norm;
			mat.uniaxis_AFM[2] = mat.uniaxis_AFM[2] / norm;

			file >> mat.Aex_AFM1 >> mat.Aex_AFM2;													file.ignore(1000, '\n');
			file >> mat.J_AFM;																		file.ignore(1000, '\n');
			file >> mat.lamda100_AFM1 >> mat.lamda100_AFM2;											file.ignore(1000, '\n');
			file >> mat.lamda111_AFM1 >> mat.lamda111_AFM2;											file.ignore(1000, '\n');

			mat.B1_AFM1 = -3. / 2. * mat.lamda100_AFM1 * (mat.c11 - mat.c12);
			mat.B1_AFM2 = -3. / 2. * mat.lamda100_AFM2 * (mat.c11 - mat.c12);
			mat.B2_AFM1 = -3. * mat.lamda111_AFM1 * mat.c44;
			mat.B2_AFM2 = -3. * mat.lamda111_AFM2 * mat.c44;

			file >> std::boolalpha >> mat.if_spin_pump; file.ignore(1000, '\n');
			file >> mat.spin_hall_angle >> mat.spin_mix_cond >> mat.spin_diffus_length; file.ignore(1000, '\n');

			//file >> std::boolalpha >> mat.if_spin_pump_AFM;											file.ignore(1000, '\n');
			//file >> mat.spin_hall_angle_AFM >> mat.spin_mix_cond_AFM >> mat.spin_diffus_length_AFM; file.ignore(1000, '\n');

			//Ferroelectric parameters
			file >> std::boolalpha >> mat.if_FE; file.ignore(1000, '\n');
			file >> mat.FE_mass >> mat.FE_damping; file.ignore(1000, '\n');
			file >> mat.a1; file.ignore(1000, '\n');
			file >> mat.a11 >> mat.a12; file.ignore(1000, '\n');
			file >> mat.a111 >> mat.a112 >> mat.a123; file.ignore(1000, '\n');
			file >> mat.a1111 >> mat.a1112 >> mat.a1122 >> mat.a1123; file.ignore(1000, '\n');
			file >> mat.G11; file.ignore(1000, '\n');
			file >> mat.Q11 >> mat.Q12 >> mat.Q44; file.ignore(1000, '\n');
			file >> mat.f11 >> mat.f12 >> mat.f44; file.ignore(1000, '\n');

			mat.F11 = ((mat.c11 + mat.c12) * mat.f11 - 2. * mat.c12 * mat.f12) / (mat.c11 + 2. * mat.c12) / (mat.c11 - mat.c12);
			mat.F12 = (mat.c11 * mat.f12 - mat.c12 * mat.f11) / (mat.c11 + 2. * mat.c12) / (mat.c11 - mat.c12);
			mat.F44 = mat.f44 / mat.c44;

			mat.tv1 = 2. * (mat.c11 * mat.Q11 + 2. * mat.c12 * mat.Q12);
			mat.tv2 = 2. * (mat.c11 * mat.Q12 + mat.c12 * mat.Q11 + mat.c12 * mat.Q12);
			mat.tv3 = 4. * mat.c44 * mat.Q44;
			mat.tv4 = 2. * (mat.c11 * (mat.Q11 * mat.F11 + 2. * mat.Q12 * mat.F12) + \
				2. * mat.c12 * (mat.Q11 * mat.F12 + mat.Q12 * mat.F11 + mat.Q12 * mat.F12));
			mat.tv5 = 2. * (mat.c11 * (mat.Q11 * mat.F12 + mat.Q12 * mat.F11 + mat.Q12 * mat.F12) + \
				mat.c12 * (mat.Q11 * mat.F11 + mat.Q11 * mat.F12 + mat.Q12 * mat.F11 + 3. * mat.Q12 * mat.F12));
			mat.tv6 = 2. * mat.c44 * mat.Q44 * mat.F44;
			mat.tv7 = mat.c11 * (mat.F11 * mat.F11 + 2. * mat.F12 * mat.F12) + \
				2. * mat.c12 * (mat.F12 * mat.F12 + 2. * mat.F11 * mat.F12);
			mat.tv8 = mat.c44 * mat.F44 * mat.F44;

			file >> std::boolalpha >> mat.if_PEC; file.ignore(1000, '\n');
			file >> mat.r_permittivity; file.ignore(1000, '\n');
			file >> mat.conductivity; file.ignore(1000, '\n');
			file >> mat.comp_n1 >> mat.omega_plasma_n1 >> mat.tao_e_n1; file.ignore(1000, '\n');
			file >> mat.comp_n2 >> mat.omega_plasma_n2 >> mat.tao_e_n2; file.ignore(1000, '\n');
			file >> mat.comp_n3 >> mat.omega_plasma_n3 >> mat.tao_e_n3; file.ignore(1000, '\n');

			check_material(mat);
			get_cijkl(mat);

			material_parameters[i] = mat;
			//material_parameters.push_back(mat);
		}
	}
	file.close();
}

void geometry_parameters::readgeo() {
	std::ifstream file("system_setting.in");
	if (file.is_open()) {
		file >> nx_system >> ny_system >> nz_system; file.ignore(1000, '\n');
		file >> dx >> dy >> dz; file.ignore(1000, '\n');
		file >> std::boolalpha >> periodicX >> periodicY >> periodicZ; file.ignore(1000, '\n');
		file >> std::boolalpha >> if_read_struct; file.ignore(1000, '\n');
		file >> nx_work >> ny_work; file.ignore(1000, '\n');
		file >> num_layer_unit; file.ignore(1000, '\n');
		file >> num_periods >> id_firstlayer >> id_lastlayer_unit; file.ignore(1000, '\n');

		pt_nz_layer_unit = new unsigned int[num_layer_unit];

		for (long int i = 0; i < num_layer_unit; i++) {
			file >> pt_nz_layer_unit[i];
		}
	}
	file.close();

	set_geometry();
}

void global_parameters::read_prescribe_Eext() {
	std::ifstream file("Eext.in");
	if (file.is_open()) {
		for (long int i = 0; i < num_prescribe_Eext; i++) {
			file >> prescribe_Eext_time[i] >> prescribe_Ex[i] >> prescribe_Ey[i] >> prescribe_Ez[i];
		}
	}
	file.close();
}

void global_parameters::read_prescribe_Hext() {
	std::ifstream file("Hext.in");
	if (file.is_open()) {
		for (long int i = 0; i < num_prescribe_Hext; i++) {
			file >> prescribe_Hext_time[i] >> prescribe_Hx[i] >> prescribe_Hy[i] >> prescribe_Hz[i];
		}
	}
	file.close();
}

void global_parameters::read_Hext_nonunif() {
	unsigned long x = 0;
	unsigned long y = 0;
	unsigned long z = 0;
	double Hx = 0.;
	double Hy = 0.;
	double Hz = 0.;

	std::ifstream file("Hext_nonunif.in");
	if (file.is_open()) {
		for (unsigned long i = 0; i < nx * ny * nz; i++)
		{
			file >> x >> y >> z >> Hx >> Hy >> Hz;
			Hextx_nonunif(x - 1, y - 1, z - 1) = Hx;
			Hexty_nonunif(x - 1, y - 1, z - 1) = Hy;
			Hextz_nonunif(x - 1, y - 1, z - 1) = Hz;
		}
	}
	file.close();
}

void global_parameters::read_struct() {
	unsigned long x = 0;
	unsigned long y = 0;
	unsigned long z = 0;
	unsigned long id = 0;

	std::ifstream file("struct.in");
	if (file.is_open()) {
		for (unsigned long i = 0; i < nx * ny * nz; i++)
		{
			file >> x >> y >> z >> id;
			material_cell(x - 1, y - 1, z - 1) = id;
		}
	}
	file.close();
}


