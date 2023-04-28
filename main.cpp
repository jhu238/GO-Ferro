#define _USE_SINGLETON_MACRO
#pragma once
#include "_MACRO_SINGLETON.h"
#include"geometry.h"
#include"global.h"
#include"mathlib.h"
#include"ferroelectric_system.h"
#include"magnetic_system.h"
#include"elastic_system.h"
#include"EMdynamic_system.h"
#include"inoutput.h"
#include<iostream>

#include <chrono>
#include <ctime>

int main() {
	unsigned long long int nstep = 0;

    // read geometry and copy to device
	_GEO.readgeo();
	_GEO.copy_to_device();
	// read simulation setting
	_GLB.read_global();
	_GLB.copy_to_device();
	// set system size
	_IO.get_dimension(_GEO.nx_system, _GEO.ny_system, _GEO.nz_system, &(global_parameters::glb));
	// copy math library to device
	_MLB.copy_to_device();
	
	// define different systems
	ferroelectric_system fe; //KG equation
	magnetic_system mag; //LLG
	elastic_system elasto; //elastodynamic and elastostatic
	EMdynamic_system em; //electrodynamics

	// system initialize
	fe.initialize_host(/*&elasto, &em*/);
	mag.initialize_host(/*&elasto, &em*/);
	elasto.initialize_host(/*&mag, &fe*/);
	em.initialize_host(/*&mag, &fe*/);

	//------------Physical quantities input-------------//
	{
		// magnetization
		if (_GLB.if_input_m == true) {
			_IO.input_m(&mag);
		}
		// order paramters of AFM
		if (_GLB.if_input_AFMm == true) {
			_IO.input_AFMm(&mag);
		}
		// demagnetization field
		if (_GLB.if_input_Hstat == true) {
			_IO.input_Hstat(&mag);
		}
		// polarization and partialp/partialt
		if (_GLB.if_input_pandq == true) {
			_IO.input_pandq(&fe);
		}
		// electrostatic field
		if (_GLB.if_input_Estat == true) {
			_IO.input_Estat(&fe);
		}
		// 
		if (_GLB.if_input_em == true) {
			_IO.input_em_Yee(&em);
		}
		// displacement and velocity
		if (_GLB.if_input_uandv == true) {
			_IO.input_uandv(&elasto);
		}
		// elasto force
		if (_GLB.if_input_elastoforce == true) {
			_IO.input_elastoforce(&elasto);
		}
		// static strain(obtained by elastostatic equation
		if (_GLB.if_input_straint0 == true) {
			_IO.input_straint0(&elasto);
		}
		// eigenstrain (static: with 0; dynamic: without 0)
		if (_GLB.if_input_eigenstraint0_crt == true) {
			_IO.input_eigenstraint0(&elasto);
		}
		// 
		if (_GLB.if_input_Jp == true) {
			_IO.input_Jp(&em);
		}
	}

	fe.copy_to_device();
	mag.copy_to_device();
	elasto.copy_to_device();
	em.copy_to_device();

#pragma acc serial default(present)
	{
		// left: point; right: variables
		// pt -> pointer
		fe.pt_elas = &(elasto);
		fe.pt_EM = &(em);

		mag.pt_elas = &(elasto);
		mag.pt_EM = &(em);

		elasto.pt_mag = &(mag);
		elasto.pt_fe = &(fe);

		em.pt_mag = &(mag);
		em.pt_fe = &(fe);
	}

	mag.initialize_device();
	fe.initialize_device();
	em.initialize_device();
	elasto.initialize_device();

#pragma acc wait

	//------------------------Calculation for initial state output------------------//
	{
		if (_GLB.if_FE_all == true && _GLB.if_input_Estat == false) {
			if (_GLB.if_elec_1Dmodel == true) {
#pragma acc parallel default(present)
				{
					fe.get_E_static_1D();
				}
			}
			else {
				fe.get_E_static();
			}
			fe.copyE_from_device();
		}

		if ((_GLB.if_FM_all == true || _GLB.if_AFM_all == true) && _GLB.if_input_Hstat == false && _GLB.if_magnetostatics == true) {
			if (_GLB.if_mag_1Dmodel == true) {
#pragma acc parallel default(present)
				{
					mag.get_H_static_1D();
				}
			}
			else {
				mag.get_H_static();
			}
			mag.copyH_from_device();
		}

		if (_GLB.if_elastodynamics == true && _GLB.if_input_straint0 == false && _GLB.if_input_eigenstraint0_crt == false) {
			elasto.get_strain_static();
#pragma acc parallel default(present)
			{
				elasto.get_strain_static_glb();
			}

			elasto.copy_straint0_from_device();
			elasto.copy_eigenstraint0_from_device();
		}

		if (_GLB.if_FE_all == true && _GLB.if_flexo == true) {
#pragma acc parallel default(present)
			{
				fe.get_p_gradient();
			}
		}

#pragma acc wait
	}

	//------------------------Output for initial state------------------//
	{
		if (_GLB.if_output_ave == true) {
			if (_GLB.if_FM_all == true) {
				mag.get_averagem();
				_IO.output_averagem(nstep, &mag);
			}
			if (_GLB.if_AFM_all == true) {
				mag.get_averagem_AFM();
				_IO.output_averageAFMm(nstep, &mag);
			}
			if (_GLB.if_FE_all == true) {
				fe.get_averagep();
				_IO.output_averagep(nstep, &fe);
			}
		}

		if (_GLB.if_output_Estat == true) {
			_IO.output_Estat(nstep, &fe);
		}
		if (_GLB.if_output_Hstat == true) {
			_IO.output_Hstat(nstep, &mag);
		}

		if (_GLB.if_elastodynamics == true && _GLB.if_output_strain == true) {
			_IO.output_straint0(nstep, &elasto);
		}
		if (_GLB.if_output_eigenstraint0_crt == true) {
			_IO.output_eigenstraint0_crt(nstep, &elasto);
		}

		//_GLB.Hext[2] = 1.e7;
		//_GLB.update_to_device();
	}
#pragma acc wait

	//------------------------------//
	//		Iteration starts		//
	//------------------------------//

	auto start = std::chrono::system_clock::now();
	std::time_t start_time = std::chrono::system_clock::to_time_t(start);
	std::cout << "started computation at " << std::ctime(&start_time);

	if (_GLB.if_elastodynamics == false) {
		elasto.get_strain_static();
	}

	if (_GLB.if_EMdynamic == false) {
		if (_GLB.if_FE_all == true) {
			if (_GLB.if_elec_1Dmodel == true) {
#pragma acc parallel default(present)
				{
					fe.get_E_static_1D();
				}
			}
			else {
				fe.get_E_static();
			}
		}
		if ((_GLB.if_FM_all == true || _GLB.if_AFM_all == true) && _GLB.if_magnetostatics == true) {
			if (_GLB.if_mag_1Dmodel == true) {
#pragma acc parallel default(present)
				{
					mag.get_H_static_1D();
				}
			}
			else {
				mag.get_H_static();
			}
		}
	}
#pragma acc wait

	for (nstep = 1; nstep < _GLB.step + 1; ++nstep)
	{
#pragma acc serial default(present)
		{
			// simulation timestep (in serial step)
			_GLB.time_marching_device_serial();
		}
#pragma acc wait

		//--------------------Stage 1 for RK4---------------------//
		{
			if (_GLB.if_elastodynamics == true) {
				elasto.get_dv_RK1();
				elasto.get_du_RK1();
			}

			if (_GLB.if_FM_all == true) {
				mag.get_dm_RK1();
			}
			if (_GLB.if_AFM_all == true) {
				mag.get_dm_AFM_RK1();
			}
			if (_GLB.if_FE_all == true) {
				fe.get_dq_RK1();
				fe.get_dp_RK1();
			}

			if (_GLB.if_EMdynamic == true) {
				em.get_dE_RK1();
				em.get_dH_RK1();
			}
#pragma acc wait

			if (_GLB.if_elastodynamics == true) {
				elasto.update_v_RK1();
				elasto.update_u_RK1();
			}
			if (_GLB.if_FM_all == true) {
				mag.update_m_RK1();
			}
			if (_GLB.if_AFM_all == true) {
				mag.update_m_AFM_RK1();
			}
			if (_GLB.if_FE_all == true) {
				fe.update_q_RK1();
				fe.update_p_RK1();
			}

			if (_GLB.if_EMdynamic == true) {
				em.update_DE_RK1();
				em.update_DH_RK1();
			}
#pragma acc wait
		}
		//--------------------END Stage 1 for RK4---------------------//

		if (_GLB.if_EMdynamic == false) {
			if (_GLB.if_FE_all == true && _GLB.if_elec_1Dmodel == false) {
				fe.get_E_static();
			}
			if ((_GLB.if_FM_all == true || _GLB.if_AFM_all == true) && _GLB.if_mag_1Dmodel == false && _GLB.if_magnetostatics == true) {
				mag.get_H_static();
			}
		}
#pragma acc wait

		//--------------------Stage 2 for RK4---------------------//
		{
			if (_GLB.if_elastodynamics == true) {
				elasto.get_dv_RK2();
				elasto.get_du_RK2();
			}

			if (_GLB.if_FM_all == true) {
				mag.get_dm_RK2();
			}
			if (_GLB.if_AFM_all == true) {
				mag.get_dm_AFM_RK2();
			}
			if (_GLB.if_FE_all == true) {
				fe.get_dq_RK2();
				fe.get_dp_RK2();
			}

			if (_GLB.if_EMdynamic == true) {
				em.get_dE_RK2();
				em.get_dH_RK2();
			}
#pragma acc wait

			if (_GLB.if_elastodynamics == true) {
				elasto.update_v_RK2();
				elasto.update_u_RK2();
			}
			if (_GLB.if_FM_all == true) {
				mag.update_m_RK2();
			}
			if (_GLB.if_AFM_all == true) {
				mag.update_m_AFM_RK2();
			}
			if (_GLB.if_FE_all == true) {
				fe.update_q_RK2();
				fe.update_p_RK2();
			}

			if (_GLB.if_EMdynamic == true) {
				em.update_DE_RK2();
				em.update_DH_RK2();
			}
#pragma acc wait
		}
		//--------------------END Stage 2 for RK4---------------------//

		if (_GLB.if_EMdynamic == false) {
			if (_GLB.if_FE_all == true && _GLB.if_elec_1Dmodel == false) {
				fe.get_E_static();
			}
			if ((_GLB.if_FM_all == true || _GLB.if_AFM_all == true) && _GLB.if_mag_1Dmodel == false && _GLB.if_magnetostatics == true) {
				mag.get_H_static();
			}
		}
#pragma acc wait

		//--------------------Stage 3 for RK4---------------------//
		{
			if (_GLB.if_elastodynamics == true) {
				elasto.get_dv_RK3();
				elasto.get_du_RK3();
			}

			if (_GLB.if_FM_all == true) {
				mag.get_dm_RK3();
			}
			if (_GLB.if_AFM_all == true) {
				mag.get_dm_AFM_RK3();
			}
			if (_GLB.if_FE_all == true) {
				fe.get_dq_RK3();
				fe.get_dp_RK3();
			}

			if (_GLB.if_EMdynamic == true) {
				em.get_dE_RK3();
				em.get_dH_RK3();
			}
#pragma acc wait

			if (_GLB.if_elastodynamics == true) {
				elasto.update_v_RK3();
				elasto.update_u_RK3();
			}
			if (_GLB.if_FM_all == true) {
				mag.update_m_RK3();
			}
			if (_GLB.if_AFM_all == true) {
				mag.update_m_AFM_RK3();
			}
			if (_GLB.if_FE_all == true) {
				fe.update_q_RK3();
				fe.update_p_RK3();
			}

			if (_GLB.if_EMdynamic == true) {
				em.update_DE_RK3();
				em.update_DH_RK3();
			}
#pragma acc wait
		}
		//--------------------END Stage 3 for RK4---------------------//

		if (_GLB.if_EMdynamic == false) {
			if (_GLB.if_FE_all == true && _GLB.if_elec_1Dmodel == false) {
				fe.get_E_static();
			}
			if ((_GLB.if_FM_all == true || _GLB.if_AFM_all == true) && _GLB.if_mag_1Dmodel == false && _GLB.if_magnetostatics == true) {
				mag.get_H_static();
			}
		}
#pragma acc wait

		//--------------------Stage 4 for RK4---------------------//
		{
			if (_GLB.if_elastodynamics == true) {
				elasto.get_dv_RK4();
				elasto.get_du_RK4();
			}

			if (_GLB.if_FM_all == true) {
				mag.get_dm_RK4();
				mag.get_dm();
			}
			if (_GLB.if_AFM_all == true) {
				mag.get_dm_AFM_RK4();
				mag.get_dm_AFM();
			}
			if (_GLB.if_FE_all == true) {
				fe.get_dq_RK4();
				fe.get_dp_RK4();
			}

			if (_GLB.if_EMdynamic == true) {
				em.get_dE_RK4();
				em.get_dH_RK4();
			}
#pragma acc wait

			if (_GLB.if_elastodynamics == true) {
				elasto.update_v();
				elasto.update_u();
			}
			if (_GLB.if_FM_all == true) {
				mag.update_m();
			}
			if (_GLB.if_AFM_all == true) {
				mag.update_m_AFM();
			}
			if (_GLB.if_FE_all == true) {
				fe.update_q();
				fe.update_p();
			}

			if (_GLB.if_EMdynamic == true) {
				em.update_DE();
				em.update_DH();
			}
#pragma acc wait
		}
		//--------------------END Stage 4 for RK4---------------------//

		if (_GLB.if_elastodynamics == false) {
			elasto.get_strain_static();
		}

		if (_GLB.if_EMdynamic == false) {
			if (_GLB.if_FE_all == true && _GLB.if_elec_1Dmodel == false) {
				fe.get_E_static();
			}
			if ((_GLB.if_FM_all == true || _GLB.if_AFM_all == true) && _GLB.if_mag_1Dmodel == false && _GLB.if_magnetostatics == true) {
				mag.get_H_static();
			}
		}
#pragma acc wait

		//------------------------------//
		//				OUTPUT			//
		//------------------------------//

		//--------------------Copy from device-----------------//
		{
			if ((_GLB.if_output_m == true && nstep % _GLB.output_step_m == 0) || \
				(_GLB.if_output_only_magcell == true && nstep % _GLB.output_step_magcell == 0)) {				
					mag.copym_from_device();				
			}

			if (_GLB.if_output_AFMm == true) {
				if (nstep % _GLB.output_step_AFMm == 0) {
					mag.copym_AFM_from_device();
				}
			}

			if (_GLB.if_output_p == true) {
				if (nstep % _GLB.output_step_p == 0) {
					fe.copyp_from_device();
				}
			}

			if (_GLB.if_output_q == true) {
				if (nstep % _GLB.output_step_q == 0) {
					fe.copyq_from_device();
				}
			}

			if ((_GLB.if_output_em == true && nstep % _GLB.output_step_em == 0) || \
				(_GLB.if_output_em_onecell == true && nstep % _GLB.output_step_em_onecell == 0)) {
					em.copy_from_device();
			}

			if (_GLB.if_output_emYee == true) {
				if (nstep % _GLB.output_step_emYee == 0) {
					em.copyYee_from_device();
				}
			}

			if (_GLB.if_output_Estat == true) {
				if (nstep % _GLB.output_step_Estat == 0) {
					fe.copyE_from_device();
				}
			}

			if (_GLB.if_output_Hstat == true) {
				if (nstep % _GLB.output_step_Hstat == 0) {
					mag.copyH_from_device();
				}
			}

			if (_GLB.if_output_uandv == true) {
				if (nstep % _GLB.output_step_uandv == 0) {
					elasto.copy_uandv_from_device();
				}
			}

			if (_GLB.if_output_elastoforce == true) {
				if (nstep % _GLB.output_step_elastoforce == 0) {
					elasto.copy_elastoforce_from_device();
				}
			}

			if (_GLB.if_output_eigenstraint0_crt == true) {
				if (nstep % _GLB.output_step_eigenstraint0_crt == 0) {
					elasto.copy_eigenstraint0_from_device();
				}
			}

			if (_GLB.if_output_Jp == true) {
				if (nstep % _GLB.output_step_Jp == 0) {
					em.copyJp_from_device();
				}
			}

			if (_GLB.if_output_Jishe == true) {
				if (nstep % _GLB.output_step_Jishe == 0) {
					mag.copyJishe_from_device();
				}
			}

			if (_GLB.if_output_ave == true) {
				if (nstep % _GLB.output_step_ave == 0) {
					if (_GLB.if_FM_all == true) {
						mag.get_averagem();
					}
					if (_GLB.if_AFM_all == true) {
						mag.get_averagem_AFM();
					}
					if (_GLB.if_FE_all == true) {
						fe.get_averagep();
					}
				}
			}

			if (_GLB.if_output_strain == true) {
				if (nstep % _GLB.output_step_strain == 0) {
					if (_GLB.if_elastodynamics == true) {
						elasto.copy_Dstrain_from_device();
					}
					else {
						if (_GLB.if_elastostatic == true) {
#pragma acc parallel default(present)
							{
								elasto.get_strain_static_glb();
							}
							elasto.copy_straint0_from_device();
						}
					}
				}
			}
		}
#pragma acc wait
		//------------------------File Writing-------------------------------------//
		{
			if (_GLB.if_output_ave == true) {
				if (nstep % _GLB.output_step_ave == 0) {
					if (_GLB.if_FM_all == true) {
						_IO.output_averagem(nstep, &mag);
					}
					if (_GLB.if_AFM_all == true) {
						_IO.output_averageAFMm(nstep, &mag);
					}
					if (_GLB.if_FE_all == true) {
						_IO.output_averagep(nstep, &fe);
					}
				}
			}

			if (_GLB.if_output_m == true) {
				if (nstep % _GLB.output_step_m == 0) {
					_IO.output_m(nstep, &mag);
				}
			}

			if (_GLB.if_output_only_magcell == true) {
				if (nstep % _GLB.output_step_magcell == 0) {
					_IO.output_magcell(nstep, &mag);
				}
			}

			if (_GLB.if_output_AFMm == true) {
				if (nstep % _GLB.output_step_AFMm == 0) {
					_IO.output_AFMm(nstep, &mag);
				}
			}

			if (_GLB.if_output_p == true) {
				if (nstep % _GLB.output_step_p == 0) {
					_IO.output_p(nstep, &fe);
				}
			}

			if (_GLB.if_output_q == true) {
				if (nstep % _GLB.output_step_q == 0) {
					_IO.output_q(nstep, &fe);
				}
			}

			if (_GLB.if_output_em == true) {
				if (nstep % _GLB.output_step_em == 0) {
					_IO.output_Eem(nstep, &em);
					_IO.output_Hem(nstep, &em);
				}
			}

			if (_GLB.if_output_em_onecell == true) {
				if (nstep % _GLB.output_step_em_onecell == 0) {
					_IO.output_em_onecell(nstep, &em);
				}
			}

			if (_GLB.if_output_emYee == true) {
				if (nstep % _GLB.output_step_emYee == 0) {
					_IO.output_EemYee(nstep, &em);
					_IO.output_HemYee(nstep, &em);
				}
			}

			if (_GLB.if_output_Estat == true) {
				if (nstep % _GLB.output_step_Estat == 0) {
					_IO.output_Estat(nstep, &fe);
				}
			}

			if (_GLB.if_output_Hstat == true) {
				if (nstep % _GLB.output_step_Hstat == 0) {
					_IO.output_Hstat(nstep, &mag);
				}
			}

			if (_GLB.if_output_strain == true) {
				if (nstep % _GLB.output_step_strain == 0) {
					_IO.output_strain(nstep, &elasto);
				}
			}

			if (_GLB.if_output_uandv == true) {
				if (nstep % _GLB.output_step_uandv == 0) {
					_IO.output_uandv(nstep, &elasto);
				}
			}

			if (_GLB.if_output_elastoforce == true) {
				if (nstep % _GLB.output_step_elastoforce == 0) {
					_IO.output_elastoforce(nstep, &elasto);
				}
			}

			if (_GLB.if_output_eigenstraint0_crt == true) {
				if (nstep % _GLB.output_step_eigenstraint0_crt == 0) {
					_IO.output_eigenstraint0_crt(nstep, &elasto);
				}
			}

			if (_GLB.if_output_Jp == true) {
				if (nstep % _GLB.output_step_Jp == 0) {
					_IO.output_Jp(nstep, &em);
				}
			}

			if (_GLB.if_output_Jishe == true) {
				if (nstep % _GLB.output_step_Jishe == 0) {
					_IO.output_Jishe(nstep, &mag);
				}
			}
		}
	}

	//------------------------------//
	//		Iteration ends			//
	//------------------------------//

	auto end = std::chrono::system_clock::now();
	std::chrono::duration<double> elapsed_seconds = end - start;

	std::time_t end_time = std::chrono::system_clock::to_time_t(end);

	std::cout << "finished computation at " << std::ctime(&end_time)
		<< "elapsed time: " << elapsed_seconds.count() << "s\n";

	return 0;
}

//Old Output code 

//if (_GLB.if_output_ave == true) {
//	if (nstep % _GLB.output_step_ave == 0) {
//		if (_GLB.if_FM_all == true) {
//			mag.get_averagem();
//			_IO.output_averagem(nstep, &mag);
//		}
//		if (_GLB.if_FE_all == true) {
//			fe.get_averagep();
//			_IO.output_averagep(nstep, &fe);
//		}
//	}
//}
//
//if (_GLB.if_output_m == true) {
//	if (nstep % _GLB.output_step_m == 0) {
//		mag.copym_from_device();
//#pragma acc wait
//		_IO.output_m(nstep, &mag);
//	}
//}
//
//if (_GLB.if_output_p == true) {
//	if (nstep % _GLB.output_step_p == 0) {
//		fe.copyp_from_device();
//#pragma acc wait
//		_IO.output_p(nstep, &fe);
//	}
//}
//
//if (_GLB.if_output_q == true) {
//	if (nstep % _GLB.output_step_q == 0) {
//		fe.copyq_from_device();
//#pragma acc wait
//		_IO.output_q(nstep, &fe);
//	}
//}
//
//if (_GLB.if_output_em == true) {
//	if (nstep % _GLB.output_step_em == 0) {
//		em.copy_from_device();
//#pragma acc wait
//		_IO.output_Eem(nstep, &em);
//		_IO.output_Hem(nstep, &em);
//	}
//}
//
//if (_GLB.if_output_emYee == true) {
//	if (nstep % _GLB.output_step_emYee == 0) {
//		em.copyYee_from_device();
//#pragma acc wait
//		_IO.output_EemYee(nstep, &em);
//		_IO.output_HemYee(nstep, &em);
//	}
//}
//
//if (_GLB.if_output_Estat == true) {
//	if (nstep % _GLB.output_step_Estat == 0) {
//		fe.copyE_from_device();
//#pragma acc wait
//		_IO.output_Estat(nstep, &fe);
//	}
//}
//
//if (_GLB.if_output_Hstat == true) {
//	if (nstep % _GLB.output_step_Hstat == 0) {
//		mag.copyH_from_device();
//#pragma acc wait
//		_IO.output_Hstat(nstep, &mag);
//	}
//}
//
//if (_GLB.if_output_strain == true) {
//	if (nstep % _GLB.output_step_strain == 0) {
//		if (_GLB.if_elastodynamics == true) {
//			elasto.copy_Dstrain_from_device();
//		}
//		else {
//			if (_GLB.if_elastostatic == true) {
//				elasto.get_strain_static_glb();
//				elasto.copy_straint0_from_device();
//			}
//		}
//#pragma acc wait
//		_IO.output_strain(nstep, &elasto);
//	}
//}
//
//if (_GLB.if_output_uandv == true) {
//	if (nstep % _GLB.output_step_uandv == 0) {
//		elasto.copy_uandv_from_device();
//#pragma acc wait
//		_IO.output_uandv(nstep, &elasto);
//	}
//}
//
//if (_GLB.if_output_eigenstraint0_crt == true) {
//	if (nstep % _GLB.output_step_eigenstraint0_crt == 0) {
//		elasto.copy_eigenstraint0_from_device();
//#pragma acc wait
//		_IO.output_eigenstraint0_crt(nstep, &elasto);
//	}
//}
//
//if (_GLB.if_output_Jp == true) {
//	if (nstep % _GLB.output_step_Jp == 0) {
//		em.copyJp_from_device();
//#pragma acc wait
//		_IO.output_Jp(nstep, &em);
//	}
//}
//
//if (_GLB.if_output_Jishe == true) {
//	if (nstep % _GLB.output_step_Jishe == 0) {
//		mag.copyJishe_from_device();
//#pragma acc wait
//		_IO.output_Jishe(nstep, &mag);
//	}
//}
