#include "geometry.h"

void geometry_parameters::set_geometry() {
	long int num_unitlayers;

	num_unitlayers = id_lastlayer_unit - id_firstlayer + 1;
	id_lastlayer = id_firstlayer - 1 + num_periods * num_unitlayers;
	num_layer = num_layer_unit + (num_periods - 1) * num_unitlayers;

	pt_nz_layer = new unsigned int[num_layer];
	long int j = id_firstlayer;
	for (long int i = 0; i < num_layer; i++) {
		if (i < id_firstlayer - 1) {
			pt_nz_layer[i] = pt_nz_layer_unit[i];
		}
		else if (i > id_lastlayer - 1) {
			pt_nz_layer[i] = pt_nz_layer_unit[i - num_layer + num_layer_unit];
		}
		else {
			pt_nz_layer[i] = pt_nz_layer_unit[j - 1];
			j = j + 1;
			if (j > id_lastlayer_unit) {
				j = id_firstlayer;
			}
		}
	}

	nz_work = 0;
	for (long int i = 0; i < num_layer; i++) {
		nz_work = nz_work + pt_nz_layer[i];
	}

	idxi_work = static_cast<long int>((nx_system - nx_work) / 2) + 1;
	idxf_work = idxi_work + nx_work - 1;
	idyi_work = static_cast<long int>((ny_system - ny_work) / 2) + 1;
	idyf_work = idyi_work + ny_work - 1;
	//idzi_work = static_cast<long int>((nz_system - nz_work) / 2) + 1;
	//idzf_work = idzi_work + nz_work - 1;
	idzi_work = 1;
	idzf_work = idzi_work + nz_work - 1;

	idxi_work = idxi_work - 1; idyi_work = idyi_work - 1; idzi_work = idzi_work - 1;
	idxf_work = idxf_work - 1; idyf_work = idyf_work - 1; idzf_work = idzf_work - 1;
}

void geometry_parameters::copy_to_device() {
#pragma acc enter data copyin(this)
//#pragma acc enter data copyin(this->pt_nz_layer)
//#pragma acc enter data copyin(this->pt_nz_layer[0:num_layer])
}

geometry_parameters geometry_parameters::geo;
