 #pragma once

class geometry_parameters {
private:
	unsigned int* pt_nz_layer_unit = 0;
	long int num_periods; /*number of repeating unit in the superlattice*/
	long int nx_work, ny_work, nz_work;

private:
	void set_geometry();

public:
	bool periodicX, periodicY, periodicZ;
	long int num_layer_unit;
	long int id_lastlayer_unit;
	long int id_firstlayer, id_lastlayer; /*id of the first and the last layer in the superlattice*/

	long int num_layer;
	unsigned int* pt_nz_layer = 0;

	long int nx_system, ny_system, nz_system;

	long int idxi_work = 0, idyi_work = 0, idzi_work = 0;
	long int idxf_work = 0, idyf_work = 0, idzf_work = 0;

	double dx, dy, dz;

	bool if_read_struct;
	static geometry_parameters geo;

public:
	void readgeo();
	void copy_to_device();
};
