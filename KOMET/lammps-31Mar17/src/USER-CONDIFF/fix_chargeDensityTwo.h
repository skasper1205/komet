/*
 * fix_condiff.h
 *
 *  Created on: Jan 16, 2017
 *      Author: skasper
 */

#ifdef FIX_CLASS

FixStyle(chargeDensityTwo, FixChargeDensityTwo)

#else

#ifndef LMP_FIX_CHARGEDENSITYTWO_H
#define LMP_FIX_CHARGEDENSITYTWO_H

#include "lmptype.h"
#include <mpi.h>

#ifdef FFT_SINGLE
typedef float FFT_SCALAR;
#define MPI_FFT_SCALAR MPI_FLOAT
#else
typedef double FFT_SCALAR;
#define MPI_FFT_SCALAR MPI_DOUBLE
#endif

#include "fix.h"

namespace LAMMPS_NS {

class FixChargeDensityTwo : public Fix {

  public:
	FixChargeDensityTwo(class LAMMPS *, int, char **);
	~FixChargeDensityTwo();
	int setmask();
	void end_of_step();
	void calc_cm();
	void calc_charge_density();

  protected:
	int jgroup, kgroup, groupbit_coions, groupbit_counterions;
        int me,nprocs;
	double masstotal;
       	FILE * file;
	double *cm;
	int nevery;
	
};


}


#endif
#endif
