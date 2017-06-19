/*
 * fix_dipoleMoment.h
 *
 *  Created on: Apr 11, 2017
 *      Author: skasper
 */

#ifdef FIX_CLASS

FixStyle(dipoleMomentTwo, FixDipoleMomentTwo)

#else

#ifndef LMP_FIX_DIPOLEMOMENTTWO_H
#define LMP_FIX_DIPOLEMOMENTTWO_H

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

class FixDipoleMomentTwo : public Fix {

  public:
	FixDipoleMomentTwo(class LAMMPS *, int, char **);
	~FixDipoleMomentTwo();
	int setmask();
	void end_of_step();
	void calc_cm();
	void calc_dipole_moment();

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
