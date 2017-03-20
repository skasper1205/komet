/*
 * fix_condiff.h
 *
 *  Created on: Oct 24, 2016
 *      Author: skasper
 */

#ifdef FIX_CLASS

FixStyle(condiff, FixCondiff)

#else

#ifndef LMP_FIX_CONDIFF_H
#define LMP_FIX_CONDIFF_H

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

class FixCondiff : public Fix {

  public:
	FixCondiff(class LAMMPS *, int, char **);
	~FixCondiff();
	int setmask();
	void post_force(int);
	void kspace_check();
    void pppm_check();
    void apply_boundary_conditions();


  protected:
    void setup();
    virtual void particle_map();
    virtual void make_rho();
    void reverse_make_rho();
    void rho_decompose();
    void deallocate();
    void allocate();
    void compute_rho1d(const FFT_SCALAR &, const FFT_SCALAR &,
            const FFT_SCALAR &);
    void compute_rho_coeff();
    void compute_drho1d(const FFT_SCALAR &, const FFT_SCALAR &, const
            FFT_SCALAR &);
    void set_grid_local();
    void setup_grid();

    FFT_SCALAR ***density_brick_velocity_x_loc;
    FFT_SCALAR ***density_brick_velocity_y_loc;
    FFT_SCALAR ***density_brick_velocity_z_loc;
    
    FFT_SCALAR ***density_brick_velocity_x_glo;
    FFT_SCALAR ***density_brick_velocity_y_glo;
    FFT_SCALAR ***density_brick_velocity_z_glo;

    FFT_SCALAR ***density_brick_force_x_loc;
    FFT_SCALAR ***density_brick_force_y_loc;
    FFT_SCALAR ***density_brick_force_z_loc;
    
    FFT_SCALAR ***density_brick_force_x_glo;
    FFT_SCALAR ***density_brick_force_y_glo;
    FFT_SCALAR ***density_brick_force_z_glo;
    
    FFT_SCALAR ***density_brick_counter_loc;
    FFT_SCALAR ***density_brick_counter_glo;

    double *rand;

    FFT_SCALAR **rho1d,**rho_coeff,**drho1d,**drho_coeff; //compute_rho1d

    class GridComm *cg;

    int order, minorder, order_allocated;

    int jgroup, groupbit_condiff;

    double T;
    double D;
    int seed;

    int nx_pppm, ny_pppm, nz_pppm; //particle_map

    double *boxlo;
    double shift;

    double dt;

    double wienerConst;


    int me,nprocs;

    double delxinv,delyinv,delzinv,delvolinv; //setup
    double volume;

    double shiftone; //make_rho
    int nlower, nupper;

    int nxlo_in,nylo_in,nzlo_in,nxhi_in,nyhi_in,nzhi_in;
    int nxlo_out,nylo_out,nzlo_out,nxhi_out,nyhi_out,nzhi_out;
    int nxlo_ghost,nxhi_ghost,nylo_ghost,nyhi_ghost,nzlo_ghost,nzhi_ghost;
    int nxlo_fft,nylo_fft,nzlo_fft,nxhi_fft,nyhi_fft,nzhi_fft;
    int ngrid,nfft,nfft_both;

    class RanMars *random;

    int nmax;

    double time1;
    double time2;
};


}


#endif
#endif
