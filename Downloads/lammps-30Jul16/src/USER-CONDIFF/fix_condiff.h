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
	void end_of_step();
	void kspace_check();
    void pppm_check();
    void apply_boundary_conditions();


  protected:
    void setup();
    virtual void particle_map();
    virtual void make_rho();
    void reverse_make_rho();
    void deallocate();
    void allocate();
    void compute_rho1d(const FFT_SCALAR &, const FFT_SCALAR &,
            const FFT_SCALAR &);
    void compute_rho_coeff();
    void compute_drho1d(const FFT_SCALAR &, const FFT_SCALAR &, const
            FFT_SCALAR &);
    void set_grid_local();
    void setup_grid();

    int nx_pppm, ny_pppm, nz_pppm; //particle_map
    int order, minorder, natoms;
    double *boxlo;
    double shift;
    int **part2grid;
    double dt;
    double D;
    double T;
    double wienerConst;
    double *rand;

    int me,nprocs;

    double delxinv,delyinv,delzinv,delvolinv; //setup
    double volume;

    double shiftone; //make_rho
    int nlower, nupper;


    FFT_SCALAR **rho1d,**rho_coeff,**drho1d,**drho_coeff; //compute_rho1d

    int nxlo_in,nylo_in,nzlo_in,nxhi_in,nyhi_in,nzhi_in;
    int nxlo_out,nylo_out,nzlo_out,nxhi_out,nyhi_out,nzhi_out;
    int nxlo_ghost,nxhi_ghost,nylo_ghost,nyhi_ghost,nzlo_ghost,nzhi_ghost;
    int nxlo_fft,nylo_fft,nzlo_fft,nxhi_fft,nyhi_fft,nzhi_fft;
    int ngrid,nfft,nfft_both;

    class GridComm *cg;
    class GridComm *cg_peratom;

    char *id_temp;
    class Compute *temperature;
    int tflag;
    int which;

    class RanMars *random;

    FFT_SCALAR ****density_brick;
    FFT_SCALAR ****density_brick_counter;
    FFT_SCALAR ****density_brick_force;
    double **help_v;
    double **help_f;
    double **help_x;
    /*FFT_SCALAR ***vdx_brick,***vdy_brick,***vdz_brick;
    FFT_SCALAR ***u_brick;
    FFT_SCALAR ***v0_brick,***v1_brick,***v2_brick;
    FFT_SCALAR ***v3_brick,***v4_brick,***v5_brick;
    double *greensfn;
    double **vg;
    double *fkx,*fky,*fkz;
    FFT_SCALAR *density_fft;
    FFT_SCALAR *work1,*work2;*/

    int order_allocated;
    double *gf_b;

    int jgroup,groupbit_condiff;


};


}


#endif
#endif
