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

    int me,nprocs;

    double delxinv,delyinv,delzinv,delvolinv; //setup
    double volume;

    double shiftone; //make_rho
    int nlower, nupper;


    FFT_SCALAR **rho1d,**rho_coeff,**drho1d,**drho_coeff; //compute_rho1d

    FILE *fp;

    int nxlo_in,nylo_in,nzlo_in,nxhi_in,nyhi_in,nzhi_in;
    int nxlo_out,nylo_out,nzlo_out,nxhi_out,nyhi_out,nzhi_out;
    int nxlo_ghost,nxhi_ghost,nylo_ghost,nyhi_ghost,nzlo_ghost,nzhi_ghost;
    int nxlo_fft,nylo_fft,nzlo_fft,nxhi_fft,nyhi_fft,nzhi_fft;
    int ngrid,nfft,nfft_both;

    class FFT3d *fft1,*fft2;
    class Remap *remap;
    class GridComm *cg;
    class GridComm *cg_peratom;

    FFT_SCALAR ***density_brick;
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
};


}


#endif
#endif
