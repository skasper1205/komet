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
	void final_integrate();
	void kspace_check();
    void pppm_check();
    void init_condiff();
      void setup();
      void setup_grid();
      void compute();
      int timing_1d(int, double &);
      int timing_3d(int, double &);
      double memory_usage();

      void compute_group_group(int, int, int);

      int pppmflag;
      int group_group_enable;
      double scale,qqrd2e;
      double accuracy;
      double accuracy_absolute;
      double accuracy_relative;
      int order, order_6,order_allocated;;
      int stagger_flag;
      int differentiation_flag;
      double slab_volfactor;
      int nx_pppm, ny_pppm, nz_pppm;
      double g_ewald,g_ewald_6;
      int gridflag;
      double e2group;
      double f2group[3];
      double two_charge_force;
      int collective_flag;


     protected:
      int me,nprocs;
      int nfactors;
      int *factors;
      double cutoff;
      double volume;
      double delxinv,delyinv,delzinv,delvolinv;
      double h_x,h_y,h_z;
      double shift,shiftone;
      int peratom_allocate_flag;

      int nxlo_in,nylo_in,nzlo_in,nxhi_in,nyhi_in,nzhi_in;
      int nxlo_out,nylo_out,nzlo_out,nxhi_out,nyhi_out,nzhi_out;
      int nxlo_ghost,nxhi_ghost,nylo_ghost,nyhi_ghost,nzlo_ghost,nzhi_ghost;
      int nxlo_fft,nylo_fft,nzlo_fft,nxhi_fft,nyhi_fft,nzhi_fft;
      int nlower,nupper;
      int ngrid,nfft,nfft_both;

      FFT_SCALAR ***density_brick;
      FFT_SCALAR ***vdx_brick,***vdy_brick,***vdz_brick;
      FFT_SCALAR ***u_brick;
      FFT_SCALAR ***v0_brick,***v1_brick,***v2_brick;
      FFT_SCALAR ***v3_brick,***v4_brick,***v5_brick;
      double *greensfn;
      double **vg;
      double *fkx,*fky,*fkz;
      FFT_SCALAR *density_fft;
      FFT_SCALAR *work1,*work2;

      double *gf_b;
      FFT_SCALAR **rho1d,**rho_coeff,**drho1d,**drho_coeff;
      double *sf_precoeff1, *sf_precoeff2, *sf_precoeff3;
      double *sf_precoeff4, *sf_precoeff5, *sf_precoeff6;
      double sf_coeff[6];          // coefficients for calculating ad self-forces
      double **acons;

      // group-group interactions

      int group_allocate_flag;
      FFT_SCALAR ***density_A_brick,***density_B_brick;
      FFT_SCALAR *density_A_fft,*density_B_fft;

      class FFT3d *fft1,*fft2;
      class Remap *remap;
      class GridComm *cg;
      class GridComm *cg_peratom;

      int **part2grid;             // storage for particle -> grid mapping
      int nmax;

      double *boxlo;
                                   // TIP4P settings
      int typeH,typeO;             // atom types of TIP4P water H and O atoms
      double qdist;                // distance from O site to negative charge
      double alpha;                // geometric factor

      void set_grid_global();
      void set_grid_local();
      void adjust_gewald();
      double newton_raphson_f();
      double derivf();
      double final_accuracy();

      void allocate();
      void allocate_peratom();
      void deallocate();
      void deallocate_peratom();
      int factorable(int);
      double compute_df_kspace();
      double estimate_ik_error(double, double, bigint);
      double compute_qopt();
      void compute_gf_denom();
      void compute_gf_ik();
      void compute_gf_ad();
      void compute_sf_precoeff();

      void particle_map();
      void make_rho();
      void brick2fft();

      void poisson();
      void poisson_ik();
      void poisson_ad();

      void fieldforce();
      void fieldforce_ik();
      void fieldforce_ad();

      void poisson_peratom();
      void fieldforce_peratom();
      void procs2grid2d(int,int,int,int *, int*);
      void compute_rho1d(const FFT_SCALAR &, const FFT_SCALAR &,
                         const FFT_SCALAR &);
      void compute_drho1d(const FFT_SCALAR &, const FFT_SCALAR &,
                         const FFT_SCALAR &);
      void compute_rho_coeff();
      void slabcorr();

      // grid communication

      void pack_forward(int, FFT_SCALAR *, int, int *);
      void unpack_forward(int, FFT_SCALAR *, int, int *);
      void pack_reverse(int, FFT_SCALAR *, int, int *);
      void unpack_reverse(int, FFT_SCALAR *, int, int *);

      // triclinic

      int triclinic;               // domain settings, orthog or triclinic
      void setup_triclinic();
      void compute_gf_ik_triclinic();
      void poisson_ik_triclinic();
      void poisson_groups_triclinic();

      // group-group interactions

      void allocate_groups();
      void deallocate_groups();
       void make_rho_groups(int, int, int);
      void poisson_groups(int);
      void slabcorr_groups(int,int,int);

    /* ----------------------------------------------------------------------
       denominator for Hockney-Eastwood Green's function
         of x,y,z = sin(kx*deltax/2), etc

                inf                 n-1
       S(n,k) = Sum  W(k+pi*j)**2 = Sum b(l)*(z*z)**l
               j=-inf               l=0

              = -(z*z)**n /(2n-1)! * (d/dx)**(2n-1) cot(x)  at z = sin(x)
       gf_b = denominator expansion coeffs
    ------------------------------------------------------------------------- */

      inline double gf_denom(const double &x, const double &y,
                             const double &z) const {
        double sx,sy,sz;
        sz = sy = sx = 0.0;
        for (int l = order-1; l >= 0; l--) {
          sx = gf_b[l] + sx*x;
          sy = gf_b[l] + sy*y;
          sz = gf_b[l] + sz*z;
        }
        double s = sx*sy*sz;
        return s*s;
      };
    };

    }

    #endif
    #endif
