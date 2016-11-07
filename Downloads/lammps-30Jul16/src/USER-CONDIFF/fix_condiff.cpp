/*
 * fix_condiff.cpp
 *
 *  Created on: Oct 24, 2016
 *      Author: skasper
 */


#include "fix_condiff.h"
#include <mpi.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "pppm.h"
#include "atom.h"
#include "comm.h"
#include "gridcomm.h"
#include "neighbor.h"
#include "force.h"
#include "pair.h"
#include "bond.h"
#include "angle.h"
#include "domain.h"
#include "fft3d_wrap.h"
#include "remap_wrap.h"
#include "memory.h"
#include "error.h"

#include "math_const.h"
#include "math_special.h"


using namespace LAMMPS_NS;
using namespace FixConst;
using namespace MathConst;
using namespace MathSpecial;

#define OFFSET 16384

#ifdef FFT_SINGLE
#define ZEROF 0.0f
#define ONEF  1.0f
#else
#define ZEROF 0.0
#define ONEF  1.0
#endif

//The class constructor.
FixCondiff::FixCondiff(LAMMPS *lmp, int narg, char **arg) :
		Fix(lmp, narg, arg)
{
	if (narg < 7) error->all(FLERR, "Illegal fix condiff command");

	kspace_check();
	pppm_check();


	 MPI_Comm_rank(world,&me);
	 MPI_Comm_size(world,&nprocs);
	/*int n = strlen(arg[3]) + 1; //saves 4th argument of script in a string
	string = new char[n];
	strcpy(string,arg[3]);*/



}

FixCondiff::~FixCondiff(){

}

int FixCondiff::setmask(){ //where algorithm steps in
	int mask = 0;
	mask |= FINAL_INTEGRATE;
	return mask;
}

void FixCondiff::final_integrate(){
	//setup();
	//init();
	//particle_map();

}

void FixCondiff::kspace_check(){   //check if pppm computation is used

	if (!force->kspace)
		error->all(FLERR, "Not using kspace computations");

}

void FixCondiff::pppm_check(){   //check if pppm computation is used

	if (!force->pair->pppmflag)
		error->all(FLERR, "Not using pppm computations");

}


void FixCondiff::setup()
{
	nx_pppm = force->kspace->nx_pppm;
	ny_pppm = force->kspace->ny_pppm;
	nz_pppm = force->kspace->nz_pppm;

	order =  force->kspace->order;
	minorder=2;
	order_allocated = order;

	nlower = -(order-1)/2;
	nupper = order/2;

	if (order % 2) shift = OFFSET + 0.5;
	else shift = OFFSET;
	if (order % 2) shiftone = 0.0;
	else shiftone = 0.5;


	double *prd;
	prd = domain->prd;

	double xprd = prd[0];
	double yprd = prd[1];
	double zprd = prd[2];
	double zprd_slab = zprd;
	volume = xprd * yprd * zprd_slab;

	delxinv = nx_pppm/xprd;
	delyinv = ny_pppm/yprd;
	delzinv = nz_pppm/zprd_slab;

	delvolinv = delxinv*delyinv*delzinv;
	//compute_rho_coeff();
}

void FixCondiff::particle_map()
{
	memory->destroy(part2grid);
	natoms = atom->nlocal;
	memory->create(part2grid,natoms,3,"condiff:part2grid");


	int nx,ny,nz;
    boxlo = domain->boxlo;

	double **x = atom->x;

	for (int i = 0; i < natoms; i++) {

    // (nx,ny,nz) = global coords of grid pt to "lower left" of charge
    // current particle coord can be outside global and local box
    // add/subtract OFFSET to avoid int(-0.75) = 0 when want it to be -1

		nx = static_cast<int> ((x[i][0]-boxlo[0])*delxinv+shift) - OFFSET;
		ny = static_cast<int> ((x[i][1]-boxlo[1])*delyinv+shift) - OFFSET;
		nz = static_cast<int> ((x[i][2]-boxlo[2])*delzinv+shift) - OFFSET;

	    part2grid[i][0] = nx;
		part2grid[i][1] = ny;
		part2grid[i][2] = nz;


	}

}

void FixCondiff::make_rho()
{
	int l,m,n,nx,ny,nz,mx,my,mz;
	FFT_SCALAR dx,dy,dz,x0,y0,z0;

  // clear 3d density array

	//memset(&(density_brick[nzlo_out][nylo_out][nxlo_out]),0,
    //     ngrid*sizeof(FFT_SCALAR));

  // loop over my velocities, add their contribution to nearby grid points
  // (nx,ny,nz) = global coords of grid pt to "lower left" of charge
  // (dx,dy,dz) = distance to "lower left" grid pt
  // (mx,my,mz) = global coords of moving stencil pt

	double *q = atom->q;
	double **x = atom->x;
	int nlocal = atom->nlocal;

	for (int i = 0; i < nlocal; i++) {

		nx = part2grid[i][0];
		ny = part2grid[i][1];
		nz = part2grid[i][2];
		dx = nx+shiftone - (x[i][0]-boxlo[0])*delxinv;
		dy = ny+shiftone - (x[i][1]-boxlo[1])*delyinv;
		dz = nz+shiftone - (x[i][2]-boxlo[2])*delzinv;



		//compute_rho1d(dx,dy,dz);

		/*z0 = delvolinv * q[i];
		for (n = nlower; n <= nupper; n++) {
			mz = n+nz;
			y0 = z0*rho1d[2][n];
			for (m = nlower; m <= nupper; m++) {
				my = m+ny;
				x0 = y0*rho1d[1][m];
				for (l = nlower; l <= nupper; l++) {
					mx = l+nx;
					density_brick[mz][my][mx] += x0*rho1d[0][l];


				}

			}
		}*/
	}
}

void FixCondiff::compute_rho1d(const FFT_SCALAR &dx, const FFT_SCALAR &dy,
                         const FFT_SCALAR &dz)
{
  int k,l;
  FFT_SCALAR r1,r2,r3;

  for (k = (1-order)/2; k <= order/2; k++) {
    r1 = r2 = r3 = ZEROF;

    for (l = order-1; l >= 0; l--) {
      r1 = rho_coeff[l][k] + r1*dx;
      r2 = rho_coeff[l][k] + r2*dy;
      r3 = rho_coeff[l][k] + r3*dz;
    }
    rho1d[0][k] = r1;
    rho1d[1][k] = r2;
    rho1d[2][k] = r3;
  }
}

void FixCondiff::compute_rho_coeff()
{
  int j,k,l,m;
  FFT_SCALAR s;

  FFT_SCALAR **a;

  memory->create2d_offset(a,order,-order,order,"condiff:a");


  for (k = -order; k <= order; k++)
    for (l = 0; l < order; l++)
      a[l][k] = 0.0;


  a[0][0] = 1.0;
  for (j = 1; j < order; j++) {
    for (k = -j; k <= j; k += 2) {
      s = 0.0;
      for (l = 0; l < j; l++) {
        a[l+1][k] = (a[l][k+1]-a[l][k-1]) / (l+1);
#ifdef FFT_SINGLE
        s += powf(0.5,(float) l+1) *
          (a[l][k-1] + powf(-1.0,(float) l) * a[l][k+1]) / (l+1);
#else
        s += pow(0.5,(double) l+1) *
          (a[l][k-1] + pow(-1.0,(double) l) * a[l][k+1]) / (l+1);
#endif
      }
      a[0][k] = s;
    }
  }

  m = (1-order)/2;
  for (k = -(order-1); k < order; k += 2) {
    for (l = 0; l < order; l++)
      rho_coeff[l][m] = a[l][k];
    for (l = 1; l < order; l++)
      drho_coeff[l-1][m] = l*a[l][k];
    m++;
  }

  memory->destroy2d_offset(a,-order);
}

void FixCondiff::compute_drho1d(const FFT_SCALAR &dx, const FFT_SCALAR &dy,
                          const FFT_SCALAR &dz)
{
  int k,l;
  FFT_SCALAR r1,r2,r3;

  for (k = (1-order)/2; k <= order/2; k++) {
    r1 = r2 = r3 = ZEROF;

    for (l = order-2; l >= 0; l--) {
      r1 = drho_coeff[l][k] + r1*dx;
      r2 = drho_coeff[l][k] + r2*dy;
      r3 = drho_coeff[l][k] + r3*dz;
    }
    drho1d[0][k] = r1;
    drho1d[1][k] = r2;
    drho1d[2][k] = r3;
  }
}

void FixCondiff::set_grid_local()
{
  // global indices of PPPM grid range from 0 to N-1
  // nlo_in,nhi_in = lower/upper limits of the 3d sub-brick of
  //   global PPPM grid that I own without ghost cells
  // for slab PPPM, assign z grid as if it were not extended

  nxlo_in = static_cast<int> (comm->xsplit[comm->myloc[0]] * nx_pppm);
  nxhi_in = static_cast<int> (comm->xsplit[comm->myloc[0]+1] * nx_pppm) - 1;

  nylo_in = static_cast<int> (comm->ysplit[comm->myloc[1]] * ny_pppm);
  nyhi_in = static_cast<int> (comm->ysplit[comm->myloc[1]+1] * ny_pppm) - 1;

  nzlo_in = static_cast<int>
      (comm->zsplit[comm->myloc[2]] * nz_pppm);
  nzhi_in = static_cast<int>
      (comm->zsplit[comm->myloc[2]+1] * nz_pppm) - 1;

  // nlower,nupper = stencil size for mapping particles to PPPM grid

  nlower = -(order-1)/2;
  nupper = order/2;

  // shift values for particle <-> grid mapping
  // add/subtract OFFSET to avoid int(-0.75) = 0 when want it to be -1

  if (order % 2) shift = OFFSET + 0.5;
  else shift = OFFSET;
  if (order % 2) shiftone = 0.0;
  else shiftone = 0.5;

  // nlo_out,nhi_out = lower/upper limits of the 3d sub-brick of
  //   global PPPM grid that my particles can contribute charge to
  // effectively nlo_in,nhi_in + ghost cells
  // nlo,nhi = global coords of grid pt to "lower left" of smallest/largest
  //           position a particle in my box can be at
  // dist[3] = particle position bound = subbox + skin/2.0 + qdist
  //   qdist = offset due to TIP4P fictitious charge
  //   convert to triclinic if necessary
  // nlo_out,nhi_out = nlo,nhi + stencil size for particle mapping
  // for slab PPPM, assign z grid as if it were not extended

  double *prd,*sublo,*subhi;

  prd = domain->prd;
  boxlo = domain->boxlo;
  sublo = domain->sublo;
  subhi = domain->subhi;


  double xprd = prd[0];
  double yprd = prd[1];
  double zprd = prd[2];
  double zprd_slab = zprd;
  double qdist = 0.0;
  double dist[3];
  double cuthalf = 0.5*neighbor->skin + qdist;
  dist[0] = dist[1] = dist[2] = cuthalf;

  int nlo,nhi;

  nlo = static_cast<int> ((sublo[0]-dist[0]-boxlo[0]) *
                            nx_pppm/xprd + shift) - OFFSET;
  nhi = static_cast<int> ((subhi[0]+dist[0]-boxlo[0]) *
                            nx_pppm/xprd + shift) - OFFSET;
  nxlo_out = nlo + nlower;
  nxhi_out = nhi + nupper;

  nlo = static_cast<int> ((sublo[1]-dist[1]-boxlo[1]) *
                            ny_pppm/yprd + shift) - OFFSET;
  nhi = static_cast<int> ((subhi[1]+dist[1]-boxlo[1]) *
                            ny_pppm/yprd + shift) - OFFSET;
  nylo_out = nlo + nlower;
  nyhi_out = nhi + nupper;

  nlo = static_cast<int> ((sublo[2]-dist[2]-boxlo[2]) *
                            nz_pppm/zprd_slab + shift) - OFFSET;
  nhi = static_cast<int> ((subhi[2]+dist[2]-boxlo[2]) *
                            nz_pppm/zprd_slab + shift) - OFFSET;
  nzlo_out = nlo + nlower;
  nzhi_out = nhi + nupper;



   // PPPM grid pts owned by this proc, including ghosts

   ngrid = (nxhi_out-nxlo_out+1) * (nyhi_out-nylo_out+1) *
     (nzhi_out-nzlo_out+1);

}


void FixCondiff::init()
{

  triclinic = domain->triclinic;


  int itmp = 0;


  // free all arrays previously allocated

  deallocate();

  // setup FFT grid resolution and g_ewald
  // normally one iteration thru while loop is all that is required
  // if grid stencil does not extend beyond neighbor proc
  //   or overlap is allowed, then done
  // else reduce order and try again

  int (*procneigh)[2] = comm->procneigh;

  GridComm *cgtmp = NULL;
  int iteration = 0;

  while (order >= minorder) {


      set_grid_local();

    cgtmp = new GridComm(lmp,world,1,1,
                         nxlo_in,nxhi_in,nylo_in,nyhi_in,nzlo_in,nzhi_in,
                         nxlo_out,nxhi_out,nylo_out,nyhi_out,nzlo_out,nzhi_out,
                         procneigh[0][0],procneigh[0][1],procneigh[1][0],
                         procneigh[1][1],procneigh[2][0],procneigh[2][1]);
    cgtmp->ghost_notify();
    if (!cgtmp->ghost_overlap()) break;
    delete cgtmp;

    order--;
    iteration++;
  }

  if (cgtmp) delete cgtmp;


  int ngrid_max,nfft_both_max;
  MPI_Allreduce(&ngrid,&ngrid_max,1,MPI_INT,MPI_MAX,world);
  MPI_Allreduce(&nfft_both,&nfft_both_max,1,MPI_INT,MPI_MAX,world);



  // allocate K-space dependent memory
  // don't invoke allocate peratom() or group(), will be allocated when needed

  allocate();
  cg->ghost_notify();
  cg->setup();

  // pre-compute Green's function denomiator expansion
  // pre-compute 1d charge distribution coefficients


  compute_rho_coeff();
}

void FixCondiff::deallocate()
{
  memory->destroy3d_offset(density_brick,nzlo_out,nylo_out,nxlo_out);


   memory->destroy3d_offset(vdx_brick,nzlo_out,nylo_out,nxlo_out);
   memory->destroy3d_offset(vdy_brick,nzlo_out,nylo_out,nxlo_out);
   memory->destroy3d_offset(vdz_brick,nzlo_out,nylo_out,nxlo_out);


  memory->destroy(density_fft);
  memory->destroy(greensfn);
  memory->destroy(work1);
  memory->destroy(work2);
  memory->destroy(vg);


    memory->destroy1d_offset(fkx,nxlo_fft);
    memory->destroy1d_offset(fky,nylo_fft);
    memory->destroy1d_offset(fkz,nzlo_fft);


  memory->destroy(gf_b);

  memory->destroy2d_offset(rho1d,-order_allocated/2);
  memory->destroy2d_offset(drho1d,-order_allocated/2);
  memory->destroy2d_offset(rho_coeff,(1-order_allocated)/2);
  memory->destroy2d_offset(drho_coeff,(1-order_allocated)/2);


}
void FixCondiff::allocate()
{
  memory->create3d_offset(density_brick,nzlo_out,nzhi_out,nylo_out,nyhi_out,
                          nxlo_out,nxhi_out,"pppm:density_brick");

  memory->create(density_fft,nfft_both,"pppm:density_fft");
  memory->create(greensfn,nfft_both,"pppm:greensfn");
  memory->create(work1,2*nfft_both,"pppm:work1");
  memory->create(work2,2*nfft_both,"pppm:work2");
  memory->create(vg,nfft_both,6,"pppm:vg");


    memory->create1d_offset(fkx,nxlo_fft,nxhi_fft,"pppm:fkx");
    memory->create1d_offset(fky,nylo_fft,nyhi_fft,"pppm:fky");
    memory->create1d_offset(fkz,nzlo_fft,nzhi_fft,"pppm:fkz");



    memory->create3d_offset(vdx_brick,nzlo_out,nzhi_out,nylo_out,nyhi_out,
                            nxlo_out,nxhi_out,"pppm:vdx_brick");
    memory->create3d_offset(vdy_brick,nzlo_out,nzhi_out,nylo_out,nyhi_out,
                            nxlo_out,nxhi_out,"pppm:vdy_brick");
    memory->create3d_offset(vdz_brick,nzlo_out,nzhi_out,nylo_out,nyhi_out,
                            nxlo_out,nxhi_out,"pppm:vdz_brick");


  // summation coeffs

  order_allocated = order;
  memory->create(gf_b,order,"pppm:gf_b");
  memory->create2d_offset(rho1d,3,-order/2,order/2,"pppm:rho1d");
  memory->create2d_offset(drho1d,3,-order/2,order/2,"pppm:drho1d");
  memory->create2d_offset(rho_coeff,order,(1-order)/2,order/2,"pppm:rho_coeff");
  memory->create2d_offset(drho_coeff,order,(1-order)/2,order/2,
                          "pppm:drho_coeff");


  // create ghost grid object for rho and electric field communication

  int (*procneigh)[2] = comm->procneigh;


    cg = new GridComm(lmp,world,3,1,
                      nxlo_in,nxhi_in,nylo_in,nyhi_in,nzlo_in,nzhi_in,
                      nxlo_out,nxhi_out,nylo_out,nyhi_out,nzlo_out,nzhi_out,
                      procneigh[0][0],procneigh[0][1],procneigh[1][0],
                      procneigh[1][1],procneigh[2][0],procneigh[2][1]);
}
