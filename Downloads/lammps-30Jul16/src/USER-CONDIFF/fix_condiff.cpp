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
#include "group.h"
#include "update.h"
#include "random_mars.h"
#include "compute.h"
#include "modify.h"

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

enum{NOBIAS,BIAS};
enum{REVERSE_RHO};

//The class constructor.
FixCondiff::FixCondiff(LAMMPS *lmp, int narg, char **arg) :
		Fix(lmp, narg, arg)
{
	if (narg < 5) error->all(FLERR, "Illegal fix condiff command");

	MPI_Comm_rank(world,&me);
	MPI_Comm_size(world,&nprocs);

	density_brick  = NULL;
	density_brick_counter = NULL;
	density_brick_force = NULL;

	help_v = NULL;
	help_f = NULL;
	help_x = NULL;

	rand = NULL;

	gf_b = NULL;
	rho1d = rho_coeff = drho1d = drho_coeff = NULL;

	cg = NULL;
	cg_peratom = NULL;

	part2grid = NULL;

	order =  force->kspace->order;
	minorder=2;
	order_allocated = order;

	jgroup = group->find(arg[3]);            //set grouptbit_condiff so fix knows which particles to act on
	if (jgroup == -1) error->all(FLERR,"Could not find fix group ID");
	groupbit_condiff = group->bitmask[jgroup];

	kspace_check();
	pppm_check();

	D = 1.0;
	//optinal args
	int iarg = 5;
	while (iarg < narg) {
		D = atof(arg[iarg]);
		iarg++;
	}
	printf("diffusion coefficient = %f\n",D);
	int seed = 11111 ;
	random = new RanMars(lmp,seed + comm->me);

	int n = strlen(id) + 6;
	id_temp = new char[n];
	strcpy(id_temp,id);
	strcat(id_temp,"_temp");

	kgroup = group->find(arg[4]);

	char **newarg = new char*[3];
	newarg[0] = id_temp;
	newarg[1] = group->names[kgroup];
	newarg[2] = (char *) "temp";
	modify->add_compute(3,newarg);
	delete [] newarg;
	tflag = 1;
	//printf("group = %s\n",newarg[1]);

}

FixCondiff::~FixCondiff()
{
	  deallocate();
	  memory->destroy(part2grid);
	  memory->destroy(density_brick);
	  memory->destroy(density_brick_counter);
	  memory->destroy(density_brick_force);
	  if (tflag) modify->delete_compute(id_temp);
}

int FixCondiff::setmask()	//where algorithm steps in
{
	int mask = 0;
	mask |= POST_FORCE;
	return mask;
}

void FixCondiff::post_force(int vspace)
{
	setup();
	setup_grid();
	particle_map();
    make_rho();

    cg->reverse_comm(force->kspace,REVERSE_RHO);

    reverse_make_rho();
    apply_boundary_conditions();

}

void FixCondiff::end_of_step(){

}


void FixCondiff::setup()
{
	nx_pppm = force->kspace->nx_pppm;
	ny_pppm = force->kspace->ny_pppm;
	nz_pppm = force->kspace->nz_pppm;

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
	dt = update->dt;
	wienerConst = sqrt(2*D*dt);


	T=1;

	/*int icompute = modify->find_compute(id_temp);
	if (icompute < 0)
	  error->all(FLERR,"Temperature ID for fix temp/berendsen does not exist");
	temperature = modify->compute[icompute];

	if (temperature->tempbias) which = BIAS;
	else which = NOBIAS;

	T = temperature->compute_scalar();*/

	//printf("boxhix = %f\n",boxhi[0]);

}

void FixCondiff::setup_grid()
{
  // free all arrays previously allocated

	deallocate();
  //if (peratom_allocate_flag) deallocate_peratom();
  //if (group_allocate_flag) deallocate_groups();

  // reset portion of global grid that each proc owns

	set_grid_local();

  // reallocate K-space dependent memory
  // check if grid communication is now overlapping if not allowed
  // don't invoke allocate peratom() or group(), will be allocated when needed

	allocate();

	cg->ghost_notify();

	cg->setup();

	compute_rho_coeff();


}
void FixCondiff::particle_map()
{
	int nlocal = atom->nlocal;

	if (atom->nmax > nmax) {
	    memory->destroy(part2grid);
	    nmax = atom->nmax;
	    memory->create(part2grid,nmax,3,"condiff:part2grid");
	}

	int flag;

	int nx,ny,nz;
    boxlo = domain->boxlo;

	double **x = atom->x;

	for (int i = 0; i < nlocal; i++) {

    // (nx,ny,nz) = global coords of grid pt to "lower left" of charge
    // current particle coord can be outside global and local box
    // add/subtract OFFSET to avoid int(-0.75) = 0 when want it to be -1

		nx = static_cast<int> ((x[i][0]-boxlo[0])*delxinv+shift) - OFFSET;
		ny = static_cast<int> ((x[i][1]-boxlo[1])*delyinv+shift) - OFFSET;
		nz = static_cast<int> ((x[i][2]-boxlo[2])*delzinv+shift) - OFFSET;

	    part2grid[i][0] = nx;
		part2grid[i][1] = ny;
		part2grid[i][2] = nz;

		//if (nx+nlower < nxlo_out || nx+nupper > nxhi_out ||
		//        ny+nlower < nylo_out || ny+nupper > nyhi_out ||
		//        nz+nlower < nzlo_out || nz+nupper > nzhi_out)
		//      flag = 1;
	}
	if (flag) error->one(FLERR,"Out of range atoms - cannot compute Condiff");
}

void FixCondiff::make_rho()
{
	int l,m,n,nx,ny,nz,mx,my,mz;
	FFT_SCALAR dx,dy,dz,x0,y0,z0;

  // clear 3d density array

	//memset(&(density_brick[3][nzlo_out][nylo_out][nxlo_out]),0,
    //     3*ngrid*sizeof(FFT_SCALAR));

  // loop over my velocities, add their contribution to nearby grid points
  // (nx,ny,nz) = global coords of grid pt to "lower left" of charge
  // (dx,dy,dz) = distance to "lower left" grid pt
  // (mx,my,mz) = global coords of moving stencil pt

	double **v = atom->v;
	double **x = atom->x;
	double **f = atom->f;
	int nlocal = atom->nlocal;
	int *mask = atom->mask;

	for (int i = 0; i < nlocal; i++) {

		nx = part2grid[i][0];
		ny = part2grid[i][1];
		nz = part2grid[i][2];
		dx = nx+shiftone - (x[i][0]-boxlo[0])*delxinv;
		dy = ny+shiftone - (x[i][1]-boxlo[1])*delyinv;
		dz = nz+shiftone - (x[i][2]-boxlo[2])*delzinv;

		compute_rho1d(dx,dy,dz);

		for(int j = 0; j < 3; j++){

			for (n = nlower; n <= nupper; n++) {
				mz = n+nz;
				for (m = nlower; m <= nupper; m++) {
					my = m+ny;
					for (l = nlower; l <= nupper; l++) {
						mx = l+nx;
						density_brick[j][mz][my][mx] = 0;
						density_brick_counter[j][mz][my][mx] = 0;
						density_brick_force[j][mz][my][mx] = 0;
					}
				}
			}
		}
	}

	for (int i = 0; i < nlocal; i++) { //take velocity of dpd-particles and map them on grid

		if (mask[i] & groupbit){
			nx = part2grid[i][0];
			ny = part2grid[i][1];
			nz = part2grid[i][2];
			dx = nx+shiftone - (x[i][0]-boxlo[0])*delxinv;
			dy = ny+shiftone - (x[i][1]-boxlo[1])*delyinv;
			dz = nz+shiftone - (x[i][2]-boxlo[2])*delzinv;

			compute_rho1d(dx,dy,dz);

			for (int j = 0; j < 3; j++) {
				//f[i][j] -= help_f[i][j];
				z0 = delvolinv;
				for (n = nlower; n <= nupper; n++) {
					mz = n+nz;
					y0 = z0*rho1d[2][n];
					for (m = nlower; m <= nupper; m++) {
						my = m+ny;
						x0 = y0*rho1d[1][m];
						for (l = nlower; l <= nupper; l++) {
							mx = l+nx;
							density_brick[j][mz][my][mx] += x0*rho1d[0][l]*v[i][j];
							density_brick_counter[j][mz][my][mx]+=x0*rho1d[0][l];

						}
					}
				}
			}
		}
		if (mask[i] & groupbit_condiff){ //take force of condiff-particles and map them on grid
			nx = part2grid[i][0];
			ny = part2grid[i][1];
			nz = part2grid[i][2];
			dx = nx+shiftone - (x[i][0]-boxlo[0])*delxinv;
			dy = ny+shiftone - (x[i][1]-boxlo[1])*delyinv;
			dz = nz+shiftone - (x[i][2]-boxlo[2])*delzinv;

			compute_rho1d(dx,dy,dz);

			for (int j = 0; j < 3; j++) {
				z0 = delvolinv;
				for (n = nlower; n <= nupper; n++) {
					mz = n+nz;
					y0 = z0*rho1d[2][n];
					for (m = nlower; m <= nupper; m++) {
						my = m+ny;
						x0 = y0*rho1d[1][m];
						for (l = nlower; l <= nupper; l++) {
							mx = l+nx;
							density_brick_force[j][mz][my][mx] += x0*rho1d[0][l]*f[i][j];
							//density_brick[j][mz][my][mx] += x0*rho1d[0][l]*v[i][j];
						}
					}
				}
			}
		}
	}
}

void FixCondiff::reverse_make_rho()
{

	int l,m,n,nx,ny,nz,mx,my,mz;
	FFT_SCALAR dx,dy,dz,x0,y0,z0;

  // clear 3d density array

	//memset(&(density_brick[3][nzlo_out][nylo_out][nxlo_out]),0,
    //     3*ngrid*sizeof(FFT_SCALAR));

  // loop over my velocities, add their contribution to nearby grid points
  // (nx,ny,nz) = global coords of grid pt to "lower left" of charge
  // (dx,dy,dz) = distance to "lower left" grid pt
  // (mx,my,mz) = global coords of moving stencil pt



	double **v = atom->v;
	double **x = atom->x;
	double **f = atom->f;
	int nlocal = atom->nlocal;
	int *mask = atom->mask;

	for (int i = 0; i < nlocal; i++) {

		if (mask[i] & groupbit_condiff){ //remap velocities on condiff-particles (pseudo-ions)
			nx = part2grid[i][0];
			ny = part2grid[i][1];
			nz = part2grid[i][2];
			dx = nx+shiftone - (x[i][0]-boxlo[0])*delxinv;
			dy = ny+shiftone - (x[i][1]-boxlo[1])*delyinv;
			dz = nz+shiftone - (x[i][2]-boxlo[2])*delzinv;

			compute_rho1d(dx,dy,dz);

			for(int j = 0; j < 3; j++){

				rand[j] = random->gaussian();

				v[i][j] = 0;
				for (n = nlower; n <= nupper; n++) {
					mz = n+nz;
					z0 = rho1d[2][n];
					for (m = nlower; m <= nupper; m++) {
						my = m+ny;
						y0 = z0*rho1d[1][m];
						for (l = nlower; l <= nupper; l++) {
							mx = l+nx;
							x0 = y0*rho1d[0][l];
							if(density_brick_counter[j][mz][my][mx] !=0) {
								v[i][j] += (density_brick[j][mz][my][mx]*x0/density_brick_counter[j][mz][my][mx]);
							}
						}
					}
				}
				x[i][j] += v[i][j]*dt + f[i][j]*dt*D/T + wienerConst*rand[j];     //euler step //Temperatur soll einfach konstant uebergeben werden, nicht imer ausgerechnet
				//printf("%f\n",rand[j]);
			}
		}
		if (mask[i] & groupbit){
			nx = part2grid[i][0];
			ny = part2grid[i][1];
			nz = part2grid[i][2];
			dx = nx+shiftone - (x[i][0]-boxlo[0])*delxinv;
			dy = ny+shiftone - (x[i][1]-boxlo[1])*delyinv;
			dz = nz+shiftone - (x[i][2]-boxlo[2])*delzinv;

			compute_rho1d(dx,dy,dz);

			for(int j = 0; j < 3; j++){
				//f[i][j] -= help_f[i][j];
				help_f[i][j] = 0;
				for (n = nlower; n <= nupper; n++) {
					mz = n+nz;
					z0 = rho1d[2][n];
					for (m = nlower; m <= nupper; m++) {
						my = m+ny;
						y0 = z0*rho1d[1][m];
						for (l = nlower; l <= nupper; l++) {
							mx = l+nx;
							x0 = y0*rho1d[0][l];
							if(density_brick_counter[j][mz][my][mx] !=0){
								help_f[i][j] += (density_brick_force[j][mz][my][mx]*x0/density_brick_counter[j][mz][my][mx]);
							}
						}
					}
				}
				f[i][j] += help_f[i][j];
			}
		}
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
		for (l = 0; l < order; l++){
			rho_coeff[l][m] = a[l][k];
		}
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


  	if (force->kspace->stagger_flag) {
  	    nxhi_out++;
  	    nyhi_out++;
  	    nzhi_out++;
  	  }
   // PPPM grid pts owned by this proc, including ghosts

  	ngrid = (nxhi_out-nxlo_out+1) * (nyhi_out-nylo_out+1) *
     (nzhi_out-nzlo_out+1);

}



void FixCondiff::deallocate()
{
	memory->destroy4d_offset(density_brick,nzlo_out,nylo_out,nxlo_out);
	memory->destroy4d_offset(density_brick_counter,nzlo_out,nylo_out,nxlo_out);
	memory->destroy4d_offset(density_brick_force,nzlo_out,nylo_out,nxlo_out);

	memory->destroy2d_offset(rho1d,-order_allocated/2);
	memory->destroy2d_offset(drho1d,-order_allocated/2);
	memory->destroy2d_offset(rho_coeff,(1-order_allocated)/2);
	memory->destroy2d_offset(drho_coeff,(1-order_allocated)/2);

	memory->destroy(help_v);
	memory->destroy(help_f);
	memory->destroy(help_x);
	memory->destroy(rand);

	delete cg;

}
void FixCondiff::allocate()
{
	memory->create4d_offset(density_brick,3,nzlo_out,nzhi_out,nylo_out,nyhi_out,
                          nxlo_out,nxhi_out,"condiff:density_brick");
	memory->create4d_offset(density_brick_counter,3,nzlo_out,nzhi_out,nylo_out,nyhi_out,
	                          nxlo_out,nxhi_out,"condiff:density_brick");
	memory->create4d_offset(density_brick_force,3,nzlo_out,nzhi_out,nylo_out,nyhi_out,
		                          nxlo_out,nxhi_out,"condiff:density_brick");

	order_allocated = order;
	memory->create2d_offset(rho1d,3,-order/2,order/2,"condiff:rho1d");
	memory->create2d_offset(drho1d,3,-order/2,order/2,"condiff:drho1d");
	memory->create2d_offset(rho_coeff,order,(1-order)/2,order/2,"condiff:rho_coeff");
	memory->create2d_offset(drho_coeff,order,(1-order)/2,order/2,
                          "condiff:drho_coeff");

	memory->create(help_v,atom->nmax,3,"condiff:help_v");
	memory->create(help_f,atom->nmax,3,"condiff:help_f");
	memory->create(help_x,atom->nmax,3,"condiff:help_x");

	memory->create(rand,3,"condiff:rand");
  // create ghost grid object for rho and electric field communication

	int (*procneigh)[2] = comm->procneigh;


    cg = new GridComm(lmp,world,3,1,
                      nxlo_in,nxhi_in,nylo_in,nyhi_in,nzlo_in,nzhi_in,
                      nxlo_out,nxhi_out,nylo_out,nyhi_out,nzlo_out,nzhi_out,
                      procneigh[0][0],procneigh[0][1],procneigh[1][0],
                      procneigh[1][1],procneigh[2][0],procneigh[2][1]);
}

void FixCondiff::kspace_check()	//check if pppm computation is used
{

	if (!force->kspace)
		error->all(FLERR, "Not using kspace computations");

}

void FixCondiff::pppm_check()	//check if pppm computation is used
{

	if (!force->pair->pppmflag)
		error->all(FLERR, "Not using pppm computations");

}

void FixCondiff::apply_boundary_conditions()
{
	double **x = atom->x;
	int nlocal = atom->nlocal;
	int *mask = atom->mask;

	double *boxlo, *boxhi;
	boxlo = domain->boxlo;
	boxhi = domain->boxhi;

	double *prd;
	prd = domain->prd;

	for (int i = 0; i < nlocal; i++) {

		if(mask[i] & groupbit_condiff){
			while(x[i][0] > boxhi[0]) x[i][0] -= prd[0];
			while(x[i][1] > boxhi[1]) x[i][1] -= prd[1];
			while(x[i][2] > boxhi[2]) x[i][2] -= prd[2];

			while(x[i][0] < boxlo[0]) x[i][0] += prd[0];
			while(x[i][1] < boxlo[1]) x[i][1] += prd[1];
			while(x[i][2] < boxlo[2]) x[i][2] += prd[2];
		}
	}
}
