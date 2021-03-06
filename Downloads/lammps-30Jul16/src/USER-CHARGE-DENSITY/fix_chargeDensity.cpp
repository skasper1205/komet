#include "fix_chargeDensity.h"
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
FixChargeDensity::FixChargeDensity(LAMMPS *lmp, int narg, char **arg) :
		Fix(lmp, narg, arg)
{
	if (narg < 6) error->all(FLERR, "Illegal fix charge density command"); //5 mandatory arguments
	nevery = force->inumeric(FLERR,arg[5]);
        if (nevery <= 0) error->all(FLERR,"Illegal fix charge density command");

       	MPI_Comm_rank(world,&me);
       	MPI_Comm_size(world,&nprocs);
	
	masstotal = 0;
	
       	jgroup = group->find(arg[3]);
	if (jgroup == -1) error->all(FLERR,"Could not find fix chargeDensity group ID");
		kgroup = group->find(arg[4]);
	if (kgroup == -1) error->all(FLERR,"Could not find fix chargeDensity group ID");
	groupbit_coions = group->bitmask[jgroup];
	groupbit_counterions = group->bitmask[kgroup];

	file = fopen("charge_density","w");
	fclose(file);
	file = fopen("charge_density","a");

}

//The class destructor
FixChargeDensity::~FixChargeDensity()
{
   fclose(file);
}

//Where algorithm steps in
int FixChargeDensity::setmask()
{
	int mask = 0;
	mask |= END_OF_STEP;
	return mask;
}

void FixChargeDensity::end_of_step()
{
  if (update->ntimestep % nevery) return;
  calc_cm();
  calc_charge_density();
  delete [] cm;
}

void FixChargeDensity::calc_cm()
{
  cm = new double[3];
  masstotal = group->mass(igroup);
  group->xcm(igroup,masstotal,cm);
   //printf("Center of mass: %f, %f, %f", cm[0], cm[1], cm[2]);

}

void FixChargeDensity::calc_charge_density()
{

    double *prd;

  	prd = domain->prd;
  
  	double xprd = prd[0];
  	double yprd = prd[1];
  	double zprd = prd[2];
	
	double **x = atom->x;
	int nlocal = atom->nlocal;
	int *mask = atom->mask;
	double *q = atom->q;
	
	double R = 2.5;
	double r = 0.25;
	double counter1;
	double counter2;
	double V;

	double x_help;
	double y_help;
	double z_help;
	double rho1, rho2;
	
	while(R<=15){
	  counter1 = 0;
	  counter2 = 0;
	  V = (4./3.)*MY_PI*(pow((R+r),3)-pow(R,3));
	  for (int i = 0; i < nlocal; i++) {

	    if(mask[i] & groupbit_coions){
	      x_help = fabs(x[i][0]-cm[0]);
	      y_help = fabs(x[i][1]-cm[1]);
	      z_help = fabs(x[i][2]-cm[2]);
	      
	      if(fabs(x[i][0]-cm[0])>xprd-(R+r)){
		x_help=fabs(xprd-(x[i][0]-cm[0]));
	      }
	      if(fabs(x[i][1]-cm[1])>yprd-(R+r)){
		y_help=fabs(yprd-(x[i][1]-cm[1]));
	      }
	      if(fabs(x[i][2]-cm[2])>zprd-(R+r)){
		z_help=fabs(zprd-(x[i][2]-cm[2]));
	      }

		/*if(x_help > (xprd/2.)){
			x_help = xprd - x_help;
		}
		if((y_help > yprd/2.)){
			y_help = yprd - y_help;
		}
		if((z_help > zprd/2.)){
			z_help = zprd - z_help;
		}*/
	      
	      if(R <= sqrt(pow(x_help,2)+pow(y_help,2)+pow(z_help,2)) &&
		 sqrt(pow(x_help,2)+pow(y_help,2)+pow(z_help,2)) < (R+r)){
		counter1++;
	      }
	    }
	      
	    if(mask[i] & groupbit_counterions){
		x_help = fabs(x[i][0]-cm[0]);
		y_help = fabs(x[i][1]-cm[1]);
		z_help = fabs(x[i][2]-cm[2]);
	      
			      if(fabs(x[i][0]-cm[0])>xprd-(R+r)){
		x_help=fabs(xprd-(x[i][0]-cm[0]));
	      }
	      if(fabs(x[i][1]-cm[1])>yprd-(R+r)){
		y_help=fabs(yprd-(x[i][1]-cm[1]));
	      }
	      if(fabs(x[i][2]-cm[2])>zprd-(R+r)){
		z_help=fabs(zprd-(x[i][2]-cm[2]));
	      }

		/*if(x_help > (xprd/2.)){
			x_help = xprd - x_help;
		}
		if((y_help > yprd/2.)){
			y_help = yprd - y_help;
		}
		if((z_help > zprd/2.)){
			z_help = zprd - z_help;
		}*/

	      
		if(R <= sqrt(pow(x_help,2)+pow(y_help,2)+pow(z_help,2)) &&
		   sqrt(pow(x_help,2)+pow(y_help,2)+pow(z_help,2)) < (R+r)){
		  counter2++;
		}
	    }

	    
	  }
	  rho1 = counter1/V;
	  rho2 = counter2/V;
	  fprintf(file,"%f\t%f\t%f\t%f\n",R,V, rho1, rho2);
	  R+=r;
   	}
}
