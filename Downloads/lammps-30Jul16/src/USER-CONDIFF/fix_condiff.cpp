/*
 * fix_condiff.cpp
 *
 *  Created on: Oct 24, 2016
 *      Author: skasper
 */


#include "fix_condiff.h"
#include "stdio.h"
#include "stdlib.h"
#include "string.h"
#include "atom.h"
#include "force.h"
#include "update.h"
#include "respa.h"
#include "error.h"
#include "memory.h"
#include "input.h"
#include "modify.h"
#include "variable.h"
#include "pair.h"
#include "pppm.h"
#include "verlet.h"
#include "kspace.h"
#include "style_kspace.h"
#include "comm.h"
#include "gridcomm.h"
#include "integrate.h"
#include "min.h"
#include "compute.h"
#include <mpi.h>
#include <ctype.h>
#include <float.h>
#include <limits.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "pair.h"
#include "atom.h"
#include "neighbor.h"
#include "neigh_list.h"
#include "domain.h"
#include "comm.h"
#include "force.h"
#include "kspace.h"
#include "update.h"
#include "modify.h"
#include "compute.h"
#include "suffix.h"
#include "atom_masks.h"
#include "memory.h"
#include "math_const.h"
#include "error.h"


using namespace LAMMPS_NS;
using namespace FixConst;

//The class constructor.
FixCondiff::FixCondiff(LAMMPS *lmp, int narg, char **arg) :
		Fix(lmp, narg, arg)
{
	if (narg < 4) error->all(FLERR, "Illegal fix condiff command");

	/*int n = strlen(arg[3]) + 1; //saves 4th argument of script in a string
	string = new char[n];
	strcpy(string,arg[3]);*/

	//fp = fopen( "out_file.txt", "w" );
	//fprintf(fp, arg[3]);       //writes that argument in file
	//fclose(fp);

}

FixCondiff::~FixCondiff()
{

}

int FixCondiff::setmask() //where should algorithm step in
{
  int mask = 0;
  mask |= END_OF_STEP;
  return mask;
}

void FixCondiff::end_of_step()
{
	kspace_check();
	pppm_check();


	p = fopen( "new.txt", "w" );
	fprintf(p, "");
	fclose(p);

	f = fopen( "grid_info.txt", "r" );
	p = fopen( "new.txt", "a" );

	fgets(data,100,f);
	n=strtok(data,"=");
	for(int i=0; i<2; i++){
	 fprintf(p,"%s\n",n);
	 n=strtok(NULL, "=");
	 fprintf(p,"%s\n",n);
     n=strtok(NULL, "=");
     fprintf(p,"%s\n",n);
     n=strtok(NULL, "=");

	}







	//puts(data);


	fclose(f);
	fclose(p);


}

void FixCondiff::kspace_check(){   //check if pppm computation is used

	if (!force->kspace)
		error->all(FLERR, "Not using kspace computations");

}

void FixCondiff::pppm_check(){   //check if pppm computation is used

	if (!force->pair->pppmflag)
		error->all(FLERR, "Not using pppm computations");

}

