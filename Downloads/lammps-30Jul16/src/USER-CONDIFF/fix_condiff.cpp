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

using namespace LAMMPS_NS;
using namespace FixConst;

//The class constructor.
FixCondiff::FixCondiff(LAMMPS *lmp, int narg, char **arg) :
		Fix(lmp, narg, arg)
{
	fp = fopen( "out_file.txt", "w" );
	fprintf(fp, "Constructor");
	fclose(fp);
}

FixCondiff::~FixCondiff()
{

}
int FixCondiff::setmask()
{
  int mask = 0;
  mask |= END_OF_STEP;
  return mask;
}



