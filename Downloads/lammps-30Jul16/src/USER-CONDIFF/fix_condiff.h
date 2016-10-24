/*
 * fix_condiff.h
 *
 *  Created on: Oct 24, 2016
 *      Author: skasper
 */

#ifdef FIX_CLASS

FixStyle(condiff, FixCondiff) //This registers this fix class with LAMMPS

#else

#ifndef LMP_FIX_CONDIFF_H     //  These are header guards.
#define LMP_FIX_CONDIFF_H

#include "fix.h"          // Must have this because our FixCONDIFF is derived from Fix...

namespace LAMMPS_NS {

class FixCondiff : public Fix {
  // The contents of your class go here i.e. member variables and methods etc...
  public:
	FixCondiff(class LAMMPS *, int, char **);
	~FixCondiff();
	int setmask();

	FILE *fp;


};


}


#endif
#endif

