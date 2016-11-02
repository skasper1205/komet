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

#include "fix.h"

namespace LAMMPS_NS {

class FixCondiff : public Fix {

  public:
	FixCondiff(class LAMMPS *, int, char **);
	~FixCondiff();
	int setmask();
	void post_force();
	void kspace_check();
    void pppm_check();

    int nx_pppm, ny_pppm, nz_pppm;
    int order;

    FILE *fp;

};


}


#endif
#endif

