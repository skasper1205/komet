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
	void end_of_step();
	void kspace_check();
    void pppm_check();

    FILE *f;
    FILE *p;
    char data[100];
    char * n;
    //int **array;
    char *string;
    int temp;

  private:
    int me;

};


}


#endif
#endif

