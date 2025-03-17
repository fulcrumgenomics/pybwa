#include "bwa.h"


// Sets the verbosity for the bwa C-API
int set_bwa_c_verbosity(int level)
{
	//extern int bwa_verbose;
	fprintf(stderr, "cur: %d new: %d\n", bwa_verbose, level);
	int retval = level == bwa_verbose ? 0 : 1;
	bwa_verbose = level;
	return retval;
}
