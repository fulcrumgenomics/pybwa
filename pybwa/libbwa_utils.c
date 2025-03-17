#include <stdio.h>

// Sets the verbosity for the bwa C-API
extern int bwa_verbose;
int set_bwa_c_verbosity(int level)
{
	fprintf(stderr, "cur: %d new: %d\n", bwa_verbose, level);
	int retval = level == bwa_verbose ? 0 : 1;
	bwa_verbose = level;
	return retval;
}
