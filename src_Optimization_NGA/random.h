#include <sys/time.h>
#include <stdio.h>
#include <math.h>
 
#define MBIG  1000000000
#define MSEED 161803398
#define FAC   (1.0 / MBIG)
 
long seed_random(long);
float knuth_random();
float get_gauss();
