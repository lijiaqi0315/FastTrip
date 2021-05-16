#include	<math.h>
#include	<stdlib.h>
#include	<malloc.h>
#include	"random.h"

#define MX 100

typedef struct{
	float *par,cost,fit;
} model;

typedef struct{
	model *m;
} deme;

typedef struct{
	deme *d;
} pop;

char* parse_line(FILE*,int);
char* kalloc(int);
float* falloc(int);
void print_pop(pop,int);

int get_info(char*,int);
int init_pop(pop*);
int shuffle(pop*);
int sort(pop*);
int mutate(pop*);
int xover(pop*);
int get_fit(pop*);
int select_pop(pop*);
int reprod(pop*);
int do_elitism(pop*);
int load_elite(pop*);
int compete(pop*);
float get_dist(model,model);

float bot[MX],top[MX],dev[MX],d_weight[MX],pc,pm,Rc;
int cflag,dflag,ndeme,nmod,nelite,ngen,npar,age;
char *Stype,*Mtype,*Xtype,*Ftype,*Dtype;
long seed;
pop Epop;
