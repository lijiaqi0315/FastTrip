#include	<stdio.h>
#include	"ga_oper.h"

main()
{
	pop p;

	if(get_info("fga.in",25)==1){
	    fprintf(stderr,"Problem reading input file\n");
	    exit(1);
	}

	if(init_pop(&p)==1){
	    fprintf(stderr,"Problem initializing population\n");
	    exit(1);
	}

        fprintf(stderr,"\nGeneration: ");
	for(age=0;age<ngen;age++){
            fprintf(stderr," %d",age);
	    if(reprod(&p)==1){
		fprintf(stderr,"problem with reproduction\n");
		exit(1);
	    }
	    print_pop(p,age);
	    if(xover(&p)==1){
		fprintf(stderr,"problem with xover\n");
		exit(1);
	    }
	    if(mutate(&p)==1){
		fprintf(stderr,"problem with mutation\n");
		exit(1);
	    }
        }
        fprintf(stderr,"\n");
	exit(0);
}
