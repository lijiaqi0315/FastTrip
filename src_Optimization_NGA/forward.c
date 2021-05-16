#include	"forward.h"

int eval_pop(pop *p)
{
	if(cflag==1){
	    return(eval_pop_mahfoud(p));
	} else if(cflag==2){
	    return(eval_pop_fortran(p));
	} else {
	    fprintf(stderr,"Unknown choice of objective function\n");
	    return(1);
	}
}

/*
* Evaulates each model, and then sorts each deme
* based on cost.
*/
int eval_pop_mahfoud(pop *p)
{
	int i,j;

	for(i=0;i<ndeme;i++){
	    for(j=0;j<nmod;j++){
		p->d[i].m[j].cost=mahfoud(p->d[i].m[j].par[0]);
	    }
	}

	if(sort_pop(p)==1){
	    fprintf(stderr,"Bad sorting\n");
	    return(1);
	}

        return(0);
}

/*
* Test objective function
*/
float mahfoud(float x)
{
        float tmp1,tmp2,phi;

        tmp1=-2.0*log(2.0)*((x-.1)/.8);
        tmp1*=((x-.1)/.8);
        tmp2=sin(5.0*mypi*x);
        phi=exp(tmp1)*pow(tmp2,6.0);
        phi=1.0-phi;

        return(phi);
}

int eval_pop_fortran(pop *p)
{
	int i,j,k,lcnt;
	float cst;
	char *iname,*oname,*script,*line;
	FILE *fpin,*fpout;

	line=kalloc(200);

	iname=".fortran_input_file";
	oname=".fortran_output_file";
	script="sh go_fort";

/*Open fortran input file*/
	if((fpin=fopen(iname,"w"))==NULL){
	    fprintf(stderr,"Can't write to %s\n",iname);
	    return(1);
	}

/*Construct fortran input file*/
	fprintf(fpin,"%d %d %d\n",npar,nmod,ndeme);
	for(i=0;i<ndeme;i++){
	    for(j=0;j<nmod;j++){
		for(k=0;k<npar;k++){
		    fprintf(fpin,"%10.3e\n",p->d[i].m[j].par[k]);
		}
	    }
	}
	fclose(fpin);

/*Execute shell command that does fortran program*/
	system(script);

/*Open fortran output file*/
	if((fpout=fopen(oname,"r"))==NULL){
	    fprintf(stderr,"Can't read %s\n",oname);
	    return(1);
	}

/*Read cost values from fortran output and check number of lines*/
	lcnt=0; i=0; j=0;
	while(fgets(line,200*sizeof(char),fpout)!=NULL){
	    sscanf(line,"%f",&cst);
	    p->d[j].m[i].cost=cst;
	      i++;
	      if(i==nmod){
		  j++;
		  i=0;
	      }
	      lcnt++;
	}
	fclose(fpout);

	if(lcnt!=nmod*ndeme){
	    fprintf(stderr,"Problem reading fortran output\n");
	    return(1);
	}

	free(line);
	return(0);
}
