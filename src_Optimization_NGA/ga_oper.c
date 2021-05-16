#include	"ga_oper.h"

void print_pop(pop p,int g)
{
	int i,j,k;
	char *fname;
	FILE *fp;
	
	fname=kalloc(200);
	sprintf(fname,"gen_sum_%03d",g);

	if((fp=fopen(fname,"w"))==NULL){
	    fprintf(stderr,"Can't write to %s\n",fname);
	    exit(1);
	}

	for(i=0;i<ndeme;i++){
	    fprintf(fp,"Deme %02d at generation %03d :\n",i,g);
	    for(j=0;j<nmod;j++){
		fprintf(fp," %03d %10.4e",j,p.d[i].m[j].cost);
		for(k=0;k<npar;k++){
		    fprintf(fp," %10.3e",p.d[i].m[j].par[k]);
	        }
		fprintf(fp,"\n");
	    }
	}
	free(fname);
	fclose(fp);
}

/*
* Mutate parameters using a Gaussian or uniform distribution,
* modulated by fractional standard deviations. If this produces
* a parameter value outside of the a'priori range then 
* then a new value is generated on the interval between the 
* violated  bound and the original value. This helps the
* ga come up with solutions close to the search boundaries,
* where most of the "volume" is for high-dimensional spaces.
*/
int mutate(pop *p)
{
	int i,j,k;
	float nval;

	if(strcmp(Mtype,"uniform")==0){
	  for(i=0;i<ndeme;i++){
	    for(j=0;j<nmod;j++){
	      for(k=0;k<npar;k++){
		if(knuth_random()<pm){
		  nval=p->d[i].m[j].par[k]+
		       dev[k]*(float)(2.0*knuth_random()-1.0);
		  if(nval<bot[k]){ 
		      p->d[i].m[j].par[k]=bot[k]+
		      (p->d[i].m[j].par[k]-bot[k])*knuth_random();
		  } else if(nval>top[k]){
		      p->d[i].m[j].par[k]=p->d[i].m[j].par[k]+
		      (top[k]-p->d[i].m[j].par[k])*knuth_random();
		  } else {
		      p->d[i].m[j].par[k]=nval;
		  }
		}
	      }
	    }
	  }
	  return(0);
	} else if(strcmp(Mtype,"gaussian")==0){
	  for(i=0;i<ndeme;i++){
	    for(j=0;j<nmod;j++){
	      for(k=0;k<npar;k++){
		if(knuth_random()<pm){
		  nval=p->d[i].m[j].par[k]+dev[k]*(float)get_gauss();
		  if(nval<bot[k]){ 
		      p->d[i].m[j].par[k]=bot[k]+
		      (p->d[i].m[j].par[k]-bot[k])*knuth_random();
		  } else if(nval>top[k]){
		      p->d[i].m[j].par[k]=p->d[i].m[j].par[k]+
		      (top[k]-p->d[i].m[j].par[k])*knuth_random();
		  } else {
		      p->d[i].m[j].par[k]=nval;
		  }
		}
	      }
	    }
	  }
	  return(0);
	} else {
	  fprintf(stderr,"Unknown mutation type: %s\n",Mtype);
	  return(1);
	}
}

int xover(pop *p)
{
	int i,j,k,b;
	float *tmod,*tmod2,a;

	tmod=falloc(npar);
	tmod2=falloc(npar);

	if(shuffle(p)==1){
	    fprintf(stderr,"Problem shuffling population\n");
	    return(1);
	}

	if(strcmp(Xtype,"spoint")==0){
	  for(i=0;i<ndeme;i++){
	    for(j=0;j<nmod;j+=2){
	      if(knuth_random()<pc){
	        b=(int)(knuth_random()*(float)(npar-1))+1;
		for(k=0;k<b;k++){
		  tmod[k]=p->d[i].m[j].par[k];
		  p->d[i].m[j].par[k]=p->d[i].m[j+1].par[k];
		  p->d[i].m[j+1].par[k]=tmod[k];
		}
	      }
	    }
	  }
	} else if(strcmp(Xtype,"lcombo")==0){
	  for(i=0;i<ndeme;i++){
	    for(j=0;j<nmod;j+=2){
	      if(knuth_random()<pc){
	        a=knuth_random();
		for(k=0;k<npar;k++){
		  tmod[k]=p->d[i].m[j].par[k];
		  tmod2[k]=p->d[i].m[j+1].par[k];
		  p->d[i].m[j].par[k]=a*tmod[k]+(1.0-a)*tmod2[k];
		  p->d[i].m[j+1].par[k]=(1.0-a)*tmod[k]+a*tmod2[k];
		}
	      }
	    }
	  }
	} else {
	    fprintf(stderr,"Unknown xover type: %s\n",Xtype);
	    return(1);
	}

	free(tmod); free(tmod2);
	return(0);
}

int reprod(pop *p)
{
	int i,j;

	if(eval_pop(p)==1){
	    fprintf(stderr,"Problem evaluating population\n");
	    return(1);
	}
	if(do_elitism(p)==1){
	    fprintf(stderr,"Problem adding elites to demes\n");
	    return(1);
	}
	if(sort_pop(p)==1){
	    fprintf(stderr,"Problem resorting population\n");
	    return(1);
	}
	if(load_elite(p)==1){
	    fprintf(stderr,"Problem loading Epop members\n");
	    return(1);
	}
	if(compete(p)==1){
	    fprintf(stderr,"Problem implementing niching\n");
	    return(1);
	}
	if(sort_pop(p)==1){
	    fprintf(stderr,"Problem re-resorting population\n");
	    return(1);
	}
	if(get_fit(p)==1){
	    fprintf(stderr,"Problem going from costs to fitnesses\n");
	    return(1);
	}
	if(select_pop(p)==1){
	    fprintf(stderr,"Problem doing selection\n");
	    return(1);
	}

	return(0);
}

/*
* Determines model progeny based on fitnesses
*/
int select_pop(pop *p)
{
	int i,j,k,l,out,in;
	float *R,*cnt,*len,sp;
	pop ptmp;
/*
* Allocate memory
*/
	R=falloc(nmod);
	cnt=falloc(nmod);
	len=falloc(nmod+1);

	if((ptmp.d=(deme*)malloc(ndeme*sizeof(deme)))==NULL){
            fprintf(stderr,"Can't allocate memory\n");
            return(1);
        }
        for(i=0;i<ndeme;i++){
            if((ptmp.d[i].m=(model*)malloc(nmod*sizeof(model)))==NULL){
                fprintf(stderr,"Can't allocate memory\n");
                return(1);
            }
            for(j=0;j<nmod;j++){
                if((ptmp.d[i].m[j].par=(float*)malloc(npar*sizeof(float)))==NULL){
                    fprintf(stderr,"Can't allocate memory\n");
                    return(1);
                }
            }
        }

	for(i=0;i<ndeme;i++){
	    for(j=0;j<nmod;j++){
	        cnt[j]=0.0;
	        R[j]=knuth_random();
	    }

	    if(strcmp("roulette",Stype)==0){
		;;
	    } else if(strcmp("universal",Stype)==0){
		sp=1.0/(float)nmod;
	        for(j=1;j<nmod;j++){
		    R[j]=R[j-1]+sp;
		    if(R[j]>1.0){
			R[j]-=1.0;
		    }
		}
	    } else {
		fprintf(stderr,"Unrecognized selection type: %s\n",Stype);
		return(1);
	    }

	    len[0]=0.0;
	    for(j=1;j<nmod;j++){
	        len[j]=len[j-1]+p->d[0].m[j].fit;
	    }
	    len[nmod]=1.0;

	    for(k=0;k<nmod;k++){
                for(j=0;j<nmod;j++){
                  if((R[k]>len[j])&&(R[k]<=len[j+1])){
                    cnt[j]+=1.0;
		    j=nmod;
                  }
                }
            }

	    out=0;in=0;
            for(l=0;l<nmod;l++){
              for(j=0;j<(int)(cnt[l]+0.5);j++){
                for(k=0;k<npar;k++){
		  ptmp.d[i].m[out+j].par[k]=
		  p->d[i].m[l].par[k];
                }
		ptmp.d[i].m[out+j].cost=p->d[i].m[l].cost;
		ptmp.d[i].m[out+j].fit=p->d[i].m[l].fit;
                in++;
              }
              out=in;
            }

	    for(j=0;j<nmod;j++){
		for(k=0;k<npar;k++){
		  p->d[i].m[j].par[k]=
		  ptmp.d[i].m[j].par[k];
	 	}
		p->d[i].m[j].cost=ptmp.d[i].m[j].cost;
		p->d[i].m[j].fit=ptmp.d[i].m[j].fit;
	    }
	}

/*
* Free ptmp, backwards from way allocated
*/
        for(i=0;i<ndeme;i++){
            for(j=0;j<nmod;j++){
		free(ptmp.d[i].m[j].par);
	    }
	    free(ptmp.d[i].m);
	}
	free(ptmp.d);

	free(R);
	free(cnt);
	free(len);
	return(0);
}

/*
* Maps costs to fitnesses (probablitites)
*/
int get_fit(pop *p)
{
	int i,j;
	float m,sum,max,min,ave,cmult,a,b;


	if(strcmp(Ftype,"rank")==0){
	    m=2.0/(float)(nmod*(nmod-1));
	    for(i=0;i<ndeme;i++){
		for(j=0;j<nmod;j++){
		    p->d[i].m[j].fit=m*(float)(nmod-j);
		}
	    }
	} else if(strcmp(Ftype,"basic")==0){
	    for(i=0;i<ndeme;i++){
		max=0.0;
		for(j=0;j<nmod;j++){
		    if(p->d[i].m[j].cost>max){
		        max=p->d[i].m[j].cost;
		    }
		}
	        sum=0.0;
		for(j=0;j<nmod;j++){
		    p->d[i].m[j].fit=max-p->d[i].m[j].cost;
		    sum+=p->d[i].m[j].fit;
		}
		for(j=0;j<nmod;j++){
		    p->d[i].m[j].fit/=sum;
		}
	    }
	} else if(strcmp(Ftype,"linear")==0){
	    cmult=2.0;
            for(i=0;i<ndeme;i++){
                max=0.0;
                for(j=0;j<nmod;j++){
                    if(p->d[i].m[j].cost>max){
                        max=p->d[i].m[j].cost;
                    }
                }
                sum=0.0;
                for(j=0;j<nmod;j++){
                    p->d[i].m[j].fit=max-p->d[i].m[j].cost;
		    sum+=p->d[i].m[j].fit;
                }
                for(j=0;j<nmod;j++){
                    p->d[i].m[j].fit/=sum;
                }
	        min=1.0e10; max=-1.0e10; ave=0.0;
                for(j=0;j<nmod;j++){
                    if(p->d[i].m[j].fit>max){ max=p->d[i].m[j].fit; }
                    if(p->d[i].m[j].fit<min){ min=p->d[i].m[j].fit; }
                    ave+=p->d[i].m[j].fit;
                }
		ave/=(float)nmod;
		if(min>((cmult*ave-max)/(cmult-1.0))){
                    a=(cmult-1.0)*ave/(max-ave);
                    b=ave*(max-cmult*ave)/(max-ave);
                } else {
                    a=ave/(ave-min);
                    b=-min*ave/(ave-min);
                }
                for(j=0;j<nmod;j++){
		    p->d[i].m[j].fit=a*p->d[i].m[j].fit+b;
                }
            }
	} else {
	    fprintf(stderr,"Unknown fitness type: %s\n",Ftype);
	    return(1);
	}

	return(0);
}

int shuffle(pop *p)
{
        int i,j,k,l;
        float *rand,*tmod,tmp;

        rand=falloc(nmod);
        tmod=falloc(npar);

/*Generates random numbers for piggybacking in sorting*/
        for(i=0;i<nmod;i++){
            rand[i]=knuth_random();
        }

/*
* sorts (low to high) array of RN's and array of models piggybacks,
* so it is in effect randomly shuffled. This sorting routine is way sucky.
*/
        for(l=0;l<ndeme;l++){
          for(i=0;i<(nmod-1);++i){
            for(j=(nmod-1);j>i;--j){
              if(rand[j-1]>rand[j]){
		tmp=rand[j];
		rand[j]=rand[j-1];
		rand[j-1]=tmp;
                for(k=0;k<npar;k++){
		    tmod[k]=p->d[l].m[j].par[k];
		    p->d[l].m[j].par[k]=p->d[l].m[j-1].par[k];
		    p->d[l].m[j-1].par[k]=tmod[k];
                }
              }
            }
          }
        }
	free(rand);
	free(tmod);
	return(0);
}

/*
* sort demes by cost
*/
int sort_pop(pop *p)
{
	int i,j,k,l;
	float tmp,*tmod;

	for(i=0;i<ndeme;i++){
	    for(j=0;j<nmod;j++){
		if(!(p->d[i].m[j].cost>-1.0e-10)){
		    fprintf(stderr,"Bad cost value (%e), can't sort\n",
			p->d[i].m[j].cost);
		    fprintf(stderr,"Deme: %02d, Model %03d\n",i,j);
		    return(1);
		}
	    }
	}

	tmod=falloc(npar);

        for(l=0;l<ndeme;l++){
          for(i=0;i<(nmod-1);++i){
            for(j=(nmod-1);j>i;--j){
              if(p->d[l].m[j-1].cost>p->d[l].m[j].cost){
                tmp=p->d[l].m[j].cost;
                p->d[l].m[j].cost=p->d[l].m[j-1].cost;
                p->d[l].m[j-1].cost=tmp;
                for(k=0;k<npar;k++){
                    tmod[k]=p->d[l].m[j].par[k];
                    p->d[l].m[j].par[k]=p->d[l].m[j-1].par[k];
                    p->d[l].m[j-1].par[k]=tmod[k];
                }
              }
            }
          }
        }

	free(tmod);
	return(0);
}

int init_pop(pop *p)
{
	int i,j,k;

	if((p->d=(deme*)malloc(ndeme*sizeof(deme)))==NULL){
	    fprintf(stderr,"Can't allocate memory\n");
	    return(1);
	}

	for(i=0;i<ndeme;i++){
	    if((p->d[i].m=(model*)malloc(nmod*sizeof(model)))==NULL){
	        fprintf(stderr,"Can't allocate memory\n");
		return(1);
	    }
	    for(j=0;j<nmod;j++){
	        if((p->d[i].m[j].par=(float*)malloc(npar*sizeof(float)))==NULL){
	            fprintf(stderr,"Can't allocate memory\n");
		    return(1);
	        }
	    }
	}

	for(i=0;i<ndeme;i++){
	    for(j=0;j<nmod;j++){
		for(k=0;k<npar;k++){
		    p->d[i].m[j].par[k]=
			bot[k]+(top[k]-bot[k])*knuth_random();
		}
	    }
	}
	return(0);
}

/*
* Alter demes so that elitist models survive automatically.
* Elitist models for each deme are stored in the global Epop
* population. Replaces worst guys in each deme.
*/
int do_elitism(pop *p)
{
	int i,j,k,cnt,flag;
/*
* Check elitist population
*/
	flag=0;
	for(i=0;i<ndeme;i++){
	    for(j=0;j<nelite;j++){
		if(Epop.d[i].m[j].cost<0.0){
		    flag++;
		}
	    }
	}

	if(flag==(ndeme*nelite)){
/*	   fprintf(stderr,"Skipping elitism of initial demes\n");*/
	   return(0);
	}
	
	for(i=0;i<ndeme;i++){
	    for(j=1;j<nelite;j++){
		if(Epop.d[i].m[j].cost<Epop.d[i].m[j-1].cost){
		    fprintf(stderr,"Problem with elitist population\n");
		    return(1);
		}
	    }
	}

	for(i=0;i<ndeme;i++){
	    cnt=0;
	    for(j=(nmod-1);j>(nmod-1-nelite);j--){
		p->d[i].m[j].cost=Epop.d[i].m[nelite-1-cnt].cost;
		p->d[i].m[j].fit=Epop.d[i].m[nelite-1-cnt].fit;
		for(k=0;k<npar;k++){
		    p->d[i].m[j].par[k]=Epop.d[i].m[nelite-1-cnt].par[k];
		}
		cnt++;
	    }
	}
	return(0);
}

/*
* Load elitist guys into global Epop.
*/
int load_elite(pop *p)
{
        int i,j,k,cnt;
/*
* Load elitist population
*/
        for(i=0;i<ndeme;i++){
            for(j=0;j<nelite;j++){
                Epop.d[i].m[j].cost=p->d[i].m[j].cost;
                Epop.d[i].m[j].fit=p->d[i].m[j].fit;
		for(k=0;k<npar;k++){
                    Epop.d[i].m[j].par[k]=p->d[i].m[j].par[k];
                }
            }
        }
	return(0);
}

/*
* Introduce competition between demes, does not
* assume each deme is sorted. Uses worst cost of
* all demes as penalty for being to similar to a lead
* model. Don't qorry about fitnesses of models because
* they will be recalculated in next step. Just keep track
* of parameters and costs.
*/
int compete(pop *p)
{
	int i,j,k,l;
	float min,max,penal_cost,dist;
	model *bestm;
/*
* Allocate array of top models
*/
	bestm=(model*)malloc(ndeme*sizeof(model));
	for(i=0;i<ndeme;i++){
	    bestm[i].par=falloc(npar);
	}
/*
* Get penalty cost for bad models
*/
	penal_cost=0.0;
	for(i=0;i<ndeme;i++){
	    for(j=0;j<nmod;j++){
		if(p->d[i].m[j].cost>penal_cost){
		    penal_cost=p->d[i].m[j].cost;
		}
	    }
	}

/*
* Find alpha model
*/
	min=1.0e10;
	for(j=0;j<nmod;j++){
	    if(p->d[0].m[j].cost<min){
		min=p->d[0].m[j].cost;
	        bestm[0].cost=min;
		for(k=0;k<npar;k++){
	            bestm[0].par[k]=p->d[0].m[j].par[k];
		}
	    }
	}
/*
* Start deme competition with second deme. Define
* beta, gamma, etc. after competition against lower
* order demes has been carried out.
*/
	for(i=1;i<ndeme;i++){
	    for(j=0;j<nmod;j++){
	        for(k=0;k<i;k++){
		    if((dist=get_dist(bestm[k],p->d[i].m[j]))<0.0){
			fprintf(stderr,"problem computing distance\n");
		 	return(1);
		    }	
		    if(dist<Rc){
			p->d[i].m[j].cost=penal_cost;
/*
		        printf("Deme %d, model %d is penalized to %f because of\n"
			       "distance of %f w.r.t to lead model of deme %d\n",
				i,j,penal_cost,dist,k);
			for(l=0;l<npar;l++){
			    printf(" %f %f\n",p->d[i].m[j].par[l],bestm[k].par[l]);
			}
*/
		    }
		}
	    }
	    min=1.0e10;
	    for(j=0;j<nmod;j++){
	        if(p->d[i].m[j].cost<min){
		    min=p->d[i].m[j].cost;
	            bestm[i].cost=min;
		    for(k=0;k<npar;k++){
	                bestm[i].par[k]=p->d[i].m[j].par[k];
		    }
	        }
	    }
	}
	for(i=0;i<ndeme;i++){
	    free(bestm[i].par);
	}
	free(bestm);
	return(0);
}

/*
* Normal distance metric, parameters weights are supplied
* by user. If weight is 0.0 then that parameter is not
* included in distance calculation, and does not affect
* the normalization.
*/
float get_dist(model m1,model m2)
{
	int i,rpar;
	float dist;

	rpar=0;
	if(strcmp("basic",Dtype)==0){
	    dist=0.0;
	    for(i=0;i<npar;i++){
		if(d_weight[i]>1.0e-6){
	            dist+=d_weight[i]*fabs(m2.par[i]-m1.par[i])/(top[i]-bot[i]);
		    rpar++;
		}
	    }
	    dist/=(float)rpar;
	    return(dist);
	} else {
	    fprintf(stderr,"Unrecognized distance type: %s\n",Dtype);
	    return(-1.0);
	}
}

int get_info(char *fname,int lim)
{
        int i,j;
        char *input;
        FILE *fp;

        input=kalloc(lim);

	Stype=kalloc(200);
	Xtype=kalloc(200);
	Mtype=kalloc(200);
	Ftype=kalloc(200);
	Dtype=kalloc(200);

        if((fp=fopen(fname,"r"))==NULL){
            fprintf(stderr,"Can't find file %s\n",fname);
            return(1);
        }

        sscanf(parse_line(fp,lim),"%d",&cflag);
        sscanf(parse_line(fp,lim),"%ld",&seed);
        sscanf(parse_line(fp,lim),"%d",&ndeme);
        sscanf(parse_line(fp,lim),"%d",&nmod);
        sscanf(parse_line(fp,lim),"%d",&nelite);
        sscanf(parse_line(fp,lim),"%d",&ngen);
        sscanf(parse_line(fp,lim),"%f",&pc);
        sscanf(parse_line(fp,lim),"%f",&pm);
        sscanf(parse_line(fp,lim),"%f",&Rc);
        sscanf(parse_line(fp,lim),"%s",Dtype);
        sscanf(parse_line(fp,lim),"%s",Xtype);
        sscanf(parse_line(fp,lim),"%s",Mtype);
        sscanf(parse_line(fp,lim),"%s",Ftype);
        sscanf(parse_line(fp,lim),"%s",Stype);

/*Initialize RNG*/
        seed=seed_random(seed);

        i=0;
        while((input=parse_line(fp,lim))!=NULL){
            sscanf(input,"%f%f%f%f",&bot[i],&top[i],&dev[i],&d_weight[i]);
	    dev[i]*=(top[i]-bot[i]);
            i++;
        }
        npar=i;

/*Initialize global elistist population*/
	if((Epop.d=(deme*)malloc(ndeme*sizeof(deme)))==NULL){
            fprintf(stderr,"Can't allocate memory\n");
            return(1);
        }
        for(i=0;i<ndeme;i++){
            if((Epop.d[i].m=(model*)malloc(nelite*sizeof(model)))==NULL){
                fprintf(stderr,"Can't allocate memory\n");
                return(1);
            }
            for(j=0;j<nelite;j++){
                if((Epop.d[i].m[j].par=(float*)malloc(npar*sizeof(float)))==NULL){
                    fprintf(stderr,"Can't allocate memory\n");
                    return(1);
                }
	        Epop.d[i].m[j].cost=-1.0;
            }
        }

        fclose(fp);
	return(0);
}

char* parse_line(FILE *fp,int len)
{
        char *string,*line;

        line=(char*)malloc(200*sizeof(char));
        string=(char*)malloc(len*sizeof(char));

        if(fgets(line,200*sizeof(char),fp)==NULL){
            return(NULL);
        } else {
            strncpy(string,line,len);
            return(string);
        }
}

float* falloc(int x)
{
        float *array;

        if((array=(float*)malloc(x*sizeof(float)))==NULL){
            fprintf(stderr,"Memory Allocation Error\n");
            exit(0);
        }

        return(array);
}

char* kalloc(int x)
{
        char *array;

        if((array=(char*)malloc(x*sizeof(char)))==NULL){
            fprintf(stderr,"Memory Allocation Error\n");
            exit(0);
        }

        return(array);
}
