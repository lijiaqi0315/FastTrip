FC 		= gfortran
FFLAG 	= -O3
#OBJDIR	= code09-t
OBJECT 	= bessj0.o bessj1.o bessj.o caxcb.o cdgemp.o cmemcpy.o\
		  four1.o getdata.o qsam2ve.o qsbsj.o qsfftinv.o\
		  qsgetinp.o qshksh.o qskern.o qslayer.o qsmain.o\
		  qsmultis.o qspsv.o qsqmodel.o qssh.o qssource.o\
		  qssublay.o qsve2am.o qswavelet.o qswaveno.o qswvint.o\
		  taper.o
EXE		= qseis

$(EXE) : $(OBJECT)
	$(FC) $(OBJECT) -o $(EXE)

$(OBJECT): %.o : %.f
	$(FC) $(FFLAG) -c $< -o $@

#.PHONY: clean

clean:
	rm -f *.o $(EXE)
