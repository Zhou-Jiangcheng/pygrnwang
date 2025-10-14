FORFILES =  edcdisc.f    edcgrn.f     edcmain.f    edcoutput.f  getdata.f  \
            okada.f      dc3d.f       dc3d0.f

OBJECTS =	$(FORFILES:%.f=%.o)

#FFLAGS =	-O3

PROGRAM =	edcmp2.0

$(PROGRAM): 	$(OBJECTS)
		$(FC) $(FFLAGS) $(OBJECTS) -o $@

edcdisc.o:edcglobal.h edcdisc.f
edcgrn.o:edcglobal.h edcgrn.f
edcmain.o:edcglobal.h edcmain.f
edcoutput.o:edcglobal.h edcoutput.f
getdata.o:getdata.f
okada.o:okada.f
dc3d.o:dc3d.f
dc3d0.o:dc3d0.f






