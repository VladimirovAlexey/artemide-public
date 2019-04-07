###########################################################################################################
#			Make file for artemide + harpy
#				modify the first section with your values
###########################################################################################################
# location of artemide
aTMDeHOME       = $(PWD)

#PUT YOUR FORTRAN COMPILER
FCompilator=f95 
#PUT HERE extra flags for compilator (put "space" if not flags requared)
Fflags= -fopenmp
#path to fortran compilator (needed for f2py)
Fpath=/usr/bin/f95

#options for COMILATOR to compile QCDinput. e.g. link to LHA
FOPT=$(shell lhapdf-config --ldflags)

################################################################### LIST OF FILES ####################################
SOURCEDIR       = $(aTMDeHOME)/src
BIN		= $(aTMDeHOME)/bin
OBJ		= $(aTMDeHOME)/obj
MOD		= $(aTMDeHOME)/mod
HDIR		= $(aTMDeHOME)/harpy

aTMDeFILES = \
$(SOURCEDIR)/LeptonCutsDY.f90 \
$(SOURCEDIR)/aTMDe_setup.f90 \
$(SOURCEDIR)/QCDinput.f90 \
$(SOURCEDIR)/EWinput.f90 \
$(SOURCEDIR)/TMDR.f90 \
$(SOURCEDIR)/uTMDPDF-MODELinterface.f90 \
$(SOURCEDIR)/uTMDPDF.f90 \
$(SOURCEDIR)/uTMDFF-MODELinterface.f90 \
$(SOURCEDIR)/uTMDFF.f90 \
$(SOURCEDIR)/TMDs.f90 \
$(SOURCEDIR)/TMDF.f90 \
$(SOURCEDIR)/TMDs_inKT.f90 \
$(SOURCEDIR)/TMDX_DY.f90 \
$(SOURCEDIR)/TMDX_SIDIS.f90 \
$(SOURCEDIR)/aTMDe_control.f90 



CommonFiles=\
$(SOURCEDIR)/CommonCode/Twist2Convolution.f90 \
$(SOURCEDIR)/CommonCode/Twist2Grid.f90 \
$(SOURCEDIR)/CommonCode/Twist2Convolution-VAR.f90 \
$(SOURCEDIR)/CommonCode/Twist2Grid-VAR.f90 

aTMDeMODEL = \
$(SOURCEDIR)/Model/TMDR_model.f90 \
$(SOURCEDIR)/Model/TMDs_model.f90 \
$(SOURCEDIR)/Model/uTMDFF_model.f90 \
$(SOURCEDIR)/Model/uTMDPDF_model.f90 

aTMDeOBJ = \
$(OBJ)/LeptonCutsDY.o \
$(OBJ)/aTMDe_setup.o \
$(OBJ)/QCDinput.o \
$(OBJ)/EWinput.o\
$(OBJ)/TMDR.o\
$(OBJ)/uTMDPDF-MODELinterface.o \
$(OBJ)/uTMDPDF.o \
$(OBJ)/uTMDFF-MODELinterface.o \
$(OBJ)/uTMDFF.o\
$(OBJ)/TMDs.o \
$(OBJ)/TMDF.o \
$(OBJ)/TMDs_inKT.o \
$(OBJ)/TMDX_DY.o \
$(OBJ)/TMDX_SIDIS.o \
$(OBJ)/aTMDe_control.o 


################################################################### COMPILATION OF ARTEMIDE ####################################
FC=$(FCompilator) $(Fflags)

.PHONY: clean default obj program test harpy harpy-signature

default: obj


obj: $(aTMDeOBJ) $(aTMDeFILES) $(aTMDeMODEL) $(CommonFiles)

$(OBJ)/uTMDPDF-MODELinterface.o: $(SOURCEDIR)/uTMDPDF-MODELinterface.f90 $(SOURCEDIR)/Model/uTMDPDF_model.f90
	$(FC) -c $(SOURCEDIR)/uTMDPDF-MODELinterface.f90
	mv *.o $(OBJ)
	mv *.mod $(MOD)
	
$(OBJ)/uTMDFF-MODELinterface.o: $(SOURCEDIR)/uTMDFF-MODELinterface.f90 $(SOURCEDIR)/Model/uTMDFF_model.f90
	$(FC) -c $(SOURCEDIR)/uTMDFF-MODELinterface.f90
	mv *.o $(OBJ)
	mv *.mod $(MOD)
	
	
$(OBJ)/LeptonCutsDY.o: $(SOURCEDIR)/LeptonCutsDY.f90
	mkdir -p obj
	mkdir -p mod
	$(FC) -c $(SOURCEDIR)/LeptonCutsDY.f90
	mv *.o $(OBJ)
	mv *.mod $(MOD)
	
$(OBJ)/fDSS.o: $(SOURCEDIR)/fDSS.f
#	mkdir -p obj
	$(FC) -c $(SOURCEDIR)/fDSS.f
	mv *.o $(OBJ)
	#mv *.mod $(MOD)
	
$(OBJ)/QCDinput.o: $(SOURCEDIR)/QCDinput.f90
#	mkdir -p obj
	$(FC) -c $(SOURCEDIR)/QCDinput.f90 $(FOPT)
	mv *.o $(OBJ)
	mv *.mod $(MOD)
	
$(OBJ)/EWinput.o: $(SOURCEDIR)/EWinput.f90
#	mkdir -p obj
	$(FC) -c $(SOURCEDIR)/EWinput.f90
	mv *.o $(OBJ)
	mv *.mod $(MOD)

$(OBJ)/uTMDPDF.o: $(SOURCEDIR)/uTMDPDF.f90 $(OBJ)/QCDinput.o $(SOURCEDIR)/Model/uTMDPDF_model.f90 $(CommonFiles)
#	mkdir -p obj
	$(FC) -c $(SOURCEDIR)/uTMDPDF.f90 -I$(MOD)
	mv *.o $(OBJ)
	mv *.mod $(MOD)
	
$(OBJ)/uTMDFF.o: $(SOURCEDIR)/uTMDFF.f90 $(OBJ)/QCDinput.o $(SOURCEDIR)/Model/uTMDFF_model.f90 $(CommonFiles) $(OBJ)/fDSS.o
#	mkdir -p obj
	$(FC) -c $(SOURCEDIR)/uTMDFF.f90 -I$(MOD)
	mv *.o $(OBJ)
	mv *.mod $(MOD)

$(OBJ)/TMDR.o: $(SOURCEDIR)/TMDR.f90 $(SOURCEDIR)/Model/TMDR_model.f90 $(OBJ)/QCDinput.o
#	mkdir -p obj
	$(FC) -c $(SOURCEDIR)/TMDR.f90 -I$(MOD)
	mv *.o $(OBJ)
	mv *.mod $(MOD)
	
$(OBJ)/TMDs.o: $(SOURCEDIR)/TMDs.f90 $(OBJ)/uTMDPDF.o $(OBJ)/uTMDFF.o $(SOURCEDIR)/Model/TMDs_model.f90 $(OBJ)/TMDR.o
#	mkdir -p obj
	$(FC) -c $(SOURCEDIR)/TMDs.f90 -I$(MOD)
	mv *.o $(OBJ)
	mv *.mod $(MOD)

$(OBJ)/TMDs_inKT.o: $(SOURCEDIR)/TMDs_inKT.f90 $(OBJ)/TMDs.o
#	mkdir -p obj
	$(FC) -c $(SOURCEDIR)/TMDs_inKT.f90 -I$(MOD)
	mv *.o $(OBJ)
	mv *.mod $(MOD)
	
$(OBJ)/TMDF.o: $(SOURCEDIR)/TMDF.f90 $(OBJ)/TMDs.o
#	mkdir -p obj
	$(FC) -c $(SOURCEDIR)/TMDF.f90 -I$(MOD)
	mv *.o $(OBJ)
	mv *.mod $(MOD)

$(OBJ)/TMDX_DY.o: $(SOURCEDIR)/TMDX_DY.f90 $(OBJ)/TMDF.o  $(OBJ)/QCDinput.o
#	mkdir -p obj
	$(FC) -c $(SOURCEDIR)/TMDX_DY.f90 -I$(MOD)
	mv *.o $(OBJ)
	mv *.mod $(MOD)
	
$(OBJ)/TMDX_SIDIS.o: $(SOURCEDIR)/TMDX_SIDIS.f90 $(OBJ)/TMDs.o $(OBJ)/QCDinput.o
#	mkdir -p obj
	$(FC) -c $(SOURCEDIR)/TMDX_SIDIS.f90 -I$(MOD)
	mv *.o $(OBJ)
	mv *.mod $(MOD)

$(OBJ)/aTMDe_setup.o: $(SOURCEDIR)/aTMDe_setup.f90
	$(FC) -c $(SOURCEDIR)/aTMDe_setup.f90 -I$(MOD)
	mv *.o $(OBJ)
	mv *.mod $(MOD)
	
$(OBJ)/aTMDe_control.o: $(SOURCEDIR)/aTMDe_control.f90 
	$(FC) -c $(SOURCEDIR)/aTMDe_control.f90 -I$(MOD)
	mv *.o $(OBJ)
	mv *.mod $(MOD)

clean: 
	$(RM) a.out
	$(RM) aTMDe-temporary
	$(RM) count *.o *.mod
	$(RM) count $(OBJ)/*.o
	$(RM) count $(MOD)/*.mod
	$(RM) $(HDIR)/*.pyc
	$(RM) $(HDIR)/artemide.so
	
program: 
	echo $(TARGET)
	$(FC) $(aTMDeHOME)/Prog/$(TARGET) $(aTMDeOBJ) $(FOPT) -I$(MOD)
	
test: 
	$(FC) $(aTMDeHOME)/Prog/test.f90 $(aTMDeOBJ) $(FOPT) -I$(MOD)
	./a.out
	
################################################  HARPY PART  #######################################

harpy-signature: 
	f2py -h $(HDIR)/artemide.pyf --overwrite-signature $(HDIR)/harpy.f90
	sed -i '3i\\' $(HDIR)/artemide.pyf	
	sed -i '3i interface' $(HDIR)/artemide.pyf
	sed -i '3i python module artemide' $(HDIR)/artemide.pyf
	sed -i '3i\\' $(HDIR)/artemide.pyf
	echo 'end interface' >> $(HDIR)/artemide.pyf
	echo 'end python module artemide' >> $(HDIR)/artemide.pyf

harpy: 
	f2py -c --f90exec=$(Fpath) --f90flags=$(Fflags) $(FOPT) -lgomp -I$(MOD) $(aTMDeFILES) $(HDIR)/harpy.f90 $(HDIR)/artemide.pyf
	mv artemide.so $(HDIR)
