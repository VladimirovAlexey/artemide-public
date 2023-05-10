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
#Fflags=  
#path to fortran compilator (needed for f2py)
Fpath=/usr/bin/f95

#options for COMILATOR to compile QCDinput. e.g. link to LHA
FOPT=$(shell lhapdf-config --ldflags)
#### for debuging -g -fbacktrace -ffpe-trap=zero,overflow,underflow
#FOPT=-L/home/vla18041/LinkData2/LHAPDF/Installation/lib -lLHAPDF -lstdc++
#FOPT=-L/home/alexey/WorkingFiles/LHAPDF/Intallation/lib -lLHAPDF -lstdc++

################################################################### LIST OF FILES ####################################
SOURCEDIR       = $(aTMDeHOME)/src
BIN		= $(aTMDeHOME)/bin
OBJ		= $(aTMDeHOME)/obj
MOD		= $(aTMDeHOME)/mod
HDIR		= $(aTMDeHOME)/harpy

aTMDeFILES = \
$(SOURCEDIR)/Code/aTMDe_Numerics.f90 \
$(SOURCEDIR)/Code/IO_functions.f90 \
$(SOURCEDIR)/Code/IntegrationRoutines.f90 \
$(SOURCEDIR)/LeptonCutsDY.f90 \
$(SOURCEDIR)/aTMDe_setup.f90 \
$(SOURCEDIR)/QCDinput.f90 \
$(SOURCEDIR)/EWinput.f90 \
$(SOURCEDIR)/TMD_AD.f90 \
$(SOURCEDIR)/Model/TMDR_model.f90 \
$(SOURCEDIR)/TMDR.f90 \
$(SOURCEDIR)/Model/uTMDPDF_model.f90 \
$(SOURCEDIR)/uTMDPDF.f90 \
$(SOURCEDIR)/Model/uTMDFF_model.f90 \
$(SOURCEDIR)/uTMDFF.f90 \
$(SOURCEDIR)/Model/lpTMDPDF_model.f90 \
$(SOURCEDIR)/lpTMDPDF.f90 \
$(SOURCEDIR)/Model/SiversTMDPDF_model.f90 \
$(SOURCEDIR)/SiversTMDPDF.f90 \
$(SOURCEDIR)/Model/wgtTMDPDF_model.f90 \
$(SOURCEDIR)/wgtTMDPDF.f90 \
$(SOURCEDIR)/TMDs.f90 \
$(SOURCEDIR)/TMDF.f90 \
$(SOURCEDIR)/TMDs_inKT.f90 \
$(SOURCEDIR)/TMDX_DY.f90 \
$(SOURCEDIR)/TMDX_SIDIS.f90 \
$(SOURCEDIR)/aTMDe_control.f90 

Twist2Files=\
$(SOURCEDIR)/Code/Twist2/Twist2Convolution.f90 \
$(SOURCEDIR)/Code/Grids/TMDGrid-B.f90 \
$(SOURCEDIR)/Code/Twist2/Twist2Convolution-VAR.f90 \
$(SOURCEDIR)/Code/Grids/TMDGrid-B-VAR.f90

Twist3Files=\
$(SOURCEDIR)/Code/Grids/TMDGrid-B-2.f90 \
$(SOURCEDIR)/Code/Grids/TMDGrid-B-VAR.f90 \
$(SOURCEDIR)/Code/Twist3/Twist3Convolution.f90 \
$(SOURCEDIR)/Code/Twist3/Twist3Convolution-VAR.f90 

TMD_ADFiles=\
$(SOURCEDIR)/Code/TMD_AD/AD_primary.f90 \
$(SOURCEDIR)/Code/TMD_AD/AD_secondary.f90 \
$(SOURCEDIR)/Code/TMD_AD/AD_atMu.f90 \
$(SOURCEDIR)/Code/TMD_AD/exactZetaLine.f90 \
$(SOURCEDIR)/Code/TMD_AD/AD_Integral.f90 

TMDRFiles=\
$(SOURCEDIR)/Code/TMDR/type1.f90 \
$(SOURCEDIR)/Code/TMDR/type2.f90 \
$(SOURCEDIR)/Code/TMDR/type3.f90

uTMDPDFFiles=\
$(SOURCEDIR)/Code/uTMDPDF/coeffFunc.f90 \
$(SOURCEDIR)/Code/uTMDPDF/convolutions.f90 \
$(SOURCEDIR)/Code/uTMDPDF/modelTest.f90

uTMDFFFiles=\
$(SOURCEDIR)/Code/uTMDFF/coeffFunc.f90 \
$(SOURCEDIR)/Code/uTMDFF/convolutions.f90 \
$(SOURCEDIR)/Code/uTMDFF/modelTest.f90

lpTMDPDFFiles=\
$(SOURCEDIR)/Code/lpTMDPDF/coeffFunc.f90 \
$(SOURCEDIR)/Code/lpTMDPDF/convolutions.f90 \
$(SOURCEDIR)/Code/lpTMDPDF/modelTest.f90

SiversTMDPDFFiles=\
$(SOURCEDIR)/Code/SiversTMDPDF/modelTest.f90 \
$(SOURCEDIR)/Code/SiversTMDPDF/convolutions.f90

wgtTMDPDFFiles=\
$(SOURCEDIR)/Code/wgtTMDPDF/coeffFunc.f90 \
$(SOURCEDIR)/Code/wgtTMDPDF/convolutions.f90 \
$(SOURCEDIR)/Code/wgtTMDPDF/modelTest.f90

TMDsFiles=\
$(SOURCEDIR)/Code/TMDs/TMD-calls.f90 

aTMDeSetupFiles=\
$(SOURCEDIR)/Code/aTMDe_setup/const-modification.f90 

ExtraFiles=\
$(SOURCEDIR)/DYcoeff-func.f90

aTMDeMODEL = \
$(SOURCEDIR)/Model/TMDR_model.f90 \
$(SOURCEDIR)/Model/TMDs_model.f90 \
$(SOURCEDIR)/Model/uTMDFF_model.f90 \
$(SOURCEDIR)/Model/uTMDPDF_model.f90 

aTMDeOBJ = \
$(OBJ)/aTMDe_Numerics.o \
$(OBJ)/IO_functions.o \
$(OBJ)/IntegrationRoutines.o \
$(OBJ)/LeptonCutsDY.o \
$(OBJ)/aTMDe_setup.o \
$(OBJ)/QCDinput.o \
$(OBJ)/EWinput.o\
$(OBJ)/TMD_AD.o\
$(OBJ)/TMDR_model.o\
$(OBJ)/TMDR.o\
$(OBJ)/uTMDPDF_model.o \
$(OBJ)/uTMDPDF.o \
$(OBJ)/uTMDFF_model.o \
$(OBJ)/uTMDFF.o\
$(OBJ)/lpTMDPDF_model.o \
$(OBJ)/lpTMDPDF.o \
$(OBJ)/SiversTMDPDF_model.o \
$(OBJ)/SiversTMDPDF.o \
$(OBJ)/wgtTMDPDF_model.o \
$(OBJ)/wgtTMDPDF.o \
$(OBJ)/TMDs.o \
$(OBJ)/TMDF.o \
$(OBJ)/TMDs_inKT.o \
$(OBJ)/TMDX_DY.o \
$(OBJ)/TMDX_SIDIS.o \
$(OBJ)/aTMDe_control.o 

#these are utility object needed to compale any artemide module
aTMDeUTILITY = \
$(OBJ)/aTMDe_Numerics.o \
$(OBJ)/IO_functions.o \
$(OBJ)/IntegrationRoutines.o


################################################################### COMPILATION OF ARTEMIDE ####################################
FC=$(FCompilator) $(Fflags)

.PHONY: clean default obj program test harpy harpy-signature

default: obj

update: $(BIN)/update-const
	./bin/update-const $(TARGET)


obj: $(aTMDeOBJ) $(aTMDeFILES) $(aTMDeMODEL) $(Twist2Files) $(TMD_ADFiles) $(TMDRFiles) $(uTMDPDFFiles) $(uTMDFFFiles) $(lpTMDFFFiles)

$(OBJ)/aTMDe_Numerics.o: $(SOURCEDIR)/Code/aTMDe_Numerics.f90
	$(FC) -c $(SOURCEDIR)/Code/aTMDe_Numerics.f90
	mv *.o $(OBJ)
	mv *.mod $(MOD)

$(OBJ)/IO_functions.o: $(SOURCEDIR)/Code/IO_functions.f90 $(OBJ)/aTMDe_Numerics.o
	$(FC) -c $(SOURCEDIR)/Code/IO_functions.f90 -I$(MOD)
	mv *.o $(OBJ)
	mv *.mod $(MOD)

$(OBJ)/IntegrationRoutines.o: $(SOURCEDIR)/Code/IntegrationRoutines.f90 $(OBJ)/aTMDe_Numerics.o $(OBJ)/IO_functions.o
	$(FC) -c $(SOURCEDIR)/Code/IntegrationRoutines.f90 -I$(MOD)
	mv *.o $(OBJ)
	mv *.mod $(MOD)

$(OBJ)/LeptonCutsDY.o: $(SOURCEDIR)/LeptonCutsDY.f90 $(aTMDeUTILITY)
#	mkdir -p obj
#	mkdir -p mod
	$(FC) -c $(SOURCEDIR)/LeptonCutsDY.f90 -I$(MOD)
	mv *.o $(OBJ)
	mv *.mod $(MOD)
	
$(OBJ)/QCDinput.o: $(SOURCEDIR)/QCDinput.f90 $(aTMDeUTILITY)
#	mkdir -p obj
	$(FC) -c $(SOURCEDIR)/QCDinput.f90 -I$(MOD) $(FOPT)
	mv *.o $(OBJ)
	mv *.mod $(MOD)
	
$(OBJ)/EWinput.o: $(SOURCEDIR)/EWinput.f90 $(aTMDeUTILITY)
#	mkdir -p obj
	$(FC) -c $(SOURCEDIR)/EWinput.f90 -I$(MOD)
	mv *.o $(OBJ)
	mv *.mod $(MOD)
	
$(OBJ)/uTMDPDF_model.o: $(SOURCEDIR)/Model/uTMDPDF_model.f90 $(aTMDeUTILITY)
#	mkdir -p obj
	$(FC) -c $(SOURCEDIR)/Model/uTMDPDF_model.f90 -I$(MOD)
	mv *.o $(OBJ)
	mv *.mod $(MOD)

$(OBJ)/uTMDPDF.o: $(SOURCEDIR)/uTMDPDF.f90 $(OBJ)/QCDinput.o $(SOURCEDIR)/Model/uTMDPDF_model.f90 $(Twist2Files) $(aTMDeUTILITY) $(uTMDPDFFiles)
#	mkdir -p obj
	$(FC) -c $(SOURCEDIR)/uTMDPDF.f90 -I$(MOD)
	mv *.o $(OBJ)
	mv *.mod $(MOD)

$(OBJ)/uTMDFF_model.o: $(SOURCEDIR)/Model/uTMDFF_model.f90 $(aTMDeUTILITY)
#	mkdir -p obj
	$(FC) -c $(SOURCEDIR)/Model/uTMDFF_model.f90 -I$(MOD)
	mv *.o $(OBJ)
	mv *.mod $(MOD)
	
$(OBJ)/uTMDFF.o: $(SOURCEDIR)/uTMDFF.f90 $(OBJ)/QCDinput.o $(SOURCEDIR)/Model/uTMDFF_model.f90 $(Twist2Files) $(aTMDeUTILITY) $(uTMDFFFiles)
#	mkdir -p obj
	$(FC) -c $(SOURCEDIR)/uTMDFF.f90 -I$(MOD)
	mv *.o $(OBJ)
	mv *.mod $(MOD)
	
$(OBJ)/lpTMDPDF_model.o: $(SOURCEDIR)/Model/lpTMDPDF_model.f90 $(aTMDeUTILITY)
#	mkdir -p obj
	$(FC) -c $(SOURCEDIR)/Model/lpTMDPDF_model.f90 -I$(MOD)
	mv *.o $(OBJ)
	mv *.mod $(MOD)
	
$(OBJ)/lpTMDPDF.o: $(SOURCEDIR)/lpTMDPDF.f90 $(OBJ)/QCDinput.o $(SOURCEDIR)/Model/lpTMDPDF_model.f90 $(Twist2Files) $(aTMDeUTILITY) $(lpTMDPDFFiles)
#	mkdir -p obj
	$(FC) -c $(SOURCEDIR)/lpTMDPDF.f90 -I$(MOD)
	mv *.o $(OBJ)
	mv *.mod $(MOD)

$(OBJ)/SiversTMDPDF_model.o: $(SOURCEDIR)/Model/SiversTMDPDF_model.f90 $(aTMDeUTILITY)
#	mkdir -p obj
	$(FC) -c $(SOURCEDIR)/Model/SiversTMDPDF_model.f90 -I$(MOD)
	mv *.o $(OBJ)
	mv *.mod $(MOD)
	
$(OBJ)/SiversTMDPDF.o: $(SOURCEDIR)/SiversTMDPDF.f90 $(OBJ)/QCDinput.o $(SOURCEDIR)/Model/SiversTMDPDF_model.f90 $(Twist3Files) $(aTMDeUTILITY) $(SiversTMDPDFFiles)
#	mkdir -p obj
	$(FC) -c $(SOURCEDIR)/SiversTMDPDF.f90 -I$(MOD)
	mv *.o $(OBJ)
	mv *.mod $(MOD)

$(OBJ)/wgtTMDPDF_model.o: $(SOURCEDIR)/Model/wgtTMDPDF_model.f90 $(aTMDeUTILITY)
#	mkdir -p obj
	$(FC) -c $(SOURCEDIR)/Model/wgtTMDPDF_model.f90 -I$(MOD)
	mv *.o $(OBJ)
	mv *.mod $(MOD)
	
$(OBJ)/wgtTMDPDF.o: $(SOURCEDIR)/wgtTMDPDF.f90 $(OBJ)/QCDinput.o $(SOURCEDIR)/Model/wgtTMDPDF_model.f90 $(Twist2Files) $(aTMDeUTILITY) $(wgtTMDPDFFiles)
#	mkdir -p obj
	$(FC) -c $(SOURCEDIR)/wgtTMDPDF.f90 -I$(MOD)
	mv *.o $(OBJ)
	mv *.mod $(MOD)
	
$(OBJ)/TMD_AD.o: $(SOURCEDIR)/TMD_AD.f90 $(aTMDeUTILITY) $(TMD_ADFiles)
#	mkdir -p obj
	$(FC) -c $(SOURCEDIR)/TMD_AD.f90 -I$(MOD)
	mv *.o $(OBJ)
	mv *.mod $(MOD)
	
$(OBJ)/TMDR_model.o: $(SOURCEDIR)/Model/TMDR_model.f90 $(aTMDeUTILITY) $(TMD_ADFiles) $(OBJ)/TMD_AD.o
#	mkdir -p obj
	$(FC) -c $(SOURCEDIR)/Model/TMDR_model.f90 -I$(MOD)
	mv *.o $(OBJ)
	mv *.mod $(MOD)
	
$(OBJ)/TMDR.o: $(SOURCEDIR)/TMDR.f90 $(SOURCEDIR)/Model/TMDR_model.f90 $(OBJ)/QCDinput.o $(OBJ)/TMD_AD.o $(aTMDeUTILITY) $(TMDRFiles)
#	mkdir -p obj
	$(FC) -c $(SOURCEDIR)/TMDR.f90 -I$(MOD)
	mv *.o $(OBJ)
	mv *.mod $(MOD)
	
$(OBJ)/TMDs.o: $(SOURCEDIR)/TMDs.f90 $(OBJ)/uTMDPDF.o $(OBJ)/uTMDFF.o $(OBJ)/lpTMDPDF.o $(OBJ)/SiversTMDPDF.o $(OBJ)/wgtTMDPDF.o $(SOURCEDIR)/Model/TMDs_model.f90 $(OBJ)/TMDR.o $(TMDsFiles) $(aTMDeUTILITY)
#	mkdir -p obj
	$(FC) -c $(SOURCEDIR)/TMDs.f90 -I$(MOD)
	mv *.o $(OBJ)
	mv *.mod $(MOD)

$(OBJ)/TMDs_inKT.o: $(SOURCEDIR)/TMDs_inKT.f90 $(OBJ)/TMDs.o $(aTMDeUTILITY)
#	mkdir -p obj
	$(FC) -c $(SOURCEDIR)/TMDs_inKT.f90 -I$(MOD)
	mv *.o $(OBJ)
	mv *.mod $(MOD)
	
$(OBJ)/TMDF.o: $(SOURCEDIR)/TMDF.f90 $(OBJ)/TMDs.o $(aTMDeUTILITY)
#	mkdir -p obj
	$(FC) -c $(SOURCEDIR)/TMDF.f90 -I$(MOD)
	mv *.o $(OBJ)
	mv *.mod $(MOD)

$(OBJ)/TMDX_DY.o: $(SOURCEDIR)/TMDX_DY.f90 $(SOURCEDIR)/DYcoeff-func.f90 $(OBJ)/TMDF.o  $(OBJ)/QCDinput.o $(aTMDeUTILITY)
#	mkdir -p obj
	$(FC) -c $(SOURCEDIR)/TMDX_DY.f90 -I$(MOD)
	mv *.o $(OBJ)
	mv *.mod $(MOD)
	
$(OBJ)/TMDX_SIDIS.o: $(SOURCEDIR)/TMDX_SIDIS.f90 $(OBJ)/TMDs.o $(OBJ)/QCDinput.o $(aTMDeUTILITY)
#	mkdir -p obj
	$(FC) -c $(SOURCEDIR)/TMDX_SIDIS.f90 -I$(MOD)
	mv *.o $(OBJ)
	mv *.mod $(MOD)

$(OBJ)/aTMDe_setup.o: $(SOURCEDIR)/aTMDe_setup.f90 $(aTMDeUTILITY) $(aTMDeSetupFiles)
	$(FC) -c $(SOURCEDIR)/aTMDe_setup.f90 -I$(MOD)
	mv *.o $(OBJ)
	mv *.mod $(MOD)
	
$(OBJ)/aTMDe_control.o: $(SOURCEDIR)/aTMDe_control.f90 $(OBJ)/aTMDe_setup.o $(aTMDeUTILITY)
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
	$(RM) $(HDIR)/*.so
	$(RM) $(HDIR)/__pycache__/*.*
	
program: 
	echo $(TARGET)
	#$(FC) $(aTMDeHOME)/Prog/$(TARGET) $(aTMDeOBJ) $(FOPT) -I$(MOD)
	$(FC) $(TARGET) $(aTMDeOBJ) $(FOPT) -I$(MOD)
	
test: 
	$(FC) $(aTMDeHOME)/Prog/test.f90 $(aTMDeOBJ) $(FOPT) -I$(MOD)
	./a.out
	
################################################ update constants part ##############################

$(BIN)/update-const: $(aTMDeHOME)/Prog/update-constants-file.f90 $(OBJ)/aTMDe_setup.o $(aTMDeUTILITY)
	$(FC) $(aTMDeHOME)/Prog/update-constants-file.f90 $(aTMDeOBJ) $(FOPT) -I$(MOD) -o update-const
	mv update-const $(BIN)/update-const
	
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
	mv artemide*.so $(HDIR)
