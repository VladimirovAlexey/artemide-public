###########################################################################################################
#			Make file for artemide + harpy
#				modify the first section with your values
###########################################################################################################
# location of artemide
aTMDeHOME       = $(PWD)

#PUT YOUR FORTRAN COMPILER
FCompilator=f95 

#PUT HERE extra flags for compilator (put "space" if not flags requared)
Fflags= -O3 -cpp -march=native  -fforce-addr -fstrength-reduce -fcaller-saves -funroll-loops -fopenmp
# This should be as string because it is how f2py eats it.
FflagsHARPY= '-O3 -cpp -march=native  -fforce-addr -fstrength-reduce -fcaller-saves -funroll-loops -fopenmp'
#### Fir debuging
#Fflags=  -O3 -cpp -march=native  -fforce-addr -fstrength-reduce -fcaller-saves -funroll-loops -Wall -fopenmp
#path to fortran compilator (needed for f2py)
Fpath=/usr/bin/f95
F77path=/usr/bin/f77

#options for COMILATOR to compile QCDinput. e.g. link to LHA
#FOPT=$(shell lhapdf-config --ldflags)
FOPT=

#### for debuging -g -fbacktrace -fcheck=all
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
$(SOURCEDIR)/Code/InverseMatrix.f90 \
$(SOURCEDIR)/Code/LHA/LHA_alpha.f90 \
$(SOURCEDIR)/LeptonCutsDY.f90 \
$(SOURCEDIR)/aTMDe_setup.f90 \
$(SOURCEDIR)/QCDinput.f90 \
$(SOURCEDIR)/EWinput.f90 \
$(SOURCEDIR)/TMD_AD.f90 \
$(SOURCEDIR)/Model/TMDR_model.f90 \
$(SOURCEDIR)/TMDR.f90 \
$(SOURCEDIR)/Model/uTMDPDF_model.f90 \
$(SOURCEDIR)/uTMDPDF_OPE.f90 \
$(SOURCEDIR)/uTMDPDF.f90 \
$(SOURCEDIR)/Model/uTMDFF_model.f90 \
$(SOURCEDIR)/uTMDFF_OPE.f90 \
$(SOURCEDIR)/uTMDFF.f90 \
$(SOURCEDIR)/Model/lpTMDPDF_model.f90 \
$(SOURCEDIR)/lpTMDPDF_OPE.f90 \
$(SOURCEDIR)/lpTMDPDF.f90 \
$(SOURCEDIR)/Model/SiversTMDPDF_model.f90 \
$(SOURCEDIR)/SiversTMDPDF_OPE.f90 \
$(SOURCEDIR)/SiversTMDPDF.f90 \
$(SOURCEDIR)/Model/wgtTMDPDF_model.f90 \
$(SOURCEDIR)/wgtTMDPDF_OPE.f90 \
$(SOURCEDIR)/wgtTMDPDF.f90 \
$(SOURCEDIR)/Model/wglTMDPDF_model.f90 \
$(SOURCEDIR)/wglTMDPDF_OPE.f90 \
$(SOURCEDIR)/wglTMDPDF.f90 \
$(SOURCEDIR)/Model/BoerMuldersTMDPDF_model.f90 \
$(SOURCEDIR)/BoerMuldersTMDPDF_OPE.f90 \
$(SOURCEDIR)/BoerMuldersTMDPDF.f90 \
$(SOURCEDIR)/eeTMDFF.f90 \
$(SOURCEDIR)/Model/eeTMDFF_model.f90 \
$(SOURCEDIR)/TMDF.f90 \
$(SOURCEDIR)/TMDF_KPC.f90 \
$(SOURCEDIR)/TMDX_DY.f90 \
$(SOURCEDIR)/TMDX_SIDIS.f90 \
$(SOURCEDIR)/aTMDe_control.f90

Twist2Files=\
$(SOURCEDIR)/Code/Twist2/Twist2Convolution.f90 \
$(SOURCEDIR)/Code/Twist2/largeX_ADs.f90 \
$(SOURCEDIR)/Code/Twist2/Twist2_WW.f90 \
$(SOURCEDIR)/Code/Twist2/Twist2LargeX.f90 \
$(SOURCEDIR)/Code/Twist2/Twist2_ChGrid.f90 \
$(SOURCEDIR)/Code/Twist2/Twist2-AS-term.f90

Twist3Files=\
$(SOURCEDIR)/Code/Twist3/placeHolder.f90

KTspaceFiles=\
$(SOURCEDIR)/Code/KTspace/Fourier_Levin.f90\
$(SOURCEDIR)/Code/KTspace/Fourier.f90\
$(SOURCEDIR)/Code/KTspace/Moment.f90\
$(SOURCEDIR)/Code/KTspace/grid_inKT.f90

TMD_ADFiles=\
$(SOURCEDIR)/Code/TMD_AD/AD_primary.f90 \
$(SOURCEDIR)/Code/TMD_AD/AD_secondary.f90 \
$(SOURCEDIR)/Code/TMD_AD/AD_atMu.f90 \
$(SOURCEDIR)/Code/TMD_AD/exactZetaLine.f90 \
$(SOURCEDIR)/Code/TMD_AD/AD_Integral.f90 

TMDRFiles=\
$(SOURCEDIR)/Code/TMDR/placeHolder.f90

uTMDPDFFiles=\
$(SOURCEDIR)/Code/uTMDPDF/coeffFunc.f90

uTMDFFFiles=\
$(SOURCEDIR)/Code/uTMDFF/coeffFunc.f90

lpTMDPDFFiles=\
$(SOURCEDIR)/Code/lpTMDPDF/coeffFunc.f90

SiversTMDPDFFiles=\
$(SOURCEDIR)/Code/SiversTMDPDF/placeHolder.f90

wgtTMDPDFFiles=\
$(SOURCEDIR)/Code/wgtTMDPDF/coeffFunc.f90 \
$(SOURCEDIR)/Code/wgtTMDPDF/coeffFunc_largeX.f90

wglTMDPDFFiles=\
$(SOURCEDIR)/Code/wglTMDPDF/coeffFunc.f90 \
$(SOURCEDIR)/Code/wglTMDPDF/coeffFunc_largeX.f90

TMDFFiles=\
$(SOURCEDIR)/Code/TMDF/Fourier_byOgata.f90

TMDKPCFiles=\
$(SOURCEDIR)/Code/TMDF_KPC/TMDpairs.f90\
$(SOURCEDIR)/Code/TMDF_KPC/KERNELpairs_DY.f90\
$(SOURCEDIR)/Code/TMDF_KPC/KERNELpairs_SIDIS.f90

aTMDeSetupFiles=\
$(SOURCEDIR)/Code/aTMDe_setup/placeHolder.f90

TMDXFiles=\
$(SOURCEDIR)/Code/TMDX/DYcoeff-func.f90

aTMDeMODEL = \
$(SOURCEDIR)/Model/TMDR_model.f90 \
$(SOURCEDIR)/Model/uTMDFF_model.f90 \
$(SOURCEDIR)/Model/uTMDPDF_model.f90 

aTMDeOBJ = \
$(OBJ)/aTMDe_Numerics.o \
$(OBJ)/IO_functions.o \
$(OBJ)/IntegrationRoutines.o \
$(OBJ)/InverseMatrix.o \
$(OBJ)/LeptonCutsDY.o \
$(OBJ)/aTMDe_setup.o \
$(OBJ)/LHA_alpha.o \
$(OBJ)/QCDinput.o \
$(OBJ)/EWinput.o\
$(OBJ)/TMD_AD.o\
$(OBJ)/TMDR_model.o\
$(OBJ)/TMDR.o\
$(OBJ)/uTMDPDF_model.o \
$(OBJ)/uTMDPDF_OPE.o \
$(OBJ)/uTMDPDF.o \
$(OBJ)/uTMDFF_model.o \
$(OBJ)/uTMDFF_OPE.o \
$(OBJ)/uTMDFF.o\
$(OBJ)/lpTMDPDF_model.o \
$(OBJ)/lpTMDPDF_OPE.o \
$(OBJ)/lpTMDPDF.o \
$(OBJ)/SiversTMDPDF_model.o \
$(OBJ)/SiversTMDPDF_OPE.o \
$(OBJ)/SiversTMDPDF.o \
$(OBJ)/wgtTMDPDF_model.o \
$(OBJ)/wgtTMDPDF_OPE.o \
$(OBJ)/wgtTMDPDF.o \
$(OBJ)/wglTMDPDF_model.o \
$(OBJ)/wglTMDPDF_OPE.o \
$(OBJ)/wglTMDPDF.o \
$(OBJ)/BoerMuldersTMDPDF_model.o \
$(OBJ)/BoerMuldersTMDPDF_OPE.o \
$(OBJ)/BoerMuldersTMDPDF.o \
$(OBJ)/eeTMDFF_model.o \
$(OBJ)/eeTMDFF.o \
$(OBJ)/TMDF.o \
$(OBJ)/TMDX_DY.o \
$(OBJ)/TMDX_SIDIS.o \
$(OBJ)/TMDF_KPC.o \
$(OBJ)/aTMDe_control.o 

#these are utility object needed to compale any artemide module
aTMDeUTILITY = \
$(OBJ)/aTMDe_Numerics.o \
$(OBJ)/IO_functions.o \
$(OBJ)/IntegrationRoutines.o\
$(OBJ)/InverseMatrix.o


################################################################### COMPILATION OF ARTEMIDE ####################################
FC=$(FCompilator) $(Fflags)

DUMMY=$(shell shopt -s nullglob)

.PHONY: clean default obj program test harpy harpy-signature

.NOTPARALLEL:

default: obj

update: $(BIN)/update-const
	./bin/update-const $(TARGET)


obj: $(aTMDeOBJ) $(aTMDeFILES) $(aTMDeMODEL) $(Twist2Files) $(TMD_ADFiles) $(TMDRFiles) $(uTMDPDFFiles) $(uTMDFFFiles) \
 $(lpTMDFFFiles) $(SiversTMDPDFFiles) $(wgtTMDPDFFiles) $(wglTMDPDFFiles) $(TMDKPCFiles)




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

$(OBJ)/InverseMatrix.o: $(SOURCEDIR)/Code/InverseMatrix.f90 $(OBJ)/aTMDe_Numerics.o $(OBJ)/IO_functions.o
	$(FC) -c $(SOURCEDIR)/Code/InverseMatrix.f90 -I$(MOD)
	mv *.o $(OBJ)
	mv *.mod $(MOD)

$(OBJ)/LeptonCutsDY.o: $(SOURCEDIR)/LeptonCutsDY.f90 $(aTMDeUTILITY)
	$(FC) -c $(SOURCEDIR)/LeptonCutsDY.f90 -I$(MOD)
	mv *.o $(OBJ)
	mv *.mod $(MOD)

$(OBJ)/LHA_alpha.o: $(SOURCEDIR)/Code/LHA/LHA_alpha.f90 $(aTMDeUTILITY)
#	mkdir -p obj
	$(FC) -c $(SOURCEDIR)/Code/LHA/LHA_alpha.f90 -I$(MOD)
	mv *.o $(OBJ)
	mv *.mod $(MOD)

$(OBJ)/QCDinput.o: $(SOURCEDIR)/QCDinput.f90 $(OBJ)/LHA_alpha.o $(aTMDeUTILITY) $(SOURCEDIR)/Code/LHA/LHA_PDF.f90
#	mkdir -p obj
	$(FC) -c $(SOURCEDIR)/QCDinput.f90 -I$(MOD) $(FOPT)
	mv *.o $(OBJ)
	mv *.mod $(MOD)
	
$(OBJ)/EWinput.o: $(SOURCEDIR)/EWinput.f90 $(aTMDeUTILITY)
#	mkdir -p obj
	$(FC) -c $(SOURCEDIR)/EWinput.f90 -I$(MOD)
	mv *.o $(OBJ)
	mv *.mod $(MOD)

$(OBJ)/TMD_AD.o: $(SOURCEDIR)/TMD_AD.f90 $(aTMDeUTILITY) $(TMD_ADFiles) $(OBJ)/QCDinput.o
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

$(OBJ)/uTMDPDF_model.o: $(SOURCEDIR)/Model/uTMDPDF_model.f90 $(aTMDeUTILITY)
#	mkdir -p obj
	$(FC) -c $(SOURCEDIR)/Model/uTMDPDF_model.f90 -I$(MOD)
	mv *.o $(OBJ)
	mv *.mod $(MOD)

$(OBJ)/uTMDPDF_OPE.o: $(SOURCEDIR)/uTMDPDF_OPE.f90 $(OBJ)/QCDinput.o $(SOURCEDIR)/Model/uTMDPDF_model.f90 $(Twist2Files) $(aTMDeUTILITY) $(uTMDPDFFiles)
#	mkdir -p obj
	$(FC) -c $(SOURCEDIR)/uTMDPDF_OPE.f90 -I$(MOD)
	mv *.o $(OBJ)
	mv *.mod $(MOD)
	
$(OBJ)/uTMDPDF.o: $(SOURCEDIR)/uTMDPDF.f90 $(OBJ)/QCDinput.o $(OBJ)/TMDR.o $(SOURCEDIR)/Model/uTMDPDF_model.f90 $(SOURCEDIR)/uTMDPDF_OPE.f90 $(KTspaceFiles) $(aTMDeUTILITY) $(uTMDPDFFiles)
#	mkdir -p obj
	$(FC) -c $(SOURCEDIR)/uTMDPDF.f90 -I$(MOD)
	mv *.o $(OBJ)
	mv *.mod $(MOD)

$(OBJ)/uTMDFF_model.o: $(SOURCEDIR)/Model/uTMDFF_model.f90 $(aTMDeUTILITY)
#	mkdir -p obj
	$(FC) -c $(SOURCEDIR)/Model/uTMDFF_model.f90 -I$(MOD)
	mv *.o $(OBJ)
	mv *.mod $(MOD)

$(OBJ)/uTMDFF_OPE.o: $(SOURCEDIR)/uTMDFF_OPE.f90 $(OBJ)/QCDinput.o $(SOURCEDIR)/Model/uTMDFF_model.f90 $(Twist2Files) $(aTMDeUTILITY) $(uTMDFFFiles)
#	mkdir -p obj
	$(FC) -c $(SOURCEDIR)/uTMDFF_OPE.f90 -I$(MOD)
	mv *.o $(OBJ)
	mv *.mod $(MOD)
	
$(OBJ)/uTMDFF.o: $(SOURCEDIR)/uTMDFF.f90 $(OBJ)/QCDinput.o $(OBJ)/TMDR.o $(SOURCEDIR)/Model/uTMDFF_model.f90 $(SOURCEDIR)/uTMDFF_OPE.f90 $(KTspaceFiles) $(aTMDeUTILITY) $(uTMDFFFiles)
#	mkdir -p obj
	$(FC) -c $(SOURCEDIR)/uTMDFF.f90 -I$(MOD)
	mv *.o $(OBJ)
	mv *.mod $(MOD)
	
$(OBJ)/lpTMDPDF_model.o: $(SOURCEDIR)/Model/lpTMDPDF_model.f90 $(aTMDeUTILITY)
#	mkdir -p obj
	$(FC) -c $(SOURCEDIR)/Model/lpTMDPDF_model.f90 -I$(MOD)
	mv *.o $(OBJ)
	mv *.mod $(MOD)

$(OBJ)/lpTMDPDF_OPE.o: $(SOURCEDIR)/lpTMDPDF_OPE.f90 $(OBJ)/QCDinput.o $(OBJ)/TMDR.o $(SOURCEDIR)/Model/lpTMDPDF_model.f90 $(Twist2Files) $(aTMDeUTILITY) $(uTMDPDFFiles)
#	mkdir -p obj
	$(FC) -c $(SOURCEDIR)/lpTMDPDF_OPE.f90 -I$(MOD)
	mv *.o $(OBJ)
	mv *.mod $(MOD)

$(OBJ)/lpTMDPDF.o: $(SOURCEDIR)/lpTMDPDF.f90 $(OBJ)/QCDinput.o $(SOURCEDIR)/Model/lpTMDPDF_model.f90 $(SOURCEDIR)/lpTMDPDF_OPE.f90 $(KTspaceFiles) $(aTMDeUTILITY) $(lpTMDPDFFiles)
#	mkdir -p obj
	$(FC) -c $(SOURCEDIR)/lpTMDPDF.f90 -I$(MOD)
	mv *.o $(OBJ)
	mv *.mod $(MOD)

$(OBJ)/SiversTMDPDF_model.o: $(SOURCEDIR)/Model/SiversTMDPDF_model.f90 $(aTMDeUTILITY)
#	mkdir -p obj
	$(FC) -c $(SOURCEDIR)/Model/SiversTMDPDF_model.f90 -I$(MOD)
	mv *.o $(OBJ)
	mv *.mod $(MOD)

$(OBJ)/SiversTMDPDF_OPE.o: $(SOURCEDIR)/SiversTMDPDF_OPE.f90 $(OBJ)/QCDinput.o $(SOURCEDIR)/Model/SiversTMDPDF_model.f90 $(Twist3Files) $(aTMDeUTILITY) $(SiversTMDPDFFiles)
#	mkdir -p obj
	$(FC) -c $(SOURCEDIR)/SiversTMDPDF_OPE.f90 -I$(MOD)
	mv *.o $(OBJ)
	mv *.mod $(MOD)

$(OBJ)/SiversTMDPDF.o: $(SOURCEDIR)/SiversTMDPDF.f90 $(OBJ)/QCDinput.o $(OBJ)/TMDR.o $(SOURCEDIR)/Model/SiversTMDPDF_model.f90 $(SOURCEDIR)/SiversTMDPDF_OPE.f90 $(KTspaceFiles) $(aTMDeUTILITY) $(SiversTMDPDFFiles)
#	mkdir -p obj
	$(FC) -c $(SOURCEDIR)/SiversTMDPDF.f90 -I$(MOD)
	mv *.o $(OBJ)
	mv *.mod $(MOD)

$(OBJ)/wgtTMDPDF_model.o: $(SOURCEDIR)/Model/wgtTMDPDF_model.f90 $(aTMDeUTILITY)
#	mkdir -p obj
	$(FC) -c $(SOURCEDIR)/Model/wgtTMDPDF_model.f90 -I$(MOD)
	mv *.o $(OBJ)
	mv *.mod $(MOD)

$(OBJ)/wgtTMDPDF_OPE.o: $(SOURCEDIR)/wgtTMDPDF_OPE.f90 $(OBJ)/QCDinput.o $(SOURCEDIR)/Model/wgtTMDPDF_model.f90 $(Twist2Files) $(Twist3Files) $(aTMDeUTILITY) $(uTMDPDFFiles)
#	mkdir -p obj
	$(FC) -c $(SOURCEDIR)/wgtTMDPDF_OPE.f90 -I$(MOD)
	mv *.o $(OBJ)
	mv *.mod $(MOD)

$(OBJ)/wgtTMDPDF.o: $(SOURCEDIR)/wgtTMDPDF.f90 $(OBJ)/QCDinput.o $(OBJ)/TMDR.o $(SOURCEDIR)/Model/wgtTMDPDF_model.f90 $(SOURCEDIR)/wgtTMDPDF_OPE.f90 $(KTspaceFiles) $(aTMDeUTILITY) $(wgtTMDPDFFiles)
#	mkdir -p obj
	$(FC) -c $(SOURCEDIR)/wgtTMDPDF.f90 -I$(MOD)
	mv *.o $(OBJ)
	mv *.mod $(MOD)

$(OBJ)/wglTMDPDF_model.o: $(SOURCEDIR)/Model/wglTMDPDF_model.f90 $(aTMDeUTILITY)
#	mkdir -p obj
	$(FC) -c $(SOURCEDIR)/Model/wglTMDPDF_model.f90 -I$(MOD)
	mv *.o $(OBJ)
	mv *.mod $(MOD)

$(OBJ)/wglTMDPDF_OPE.o: $(SOURCEDIR)/wglTMDPDF_OPE.f90 $(OBJ)/QCDinput.o $(SOURCEDIR)/Model/wglTMDPDF_model.f90 $(Twist2Files) $(Twist3Files) $(aTMDeUTILITY) $(uTMDPDFFiles) $(wglTMDPDFFiles)
#	mkdir -p obj
	$(FC) -c $(SOURCEDIR)/wglTMDPDF_OPE.f90 -I$(MOD)
	mv *.o $(OBJ)
	mv *.mod $(MOD)

$(OBJ)/wglTMDPDF.o: $(SOURCEDIR)/wglTMDPDF.f90 $(OBJ)/QCDinput.o $(OBJ)/TMDR.o $(SOURCEDIR)/Model/wglTMDPDF_model.f90 $(SOURCEDIR)/wglTMDPDF_OPE.f90 $(KTspaceFiles) $(aTMDeUTILITY) $(wglTMDPDFFiles)
#	mkdir -p obj
	$(FC) -c $(SOURCEDIR)/wglTMDPDF.f90 -I$(MOD)
	mv *.o $(OBJ)
	mv *.mod $(MOD)

$(OBJ)/BoerMuldersTMDPDF_model.o: $(SOURCEDIR)/Model/BoerMuldersTMDPDF_model.f90 $(aTMDeUTILITY)
#	mkdir -p obj
	$(FC) -c $(SOURCEDIR)/Model/BoerMuldersTMDPDF_model.f90 -I$(MOD)
	mv *.o $(OBJ)
	mv *.mod $(MOD)

$(OBJ)/BoerMuldersTMDPDF_OPE.o: $(SOURCEDIR)/BoerMuldersTMDPDF_OPE.f90 $(OBJ)/QCDinput.o $(SOURCEDIR)/Model/BoerMuldersTMDPDF_model.f90 $(Twist3Files) $(aTMDeUTILITY) $(BoerMuldersTMDPDFFiles)
#	mkdir -p obj
	$(FC) -c $(SOURCEDIR)/BoerMuldersTMDPDF_OPE.f90 -I$(MOD)
	mv *.o $(OBJ)
	mv *.mod $(MOD)

$(OBJ)/BoerMuldersTMDPDF.o: $(SOURCEDIR)/BoerMuldersTMDPDF.f90 $(OBJ)/QCDinput.o $(OBJ)/TMDR.o $(SOURCEDIR)/Model/BoerMuldersTMDPDF_model.f90 $(SOURCEDIR)/BoerMuldersTMDPDF_OPE.f90 $(KTspaceFiles) $(aTMDeUTILITY) $(BoerMuldersTMDPDFFiles)
#	mkdir -p obj
	$(FC) -c $(SOURCEDIR)/BoerMuldersTMDPDF.f90 -I$(MOD)
	mv *.o $(OBJ)
	mv *.mod $(MOD)

$(OBJ)/eeTMDFF_model.o: $(SOURCEDIR)/Model/eeTMDFF_model.f90 $(aTMDeUTILITY)
#	mkdir -p obj
	$(FC) -c $(SOURCEDIR)/Model/eeTMDFF_model.f90 -I$(MOD)
	mv *.o $(OBJ)
	mv *.mod $(MOD)

$(OBJ)/eeTMDFF.o: $(SOURCEDIR)/eeTMDFF.f90 $(OBJ)/QCDinput.o $(OBJ)/TMDR.o $(SOURCEDIR)/Model/eeTMDFF_model.f90 $(aTMDeUTILITY)
#	mkdir -p obj
	$(FC) -c $(SOURCEDIR)/eeTMDFF.f90 -I$(MOD)
	mv *.o $(OBJ)
	mv *.mod $(MOD)

$(OBJ)/TMDF.o: $(SOURCEDIR)/TMDF.f90 $(TMDFFiles) $(OBJ)/EWinput.o $(OBJ)/uTMDPDF.o $(OBJ)/uTMDFF.o $(OBJ)/lpTMDPDF.o $(OBJ)/SiversTMDPDF.o $(OBJ)/wgtTMDPDF.o $(OBJ)/BoerMuldersTMDPDF.o $(OBJ)/wglTMDPDF.o $(OBJ)/eeTMDFF.o $(aTMDeUTILITY)
#	mkdir -p obj
	$(FC) -c $(SOURCEDIR)/TMDF.f90 -I$(MOD)
	mv *.o $(OBJ)
	mv *.mod $(MOD)

$(OBJ)/TMDF_KPC.o: $(SOURCEDIR)/TMDF_KPC.f90 $(OBJ)/EWinput.o $(OBJ)/uTMDPDF.o $(OBJ)/uTMDFF.o $(OBJ)/lpTMDPDF.o $(OBJ)/SiversTMDPDF.o $(OBJ)/wgtTMDPDF.o $(OBJ)/BoerMuldersTMDPDF.o $(OBJ)/wglTMDPDF.o $(OBJ)/eeTMDFF.o $(TMDKPCFiles) $(aTMDeUTILITY)
#	mkdir -p obj
	$(FC) -c $(SOURCEDIR)/TMDF_KPC.f90 -I$(MOD)
	mv *.o $(OBJ)
	mv *.mod $(MOD)

$(OBJ)/TMDX_DY.o: $(SOURCEDIR)/TMDX_DY.f90 $(TMDXFiles) $(OBJ)/TMDF.o $(OBJ)/TMDF_KPC.o  $(OBJ)/QCDinput.o $(aTMDeUTILITY)
#	mkdir -p obj
	$(FC) -c $(SOURCEDIR)/TMDX_DY.f90 -I$(MOD)
	mv *.o $(OBJ)
	mv *.mod $(MOD)
	
$(OBJ)/TMDX_SIDIS.o: $(SOURCEDIR)/TMDX_SIDIS.f90 $(OBJ)/TMDF.o $(OBJ)/QCDinput.o $(aTMDeUTILITY)
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
	f2py -c --f90exec=$(Fpath) --f77exec=$(F77path) --f90flags=$(FflagsHARPY) $(FOPT) -lgomp -I$(MOD) $(aTMDeFILES) $(HDIR)/harpy.f90 $(HDIR)/artemide.pyf
	mv artemide*.so $(HDIR)
