CXX = mpicc
COMPILER_OPTIONS = -O3 -openmp -I./ -I./ParMETIS/ -I./ParMETIS/METISLib/
LINKER_OPTIONS = -lm -openmp ./ParMETIS/libparmetis.a ./ParMETIS/libmetis.a

TARGET = parallelLimit

SOURCES = cfdDataInitialization.c cfdFieldValuesInitialization.c cfdGeometryInitialization.c csrHmeshDualMPI.c \
          csrPartParMETISMPI.c csrProcessCuthillMcKee.c csrProcessCuthillMcKeeMultiLevel.c hm4DualGraphCSR.c \
		  hm4ElementsCenters.c hm4ElementsFaces.c hm4ElementsHeights.c hm4ElementsVolumes.c \
		  hm4MeshPreprocessorCfdMPI.c hm4PreprocessorMPI.c ioHM3D.c iohm3DMPI.c main.c \
		  mpiCfdCore.c mpiCfdCoreMain.c mpiTransferHost.c processIntegerLIST.c \
		  solverCUSP.c solverRusanov3D.c sortIcfList.c

OBJECTS = $(patsubst %.c, %.o, $(SOURCES))

%.o : %.c
	@echo "compiling $^ ..."
	@$(CXX) $(COMPILER_OPTIONS) -c $^ -o $@

all : $(OBJECTS)
	@echo "building $(TARGET) ..."
	@$(CXX) $(OBJECTS) $(LINKER_OPTIONS) -o $(TARGET).exe

clean :
	@rm -f $(TARGET).exe $(OBJECTS)
