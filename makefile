FC = gfortran
LD = gfortran

.SUFFIXES:.f90

TARGET = main

OBJ = ./obj/stdio.o\
			./obj/inputCtrl.o\
			./obj/inputMesh.o\
			./obj/inputMediumProperty.o\
			./obj/inputFluidProperty.o\
			./obj/physicalConstant.o\
			./obj/molecularProps.o\
			./obj/cardano.o\
			./obj/meshObj.o\
			./obj/diffusion.o\
			./obj/Thomas.o\
			./obj/matrixCalcTest.o\
			./obj/main.o

all: $(TARGET)

./obj/%.o: ./src/%.f90
	@mkdir -p obj
	cd obj &&$(FC) -o ./$*.o -c ../$< &&cd ../

%.mod: %.f90 %.o
  @:

$(TARGET):$(OBJ)
	@mkdir -p bin
	$(LD) -o ./bin/$@ $(OBJ)

clean:
	rm -rf ./obj *~ core

go:
	make clean && make && bin/main

plot:
	Rscript output.R && open Rplots.pdf
