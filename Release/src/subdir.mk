################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CC_SRCS += \
../src/BlockCube.cc \
../src/DataCube.cc \
../src/FragCube.cc \
../src/Handler.cc \
../src/Main.cc 

CC_DEPS += \
./src/BlockCube.d \
./src/DataCube.d \
./src/FragCube.d \
./src/Handler.d \
./src/Main.d 

OBJS += \
./src/BlockCube.o \
./src/DataCube.o \
./src/FragCube.o \
./src/Handler.o \
./src/Main.o 


# Each subdirectory must supply rules for building sources it contributes
src/%.o: ../src/%.cc src/subdir.mk
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	mpic++ -I/home/afmorielo/Codebase/ITA/Eclipse/datacube-mpi-framework/include -I/home/afmorielo/Softwares/Bibliotecas/BoostC++/1_71_0/include -O3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$@" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


clean: clean-src

clean-src:
	-$(RM) ./src/BlockCube.d ./src/BlockCube.o ./src/DataCube.d ./src/DataCube.o ./src/FragCube.d ./src/FragCube.o ./src/Handler.d ./src/Handler.o ./src/Main.d ./src/Main.o

.PHONY: clean-src

