################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../src/adaptivity.cpp \
../src/body.cpp \
../src/cont_mech.cpp \
../src/contact.cpp \
../src/correctors.cpp \
../src/derivatives.cpp \
../src/grid.cpp \
../src/johnson_cook_Sima_2010.cpp \
../src/kernel.cpp \
../src/leap_frog.cpp \
../src/logger.cpp \
../src/material.cpp \
../src/particle.cpp \
../src/plasticity.cpp \
../src/precomp_shape_functions.cpp \
../src/refine_cut_main.cpp \
../src/simulation_data.cpp \
../src/simulation_time.cpp \
../src/thermal.cpp \
../src/tool.cpp \
../src/vtk_writer.cpp 

OBJS += \
./src/adaptivity.o \
./src/body.o \
./src/cont_mech.o \
./src/contact.o \
./src/correctors.o \
./src/derivatives.o \
./src/grid.o \
./src/johnson_cook_Sima_2010.o \
./src/kernel.o \
./src/leap_frog.o \
./src/logger.o \
./src/material.o \
./src/particle.o \
./src/plasticity.o \
./src/precomp_shape_functions.o \
./src/refine_cut_main.o \
./src/simulation_data.o \
./src/simulation_time.o \
./src/thermal.o \
./src/tool.o \
./src/vtk_writer.o 

CPP_DEPS += \
./src/adaptivity.d \
./src/body.d \
./src/cont_mech.d \
./src/contact.d \
./src/correctors.d \
./src/derivatives.d \
./src/grid.d \
./src/johnson_cook_Sima_2010.d \
./src/kernel.d \
./src/leap_frog.d \
./src/logger.d \
./src/material.d \
./src/particle.d \
./src/plasticity.d \
./src/precomp_shape_functions.d \
./src/refine_cut_main.d \
./src/simulation_data.d \
./src/simulation_time.d \
./src/thermal.d \
./src/tool.d \
./src/vtk_writer.d 


# Each subdirectory must supply rules for building sources it contributes
src/%.o: ../src/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -std=c++1y -O0 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


