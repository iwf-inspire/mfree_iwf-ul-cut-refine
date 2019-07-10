################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../src/benchmarks/material_library.cpp \
../src/benchmarks/test_benches.cpp \
../src/benchmarks/test_cuttings.cpp \
../src/benchmarks/test_density.cpp 

OBJS += \
./src/benchmarks/material_library.o \
./src/benchmarks/test_benches.o \
./src/benchmarks/test_cuttings.o \
./src/benchmarks/test_density.o 

CPP_DEPS += \
./src/benchmarks/material_library.d \
./src/benchmarks/test_benches.d \
./src/benchmarks/test_cuttings.d \
./src/benchmarks/test_density.d 


# Each subdirectory must supply rules for building sources it contributes
src/benchmarks/%.o: ../src/benchmarks/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -std=c++1y -O0 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


