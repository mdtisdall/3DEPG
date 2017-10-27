
ifeq ($(strip $(OS)),Linux)
MKLROOT = $(HOME)/intel/mkl
else
MKLROOT = /opt/intel/mkl
endif

CPPFLAGS += -DMKL_ILP64 -m64 
CPPFLAGS += -I${MKLROOT}/include 

#CXXFLAGS += -g
CXXFLAGS += -O3

# Flags passed to the C++ compiler.
CXXFLAGS += -Wall -Wextra -pthread -std=c++14
CXXFLAGS += -DMKL_DIRECT_CALL_SEQ
CXXFLAGS += -I./

#LDFLAGS += -L ~/lib
ifeq ($(strip $(OS)),Linux)
	LDLIBS += -Wl,--start-group ${MKLROOT}/lib/intel64/libmkl_intel_ilp64.a ${MKLROOT}/lib/intel64/libmkl_sequential.a ${MKLROOT}/lib/intel64/libmkl_core.a -Wl,--end-group
	LDLIBS += -lpthread -lm -ldl
else
	LDLIBS += ${MKLROOT}/lib/libmkl_intel_ilp64.a
	LDLIBS += ${MKLROOT}/lib/libmkl_sequential.a
	LDLIBS += ${MKLROOT}/lib/libmkl_core.a
	LDLIBS += -lpthread -lm -ldl
endif

all : test transtest 

transtest :

test : blas_local_MKL.o mkl_ops_mkl.o

clean:
	rm *.o
	rm -rf test
	rm -rf transtest
