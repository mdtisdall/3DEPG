
MKLROOT = /opt/intel/mkl


CPPFLAGS += -DMKL_ILP64 -m64 
CPPFLAGS += -I${MKLROOT}/include 

CXXFLAGS += -g

# Flags passed to the C++ compiler.
CXXFLAGS += -Wall -Wextra -pthread -std=c++11
CXXFLAGS += -DMKL_DIRECT_CALL_SEQ

#LDFLAGS += -L ~/lib
LDLIBS += ${MKLROOT}/lib/libmkl_intel_ilp64.a
LDLIBS += ${MKLROOT}/lib/libmkl_sequential.a
LDLIBS += ${MKLROOT}/lib/libmkl_core.a
LDFLAGS += -lpthread -lm -ldl

all : test transtest 

transtest :

test : blas_local_MKL.o

clean:
	rm -rf test
	rm -rf transtest
