CXXFLAGS += -DMKL_ILP64 -m64 

CXX = g++

CXXFLAGS += -g
#CXXFLAGS += -O3

# Flags passed to the C++ compiler.
CXXFLAGS += -Wall -Wextra -pthread -std=c++14
CXXFLAGS += -DMKL_DIRECT_CALL_SEQ
CXXFLAGS += -I./

#LDFLAGS += -L ~/lib
ifeq ($(strip $(OS)),Linux)
	LDLIBS += -Wl,--start-group libmkl_intel_ilp64.a libmkl_sequential.a libmkl_core.a -Wl,--end-group
	LDLIBS += -lpthread -lm -ldl
else
	LDLIBS += -lmkl_intel_ilp64
	LDLIBS += -lmkl_sequential
	LDLIBS += -lmkl_core
	LDLIBS += -lpthread -lm -ldl
endif

all : examples/MPSSFP_example

examples/MPSSFP_example : core/blas_local_MKL.o core/mkl_ops_mkl.o

clean:
	rm -f core/*.o
	rm -rf examples/MPSSFP_example
