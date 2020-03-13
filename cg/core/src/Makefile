
.DEFAULT_GOAL := default

CXX = icpc

STATICOPTS = -static -static-intel -static-libgcc -static-libstdc++

CXXFLAGS = -O2 -I${GSL_DIR}/include -DUSE_MKL -fPIC
LDFLAGS = -lgsl -lgslcblas -mkl -lmkl_gf_lp64 -lmkl_intel_thread -lmkl_core -liomp5 $(STATICOPTS)

SRC = $(wildcard *.cpp)
DEP = $(SRC:.cpp=.d)
OBJ = $(SRC:.cpp=.o)

EXE = fm.x

default: $(EXE)
	@echo 'Test finished'

%.o: %.cpp
	@echo '=>' Compile $@ ...
	@$(CXX) $(CXXFLAGS) $(CXXFLAGS_EXT) -c $<

%.d: %.cpp
	@echo '=>' Create Dependence $@ ...
	@$(CXX) $(CXXFLAGS) $(CXXFLAGS_EXT) -M $< > $@

$(EXE): $(OBJ)
	@echo '=>' Link $@ ...
	@$(CXX) $(OBJ) $(LDFLAGS) -o $(EXE)

clean:
	@/bin/rm -rf *.o *.d *.x

test:
	time -p ./$(EXE)

SHFLAGS= -shared -fPIC

matrix.so: matrix.o
	$(CXX) $(SHFLAGS) $(LDFLAGS) matrix.o -o matrix.so

sinclude $(DEP)
