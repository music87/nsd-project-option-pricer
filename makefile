CXX = g++
SRCDIR = src
TESTDIR = test
TARGET = fast_fourier_transfrom_call_pricer
HEADERS = $(SRCDIR)/characteristic_functions.h $(SRCDIR)/global_variables.h $(SRCDIR)/utils.h
CFLAG = -O3 -Wall -shared -std=c++17 -fPIC
#PYINCLUDE = `python3 -m pybind11 --includes`
PYINCLUDE = `python3-config --includes` -Iextern/pybind11/include
PYSUFFIX = `python3-config --extension-suffix`

.PHONY: test clean

all: $(TARGET)

$(TESTDIR)/$(TARGET): $(SRCDIR)/$(TARGET).cpp $(HEADERS)
	$(CXX) $(CFLAG) $(PYINCLUDE) $< -o $@$(PYSUFFIX)

test: 
	python3 -m pytest -v

clean:
	rm -rf *.so .pytest_cache __pycache__ 

