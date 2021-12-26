CXX = g++
TARGET = fast_fourier_transfrom_call_pricer
HEADERS = characteristic_functions.h global_variables.h utils.h
CFLAG = -O3 -Wall -shared -std=c++17 -fPIC
#PYINCLUDE = `python3 -m pybind11 --includes`
PYINCLUDE = `python3-config --includes` -Iextern/pybind11/include
PYSUFFIX = `python3-config --extension-suffix`

.PHONY: run check clean

all: $(TARGET)

$(TARGET): $(TARGET).cpp 
	$(CXX) $(CFLAG) $(PYINCLUDE) $< -o $@$(PYSUFFIX)

test: 
	python3 -m pytest -v

clean:
	rm -rf *.so .pytest_cache __pycache__ 
