# CLASS is required to compile the test code
CLASSROOT=/home/marc/soft/class_v2.0.4
LFLAGS= -lblas -llapack
CXXFLAGS= -O3 -Wall

CLASS_CXXFLAGS= -I$(CLASSROOT)/include/ -Isrc/
CLASS_LFLAGS= -fopenmp $(CLASSROOT)/libclass.a

libjla.a: src/jla.o src/ini.o
	$(AR) $(ARFLAGS) $@ $^

test_jla: src/test.o libjla.a
	$(CXX) -o $@ $^ $(CXXFLAGS) $(LFLAGS) $(CLASS_LFLAGS)

src/test.o: src/test.cc
	$(CXX) -c -o $@ $^  $(CXXFLAGS) $(CLASS_CXXFLAGS)

clean:
	rm -f src/*.o
	rm -f libjla.a
