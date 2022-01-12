CPPFLAGS      = -Wall -ansi -Wextra
CXXFLAGS      = -g -O4
SYSTEM        = $(shell uname -s)
LDLIBS.Linux  = -pthread
LDLIBS.Darwin =
LDLIBS        = $(LDLIBS.$(SYSTEM))
LINK.o        = $(CXX) $(LDFLAGS)

CXXSRC = IntegrateQED.cc QEDCalculation.cc FourVector.cc\
         Spinor.cc QuasiRandom.cc

IntegrateQED: QEDCalculation.o FourVector.o Spinor.o QuasiRandom.o

clean:
	@rm -fv *.o *~ '#'*'#' .depend IntegrateQED

.depend : $(CXXSRC)
	$(CXX) -MM $(CPPFLAGS) $^ > $@

-include .depend
