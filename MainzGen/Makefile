#SETUPS = EIC1x1 EIC1x5 EIC5x5 EIC5x1 EIC25x5 EIC100x5 EIC250x5 
SETUPS = EIC275x18d

#DL45 DL32  DL45sym

JOBS = 4
EVENTS=1e7

all : 	src/IntegrateQED $(SETUPS:%=Integration_%.pdf)  $(SETUPS:%=Signal_%.pdf) 


%.pdf : %.dat %.gnuplot
	gnuplot $*.gnuplot

Integration_%.dat Integration_%.gnuplot: src/IntegrateQED setups/%
	src/IntegrateQED -p Integration_ -f setups/$* -j $(JOBS) -e $(EVENTS) 

Signal_%.dat Signal_%.gnuplot: src/IntegrateQED setups/%
	ANTISYM=0 TRIDENT=0 \
	src/IntegrateQED -p Signal_ -f setups/$* -j $(JOBS) -e $(EVENTS) 
#	./compare.sh $*

src/IntegrateQED: src/IntegrateQED.cc
	$(MAKE) -C $(dir $@) $(notdir $@)

.PRECIOUS: Integration_%.dat Integration_%.gnuplot Signal_%.dat Singal_%.gnuplot

clean:
	@ rm -fv *~ '#'* *.dat *.gnuplot *.ps
	@ $(MAKE) -C src clean
