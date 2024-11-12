CXXFLAGS += -I. $(shell root-config --cflags) -g
LDFLAGS += $(shell root-config --libs) -lPhysics -g

PROGRAMS = WaveformStudy ProcessWaveformFiles

all:	clean $(PROGRAMS)

$(PROGRAMS):
	@echo '<<compiling' $@'>>'
	@$(CXX) $@.cpp -o $@ $(CXXFLAGS) $(LDFLAGS)
	@rm -rf *.dSYM
clean:	
	rm -f $(PROGRAMS)
