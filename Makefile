
ROOTCFLAGS = `root-config --cflags`
ROOTLIBS   = `root-config --libs`

CXXFLAGS += -I. -Wall
SHAREDFLAGS =  -shared -Wl



%.o : %.c
	$(RM) $@
	$(CC) -c $(CFLAGS) -o $@ $<
	
%.o : %.cc
	$(RM) $@
	$(CXX) -c $(CXXFLAGS) $(ROOTCFLAGS) -o $@ $*.cc

WRAPPERS = py_wrapper.o

OBJS    = EarthDensity.o BargerPropagator.o mosc.o mosc3.o FullSMEPropagator.o \
          $(WRAPPERS)

LIBBASE   = ThreeProb
VER       = 3.10
TAG       = 
LIBALIAS  = $(LIBBASE)$(TAG)
LIBNAME   = $(LIBALIAS)_$(VER)

lib3p     = lib$(LIBNAME).a
lib3ps    = lib$(LIBNAME).so
LINK      = lib$(LIBBASE).so


targets = $(lib3p) probRoot probLinear probAnalytic


$(lib3p) : $(OBJS) 
	$(RM) $@
	ar clq $@ $(OBJS)
	ranlib $@


$(lib3ps) : $(OBJS)
	$(RM) $@
	$(CXX) $(SHAREDFLAGS)$@ -o $@ $(OBJS) 
	$(RM) $(LINK)
	ln -s $(lib3ps) $(LINK)

shared : $(lib3ps)
	

probRoot: probRoot.o $(lib3p) 
	$(RM) $@
	$(CXX) -o $@ $(CXXFLAGS) -L. $^ $(ROOTLIBS)


.PHONY: probRoot.o
probRoot.o: 
	$(CXX) -o probRoot.o $(ROOTCFLAGS) $(CXXFLAGS) -c probRoot.cc



probLinear: probLinear.o $(lib3p) 
	$(RM) $@
	$(CXX) -o $@ $(CXXFLAGS) -L. $^ $(ROOTLIBS)


.PHONY: probLinear.o
probLinear.o: 
	$(CXX) -o probLinear.o $(ROOTCFLAGS) $(CXXFLAGS) -c probLinear.cc


probAnalytic: probAnalytic.o $(lib3p) 
	$(RM) $@
	$(CXX) -o $@ $(CXXFLAGS) -L. $^ $(ROOTLIBS)


.PHONY: probAnalytic.o
probAnalytic.o: 
	$(CXX) -o probAnalytic.o $(ROOTCFLAGS) $(CXXFLAGS) -c probAnalytic.cc


probLV: probLV.o $(lib3p) 
	$(RM) $@
	$(CXX) -o $@ $(CXXFLAGS) -L. $^ $(ROOTLIBS)


.PHONY: probLV.o
probLV.o: 
	$(CXX) -o probLV.o $(ROOTCFLAGS) $(CXXFLAGS) -c probLV.cc

	
.PHONY: all
all: $(targets)
	
	
.PHONY: clean
clean:
	$(RM) $(targets) *.o *.so
	
emptyrule:: $(lib3p)

	
