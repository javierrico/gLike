TARGET     = gLike
PWD        = `pwd`
SRCDIR     = ./src
INCLDIR    = ./include
OUTDIR     = ./out
LIBDIR     = ./lib
CXX        = g++ -fPIC
ROOTVERSION=`root-config --version | tr -c -d [:digit:]`
MINVERSION = 62004
ROOTLIBS   = `root-config --libs` -lMinuit
# add the FITSIO flag for compilation if a version greater than 6.20/04 is found 
ifeq ($(shell test $(ROOTVERSION) -gt $(MINVERSION); echo $$?), 0)
	ROOTLIBS += -lFitsio
endif
ROOTCFLAGS = `root-config --cflags`
SRCFILES   = Lkl ParabolaLkl PoissonLkl JointLkl Iact1dUnbinnedLkl Iact1dBinnedLkl IactEventListIrf FermiTables2016Lkl TemplateLkl MIACTEventListIRF GloryDuckTables2019Lkl
SOURCES    = $(SRCFILES:%=$(SRCDIR)/%.cc)
HEADERS    = $(SRCFILES:%=$(INCLDIR)/%.h)
OBJECTS    = $(SRCFILES:%=$(OUTDIR)/%.o)
DICTSRC    = $(OUTDIR)/$(TARGET)Dict.cc
LINKDEF    = $(INCLDIR)/$(TARGET)LinkDef.h

$(LIBDIR)/lib$(TARGET).so: $(OBJECTS) $(DICTSRC:%.cc=%.o)
	@echo "Linking shared object $@..."
	@mkdir -p lib
	$(CXX) -shared $^ $(ROOTLIBS) -o $@
	@echo "Done!"
	@echo $(ROOTLIBS)

$(DICTSRC:%.cc=%.o): $(DICTSRC) $(HEADERS)
	@echo "Compiling $<..."
	$(CXX) $(ROOTCFLAGS) -c $< -o $@
	@echo "Generated $@ \n"

$(DICTSRC): $(HEADERS) $(LINKDEF)
	@echo "Generating dictionary $(@:out/%.cc=%)..."
	$(ROOTSYS)/bin/rootcint -f $@ -c -p $(HEADERS:%=$(PWD)/%) $(PWD)/$(LINKDEF)
	@echo "Generated $@ and $(@:.cc=.h) \n"

$(OUTDIR)/%.o: $(SRCDIR)/%.cc $(HEADERS)
	@echo "Compiling $<..."
	@mkdir -p out 
	$(CXX) -I$(INCLDIR) $(ROOTCFLAGS) -c $< -o $@
	@echo "Generated $@ \n"

.PHONY : clean
clean:
	@rm -f $(OUTDIR)/*.o
	@rm -f $(OUTDIR)/$(TARGET)Dict.*
	@rm -f $(LIBDIR)/lib$(TARGET).so
	@rm -rf htmldoc/* 2>&1

.PHONY : doc
doc: $(LIBDIR)/lib$(TARGET).so
	@echo "\nCreating html documentation and logfile dohtml.log..."
	@mkdir -p htmldoc 
	@rm -rf htmldoc/*
	root -b -q scripts/dohtml.C 2>&1 > htmldoc/dohtml.log
	@echo "Done!"

