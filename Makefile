#
# MERLIN Makefile -- Compiles and installs MERLIN and accessory applications
# (c) 2000-2007 Goncalo Abecasis
#

# version information
VERSION = 1.1.2
PSVERSION = 0.6.10

# default installation directory
INSTALLDIR=/usr/local/bin

# default C++ compiler
CXX=g++ 

# default compilation flags are 
#
# CFLAGS=-O2 -I./libsrc/ -I./merlin/ -I./pdf/
# 
# The following special options may also be added to the default
# 
#      Option                        Effect
#      -D_FILE_OFFSET_BITS=64        Enables support for swapfiles larger
#                                    than 2GB on supported systems (tested
#                                    on Linux)
#      -D__USE_LONG_INT              Enables support for markers with up
#                                    to 64 alleles (default is 32). Tested
#                                    on systems where gcc supports the long
#                                    long data type and on Windows.
#      -D__ZLIB_AVAILABLE__          Enables support for GZIP'ed input files
# 
CFLAGS=-O2 -I./libsrc -I./merlin -I./pdf -I./clusters -D_FILE_OFFSET_BITS=64 -D__ZLIB_AVAILABLE__ -Wall

# executable file names and locations
BINDIR = executables
MERLIN = $(BINDIR)/merlin
MERLINX = $(BINDIR)/minx
MERLINREG = $(BINDIR)/merlin-regress
MERLINOFF = $(BINDIR)/merlin-offline
MERLINXOFF = $(BINDIR)/minx-offline
PEDSTATS = $(BINDIR)/pedstats
PEDWIPE = $(BINDIR)/pedwipe
PEDMERGE = $(BINDIR)/pedmerge
HAPMAPCONVERTER = $(BINDIR)/hapmapConverter
EXECUTABLES = $(MERLIN) $(MERLINX) $(MERLINREG) $(MERLINOFF) $(MERLINXOFF) \
              $(PEDSTATS) $(PEDWIPE) $(PEDMERGE) $(HAPMAPCONVERTER)

# MERLIN File Set
MERLINBASE = merlin/AssociationAnalysis merlin/FastAssociation \
 merlin/AnalysisTask merlin/Conquer \
 merlin/ConquerHaplotyping merlin/DiseaseModel \
 merlin/ParametricLikelihood merlin/GenotypeInference \
 merlin/Houdini \
 merlin/KongAndCox merlin/Manners merlin/MerlinBitSet \
 merlin/MerlinCluster merlin/MerlinCore \
 merlin/MerlinError merlin/MerlinFamily merlin/MerlinIBD \
 merlin/InformationContent merlin/MerlinCache merlin/MerlinModel \
 merlin/MerlinKinship merlin/MerlinKinship15 \
 merlin/MerlinHaplotype merlin/MerlinMatrix merlin/MerlinParameters \
 merlin/MerlinPDF merlin/MerlinSimulator merlin/MerlinSimwalk2 \
 merlin/MerlinSort merlin/NPL-ASP merlin/NPL-QTL \
 merlin/Magic merlin/Mantra merlin/Parametric merlin/QtlModel \
 merlin/Tree \
 merlin/TreeBasics merlin/TreeIndex merlin/TreeManager \
 merlin/TreeInfo merlin/TreeFlips merlin/VarianceComponents
MERLINHDR = $(MERLINBASE:=.h) merlin/TreeNode.h
MERLINSRC = $(MERLINBASE:=.cpp) merlin/Merlin.cpp
MERLINOBJ = $(MERLINSRC:.cpp=.o)
MERLINXOBJ = $(MERLINOBJ:.o=.X.o)

# Files for dealing with clustered markers
CLUSTERS = clusters/HaploFamily clusters/HaploGraph \
 clusters/HaploSet clusters/HaploTree clusters/Likelihood \
 clusters/SparseLikelihood clusters/Unknown
CLUSTERCPP = $(CLUSTERS:=.cpp)
CLUSTERHDR = $(CLUSTERS:=.h)
CLUSTEROBJ = $(CLUSTERS:=.o)
CLUSTERXOBJ = $(CLUSTERS:=.X.o)

# Regression Analysis
REGFILES = regress/AutoFit.cpp regress/Regress.cpp \
 regress/FancyRegression.cpp regress/RegressAnalysis.cpp \
 regress/RegressKinship.cpp regress/RegressParameters.cpp \
 regress/AutoFit.h \
 regress/FancyRegression.h regress/RegressAnalysis.h \
 regress/RegressKinship.h regress/RegressParameters.h 
REGBASE = merlin/Conquer merlin/DiseaseModel \
 merlin/Houdini merlin/Manners \
 merlin/MerlinBitSet merlin/MerlinCache merlin/MerlinCluster \
 merlin/MerlinCore merlin/MerlinError \
 merlin/ParametricLikelihood merlin/MerlinSimulator \
 merlin/MerlinSort merlin/Magic merlin/Mantra merlin/MerlinPDF \
 merlin/Parametric \
 merlin/Tree merlin/TreeBasics merlin/TreeInfo merlin/TreeFlips \
 merlin/TreeManager \
 regress/AutoFit \
 regress/FancyRegression regress/RegressAnalysis \
 regress/RegressKinship regress/RegressParameters
REGHDR = $(REGBASE:=.h) merlin/TreeNode.h
REGSRC = $(REGBASE:=.cpp) regress/Regress.cpp
REGOBJ = $(REGSRC:.cpp=.o)

# Offline Analysis File Set
OFFSRC = offline/Main.cpp $(MERLINBASE:=.cpp)
OFFOBJ = $(OFFSRC:.cpp=.o)
OFFXOBJ = $(OFFSRC:.cpp=.X.o) 

# Utility Library File Set
LIBFILE = libsrc/lib-goncalo.a
LIBMAIN = libsrc/BasicHash libsrc/Error libsrc/FortranFormat \
 libsrc/GenotypeLists libsrc/InputFile libsrc/IntArray libsrc/Hash \
 libsrc/LongArray libsrc/Kinship libsrc/KinshipX  libsrc/MapFunction \
 libsrc/MathCholesky libsrc/MathDeriv libsrc/MathFloatVector \
 libsrc/MathGenMin libsrc/MathGold libsrc/MathMatrix libsrc/MathStats \
 libsrc/MathNormal libsrc/MathSVD libsrc/MathVector \
 libsrc/MemoryInfo libsrc/MiniDeflate \
 libsrc/Parameters libsrc/Pedigree libsrc/PedigreeAlleleFreq \
 libsrc/PedigreeDescription libsrc/PedigreeFamily libsrc/PedigreeGlobals \
 libsrc/PedigreePerson libsrc/QuickIndex libsrc/Random libsrc/Sort \
 libsrc/StringArray libsrc/StringBasics libsrc/StringMap \
 libsrc/StringHash libsrc/TraitTransformations
LIBPED = libsrc/PedigreeLoader libsrc/PedigreeTwin libsrc/PedigreeTrim
LIBSRC = $(LIBMAIN:=.cpp) $(LIBPED:=.cpp)
LIBHDR = $(LIBMAIN:=.h) libsrc/Constant.h \
 libsrc/MathConstant.h libsrc/PedigreeAlleles.h libsrc/LongInt.h
LIBOBJ = $(LIBSRC:.cpp=.o)
 
# PDF Library File Sets
PDFLIB = pdf/libpdf.a
PDFFILES = pdf/PDF pdf/PDFfont pdf/PDFinfo pdf/PDFpage \
 pdf/PDFchartbasics pdf/PDFchartbar pdf/PDFlinechart \
 pdf/PDFhistogram \
 pdf/PDFchartaxis pdf/PDFchartlegend pdf/PDFchartmarker \
 pdf/PDFchartline pdf/PDFchartobject
PDFSRC = $(PDFFILES:=.cpp)
PDFHDR = $(PDFFILES:=.h)
PDFOBJ = $(PDFFILES:=.o)

# private parameters
FETCHDIR=$(HOME)/code
DISTRIBDIR=$(HOME)/code/distrib/merlin-$(VERSION)

# helpful screen listing available options
help : 
	@echo "MERLIN Source Distribution"
	@echo " "
	@echo "This Makefile will compile and install merlin on your system"
	@echo " "
	@echo "Type...           To..."
	@echo "make help         Display this help screen"
	@echo "make all          Compile merlin and related tools"
	@echo "make install      Install binaries in $(INSTALLDIR)"
	@echo "make install INSTALLDIR=directory_for_binaries"
	@echo "                  Install binaries in directory_for_binaries"
	@echo "make clean        Delete temporary files"

# make everything
all : $(EXECUTABLES)

$(EXECUTABLES) : $(BINDIR)

$(BINDIR) :
	mkdir -p $(BINDIR)

# dependencies for executables
$(MERLIN) : $(LIBFILE) $(PDFLIB) $(MERLINOBJ) $(CLUSTEROBJ)
	$(CXX) $(CFLAGS) -o $@ $(MERLINOBJ) $(CLUSTEROBJ) $(PDFLIB) $(LIBFILE) -lm -lz

$(MERLINX) : $(LIBFILE) $(PDFLIB) $(MERLINXOBJ) $(CLUSTERXOBJ)
	$(CXX) $(CFLAGS) -o $@ $(MERLINXOBJ) $(CLUSTERXOBJ) $(PDFLIB) $(LIBFILE) -lm -lz

$(MERLINREG) : $(LIBFILE) $(PDFLIB) $(REGOBJ) $(CLUSTEROBJ)
	$(CXX) $(CFLAGS) -o $@ $(REGOBJ) $(CLUSTEROBJ) $(PDFLIB) $(LIBFILE) -lm -lz

$(MERLINOFF) :  $(LIBFILE) $(PDFLIB) $(OFFOBJ) $(CLUSTEROBJ)
	 $(CXX) $(CFLAGS) -o $@ $(OFFOBJ) $(CLUSTEROBJ) $(PDFLIB) $(LIBFILE) -lm -lz

$(MERLINXOFF) :  $(LIBFILE) $(PDFLIB) $(OFFXOBJ) $(CLUSTEROBJ)
	 $(CXX) $(CFLAGS) -o $@ $(OFFXOBJ) $(CLUSTERXOBJ) $(PDFLIB) $(LIBFILE) -lm -lz

$(PEDSTATS) : pedstats-$(PSVERSION).tar.gz
	gunzip -c pedstats-$(PSVERSION).tar.gz | tar -xf - 
	cd pedstats-$(PSVERSION) ; $(MAKE) executables/pedstats
	cp pedstats-$(PSVERSION)/executables/pedstats executables
	rm -rf pedstats-$(PSVERSION)

$(PEDWIPE) : $(LIBFILE) extras/pedwipe.cpp 
	$(CXX) $(CFLAGS) -o $@ extras/pedwipe.cpp $(LIBFILE) -lm -lz

$(PEDMERGE) : $(LIBFILE) extras/pedmerge.cpp
	$(CXX) $(CFLAGS) -o $@ extras/pedmerge.cpp $(LIBFILE) -lm -lz

$(HAPMAPCONVERTER) : $(LIBFILE) extras/hapmapConverter.cpp
	$(CXX) $(CFLAGS) -o $@ extras/hapmapConverter.cpp $(LIBFILE) -lm -lz

$(LIBFILE) : $(LIBOBJ) $(LIBHDR)
	ar -cr $@ $(LIBOBJ)
	ranlib $@

$(PDFLIB) : $(PDFOBJ)
	ar -cr $@ $(PDFOBJ)
	ranlib $@

$(MERLINOBJ) : $(MERLINHDR) $(CLUSTERHDR) $(LIBHDR)

$(MERLINXOBJ) : $(MERLINHDR) $(CLUSTERHDR) $(LIBHDR)

$(CLUSTEROBJ) : $(CLUSTERHDR) $(MERLINHDR) $(LIBHDR)

$(CLUSTERXOBJ) : $(CLUSTERHDR) $(MERLINHDR) $(LIBHDR)

$(REGOBJ) : $(MERLINHDR) $(LIBHDR)

$(LIBOBJ) : $(LIBHDR)

$(PDFOBJ) : $(PDFHDR)

clean :
	-rm -f */*.a */*.o $(EXECUTABLES) 

install : all $(INSTALLDIR)
	@echo " "
	@echo Installing to directory $(INSTALLDIR)
	@echo To select a different directory, run
	@echo " "
	@echo make install INSTALLDIR=your_preferred_dir
	@echo " "
	cp $(EXECUTABLES) $(INSTALLDIR)

$(INSTALLDIR) :
	@echo " "
	@echo Creating directory $(INSTALLDIR)
	@echo " "
	@mkdir -p $(INSTALLDIR)

new-version :
	mkdir -p $(DISTRIBDIR) $(DISTRIBDIR)/extras 
	mkdir -p $(DISTRIBDIR)/merlin $(DISTRIBDIR)/regress
	mkdir -p $(DISTRIBDIR)/offline
	mkdir -p $(DISTRIBDIR)/libsrc $(DISTRIBDIR)/pdf
	mkdir -p $(DISTRIBDIR)/clusters
	cp ChangeLog LICENSE.twister README $(DISTRIBDIR)
	cp Makefile $(DISTRIBDIR)
	cp -R examples $(DISTRIBDIR)
	
fetch : 
	cd $(FETCHDIR) ; cp $(MERLINSRC) $(MERLINHDR) $(DISTRIBDIR)/merlin
	cd $(FETCHDIR) ; cp merlin-offline/Main.cpp $(DISTRIBDIR)/offline
	cd $(FETCHDIR) ; cp $(LIBSRC) $(LIBHDR) $(DISTRIBDIR)/libsrc
	cd $(FETCHDIR) ; cp $(REGFILES) $(DISTRIBDIR)/regress
	cd $(FETCHDIR) ; cp $(PDFSRC) $(PDFHDR) $(DISTRIBDIR)/pdf
	cd $(FETCHDIR) ; cp $(CLUSTERHDR) $(CLUSTERCPP) $(DISTRIBDIR)/clusters
	cd $(DISTRIBDIR) ; wget -r -nd -N http://www.sph.umich.edu/csg/abecasis/pedstats/download/pedstats-$(PSVERSION).tar.gz 
	cp $(FETCHDIR)/pedwipe/pedwipe.cpp $(DISTRIBDIR)/extras
	cp $(FETCHDIR)/pedmerge/pedmerge.cpp $(DISTRIBDIR)/extras
	cp $(FETCHDIR)/hapmapConverter/hapmapConverter.cpp $(DISTRIBDIR)/extras
	cd $(DISTRIBDIR); csh ../stamp MERLIN

.c.o :
	$(CXX) $(CFLAGS) -o $@ -c $*.c

.cpp.X.o : 
	$(CXX) $(CFLAGS) -o $@ -c $*.cpp -DVERSION=\"$(VERSION)\" -D__CHROMOSOME_X__

.cpp.o : 
	$(CXX) $(CFLAGS) -o $@ -c $*.cpp -DVERSION=\"$(VERSION)\"

archive : clean
	mkdir -p merlin-$(VERSION)
	cp -R LICENSE.twister README Makefile ChangeLog extras merlin-$(VERSION)
	cp -R clusters merlin regress offline libsrc pdf examples merlin-$(VERSION)
	cp -R pedstats-$(PSVERSION).tar.gz merlin-$(VERSION)
	tar -cvf merlin-$(VERSION).tar merlin-$(VERSION) 
	gzip -f --best merlin-$(VERSION).tar
	rm -rf merlin-$(VERSION)

distrib : $(EXECUTABLES)
	mkdir -p merlin-$(VERSION)
	cp -R LICENSE.twister README ChangeLog $(EXECUTABLES) examples merlin-$(VERSION)
	tar -cvf `uname`-merlin.tar merlin-$(VERSION)
	gzip -f `uname`-merlin.tar
	rm -rf merlin-$(VERSION)

windowszip : $(EXECUTABLES)
	mkdir -p merlin-$(VERSION)
	cp -R LICENSE.twister README ChangeLog $(EXECUTABLES) examples merlin-$(VERSION)
	echo '@ECHO OFF' > merlin-$(VERSION)/STARTHERE.bat
	echo 'CMD /K "SET PATH=%CD%;%PATH%"' > merlin-$(VERSION)/STARTHERE.bat
	zip -r Windows-merlin.zip merlin-$(VERSION)
	rm -rf merlin-$(VERSION)


.SUFFIXES : .cpp .c .o .X.o $(SUFFIXES)

