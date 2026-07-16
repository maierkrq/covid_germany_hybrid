include ../../Makefile.Local
include ../../Makefile.Rules

%.o: %.cc
	$(CC) $(FLAGS) $(INCLUDE) $(CLANGINC) $< -c -o $@
	
opensuse_lapacke_fix: opensuse_lapacke_fix.o
	cp -p opensuse_lapacke_fix.o $(KASKADE7)/libs

subsystem_lapacke_fix: subsystem_lapacke_fix.o
	cp -p subsystem_lapacke_fix.o $(KASKADE7)/libs
