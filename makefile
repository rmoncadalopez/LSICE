#HOST=$(shell hostname -d)

# Modified from Harmon desktop
#CC = g++-9 -I/usr/local/include -I/usr/local/include/eigen3 -L/usr/local/lib -lmpi  Original
#CC_INCLUDE = -I/usr/include/eigen3 -I/usr/include/openmpi -I/usr/local/include Original

CC = g++-10 -I/usr/local/Cellar/open-mpi/4.0.5/include -L/usr/local/opt/libevent/lib -L/usr/local/Cellar/open-mpi/4.0.5/lib -lmpi
#CC = mpicxx
#CC = -std=c++11
#CC = /usr/local/Cellar/gcc/10.2.0/bin/g++-10;

CC_INCLUDE = -I/usr/include/Eigen -I/usr/include/openmpi -I/usr/local/include -I/usr/local/include/openmpi -I/usr/local/Cellar/libopenmp  -L/usr/local/lsmlib/lib/ -L/usr/local/lsmlib/config -L/usr/local/lsmlib/src/toolbox -I/usr/local/Cellar/open-mpi -I/usr/local/Cellar/open-mpi/4.0.5/include -L/usr/local/opt/libevent/lib -L/usr/local/Cellar/open-mpi/4.0.5/lib -lmpi -L/usr/local/lsmlib
LIBS = -L/usr/local/lsmlib/lib -llsm_serial -llsm_toolbox

TARGETS = Main Main_Periodic Main_Fluid Main_Periodic_Ocean



all: $(TARGETS)

Main: Main.cpp $(wildcard *.h)
	$(CC) $< -O3 -o $@ -O3 $(CC_INCLUDE) -fopenmp -Wall -Wextra -pedantic -Wno-long-long -Wunused $(LIBS) -Wno-int-in-bool-context -std=c++11 -std=c++17


Main_Periodic: Main_Periodic.cpp $(wildcard *.h)
	$(CC) $< -O3 -o $@ -O3 $(CC_INCLUDE) -fopenmp -Wall -Wextra -pedantic -Wno-long-long -Wunused $(LIBS) -Wno-int-in-bool-context -std=c++11 -std=c++17

Main_Fluid: Main_Fluid.cpp $(wildcard *.h)
	$(CC) $< -O3 -o $@ -O3 $(CC_INCLUDE) -fopenmp -Wall -Wextra -pedantic -Wno-long-long -Wunused $(LIBS) -Wno-int-in-bool-context -std=c++11 -std=c++17

Main_Periodic_Ocean: Main_Periodic_Ocean.cpp $(wildcard *.h)
	$(CC) $< -O3 -o $@ -O3 $(CC_INCLUDE) -fopenmp -Wall -Wextra -pedantic -Wno-long-long -Wunused $(LIBS) -Wno-int-in-bool-context -std=c++11 -std=c++17			

clean:
	rm -f $(TARGETS)

again: 
	clean $(TARGETS)