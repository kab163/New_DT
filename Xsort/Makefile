SRC3=GenMeshDT.C
SRC4=GenDelBound.C
C=g++ -std=c++11
INCLUDE=/usr/local/include/cleap-0.3.2
LIB=/usr/local/lib

default: gmdt gdelb

gmdt: $(SRC3)
	$(C) $(SRC3) -o gmdt -L$(LIB) -I$(INCLUDE) -lcleap

gdelb: $(SRC4)
	$(C) $(SRC4) -o gdelb -L$(LIB) -I$(INCLUDE) -lcleap

clean:
	rm gmdt gdelb outputMesh.off
