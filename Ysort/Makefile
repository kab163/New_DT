SRC4=slice.C
LIB=/usr/local/lib
INCLUDE=/usr/local//include/cleap-0.3.2
C=g++ -std=c++11

default: s1

s1: $(SRC4)
	$C $(SRC4) -o s1 -L$(LIB) -I$(INCLUDE) -lcleap

clean:
	rm s1 outputMesh.off
