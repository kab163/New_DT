SRC=v4X.C 
SRC2=v4Y.C
SRC3=HRC_DT.C
INCLUDE=/usr/local/include/cleap-0.3.2
LIB=/usr/local/lib

default: v4y v4x

v4y: $(SRC2)
	g++ -std=c++11 $(SRC2) -o v4y -L$(LIB) -I$(INCLUDE) -lcleap

v4x: $(SRC)
	g++ -std=c++11 $(SRC) -o v4x -L$(LIB) -I$(INCLUDE) -lcleap

clean:
	rm v4y v4x outputMesh.off
