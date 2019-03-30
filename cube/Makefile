SRC=v1cube.C 
INCLUDE=/usr/local/include/cleap-0.3.2
LIB=/usr/local/lib

default: cube

cube: $(SRC)
	g++ -std=c++11 $(SRC) -o cube -L$(LIB) -I$(INCLUDE) -lcleap

clean:
	rm cube outputMesh.off
