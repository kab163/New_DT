SRC=full.C 
SRC2=HRC_DT.C
INCLUDE=/usr/local/include/cleap-0.3.2
LIB=/usr/local/lib

default: dt

dt: $(SRC)
	g++ -std=c++11 $(SRC) -o dt -L$(LIB) -I$(INCLUDE) -lcleap

hrc: $(SRC2)
	g++ -std=c++11 $(SRC2) -o hrc -L$(LIB) -I$(INCLUDE) -lcleap

clean:
	rm dt outputMesh.off
