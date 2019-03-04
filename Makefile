SRC=DTinput.C 
SRC2=v2DTinput.C
SRC3=HRC_DT.C
SRCcu=v2DTi.cu
SRCtest=cleap_test.c
INSTALL=/usr/local

default: v3

prog: $(SRC)
	g++ -std=c++11 $(SRC) -o dt

v2: $(SRC2)
	g++ -std=c++11 $(SRC2) -o v2 -L$(INSTALL)/lib -I$(INSTALL)/include/cleap-0.3.2 -lcleap

v3: $(SRC3)
	g++ -std=c++11 $(SRC3) -o v3 -L$(INSTALL)/lib -I$(INSTALL)/include/cleap-0.3.2 -lcleap

app: $(SRCtest)
	g++ -std=c++11 $(SRCtest) -o app -L$(INSTALL)/lib -I$(INSTALL)/include/cleap-0.3.2 -lcleap

gpu: $(SRCcu)
	nvcc -arch=sm_35 $(SRCcu) -o gpu -L$(INSTALL)/lib -I$(INSTALL)/include/cleap-0.3.2 -lcleap

clean:
	rm v3
