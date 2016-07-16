app = NHprog
flags = -O2 -std=c++0x -fopenmp
compiler = icc
objects = main.o particle.o NHthermostat.o
headers = Bridge.h Constants.h Libraries.h

$(app): $(objects) $(headers)
	$(compiler) -o $(app) $(flags) *.o

main.o: main.cpp $(headers)
	$(compiler) -c $(flags) main.cpp

particle.o: particle.cpp $(headers)
	$(compiler) -c $(flags) particle.cpp

NHthermostat.o: NHthermostat.cpp $(headers)
	$(compiler) -c $(flags) NHthermostat.cpp

clean:
	@rm -rf *.o $(app)
