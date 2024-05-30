test: structures.o test.o VMesh.o distribution.o CollisionNodes.o
	g++ structures.o test.o VMesh.o distribution.o CollisionNodes.o -o test
structures.o: structures.cpp
	g++ -c structures.cpp
test.o: test.cpp
	g++ -c test.cpp
VMesh.o: VMesh.cpp
	g++ -c VMesh.cpp
distribution.o: distribution.cpp
	g++ -c distribution.cpp
CollisionNodes.o: CollisionNodes.cpp
	g++ -c CollisionNodes.cpp
clean:
	rm *.o test