g++ -Wall -I/usr/local/include -c main.cpp
g++ -L/usr/local/lib -o main.exe main.o -lgsl -lgslcblas -lm
./main.exe
