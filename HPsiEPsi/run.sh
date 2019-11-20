
mingw32-g++ -Wall -I/usr/local/include -c main.cpp
mingw32-g++ -L/usr/local/lib -o main.exe main.o -llibgsl -llibgslcblas -lm
main.exe
