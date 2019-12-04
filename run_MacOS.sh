uncrustify -c /usr/local/share/uncrustify/objc.cfg -f main.cpp -o main.cpp
if g++ -Wall -I/usr/local/include -c main.cpp;
g++ -L/usr/local/lib -o main.exe main.o -lgsl -lgslcblas -lm; then
./main.exe
gnuplot gnuPlotScript
fi
