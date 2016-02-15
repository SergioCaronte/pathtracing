#!/bin/bash 
g++ -O3 -std=c++11 -fopenmp -lstdc++  src/math/*.h src/*.h src/*.cpp  -o bin/Pathtracing
