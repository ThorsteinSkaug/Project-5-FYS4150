# Project-5-FYS4150

This is a repository for project 5 in FYS4150 (Computational Physics) at the University of Oslo. The goal of this project is to simulate the two-dimensional time-dependent Schrödinger equation. We are looking at one particle being shot towards a wall with either one, two, or three slits. One can also simulate without a wall by swithicng the wall "off". 

To run the simulation you need to run the ```main.cpp``` file, give the name of the input file and the name you want the output file to have. Each row of the output file will be the state of the system at that given timestep, all even entries being the real part and every odd entry being the imaginary part (starting from 0). The second to last enry is the time for that row, and the lastentry is how much the probability differ from 1.

The repository has four parts. One folder ```C++ files```, containing all the necessary code for running the simulation. One folder for all python plot files, containing the ```plot.py``` file with all "normal" plot codes and some calculations, and ```animation.py``` which is the code for running the animation. One folder containing all input files needed for this assignment. And finaly one folder containing all the plots.

## C++ code

The C++ source folder contains four files. A ```main.cpp``` for doing the simulation and writing the results to a .txt file. A header file, ```Schrodinger.hpp``` showing and explaining all functions used for this assignment. A ```Schrodinger.cpp``` file contianing all the functions definitions for the functions in the header file. And finaly, a ```makefile``` which makes it easy to run the code.

To run the code just type ```run All``` inn the terminal and give the input file and output file name.

### Input files

For the code to run properly you need to give it an input file on the right form. The input file need to be on the same form as the examples given in the folder ```Input files```. To simulate with other variables just change the numbers. The entry "wall" takes string arguments, either "on" or "off", "nr_slits" takes an interger value between 1 and 3 (If wall is "off" then it does not matter what this value is). All other entries takes doubles.

## Plotting

