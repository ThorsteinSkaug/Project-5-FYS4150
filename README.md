# Project-5-FYS4150

This is a repository for project 5 in FYS4150 (Computational Physics) at the University of Oslo. The goal of this project is to simulate the two-dimensional time-dependent Schrödinger equation. We are looking at one particle being shot towards a wall with either one, two, or three slits. One can also simulate without a wall by swithicng the wall "off". 

To run the simulation you need to run the ```main.cpp``` file, give the name of the input file and the name you want the output file to have. Each row of the output file will be the state of the system at that given timestep, all even entries being the real part and every odd entry being the imaginary part (starting from 0). The second to last enry is the time for that row, and the lastentry is how much the probability differ from 1.

The repository has four parts. One folder ```C++ files```, containing all the necessary code for running the simulation. One folder for all python plot files, containing the ```plot.py``` file with all "normal" plot codes and some calculations, and ```animation.py``` which is the code for running the animation. One folder containing all input files needed for this assignment. And finaly one folder containing all the plots.

## C++ code

The C++ source folder contains two files, one ```makefile``` and one file containing all necessary C++ code called ```project4_main.cpp```. To run the C++ code see the next subsection:


## Plotting
To plot the results from ```project4_main.cpp``` simply run the ```plot.py``` file. 
