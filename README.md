# multicore-n-planet-simulator - Euler Method and Runge-Kutta 2
Submitted as part of the degree of Msci Natural Sciences (3rd year) to the Board of Examiners in the Department of Computer Sciences, Durham University. 
This summative assignment was assessed and marked by the professor of the module in question:
## Grade: 1st - 90/100, 1st in year (of 78 students).
Vectorized and multicore n-body simulators written and extensively optimised in C++ for scalability to millions of particles/planets to be run on a single node of a supercomputer.
## Demo video (taken from my [portfolio page](https://github.com/shadowbourne)):
> ![Gifdemo1](https://user-images.githubusercontent.com/18665030/136667074-f1f25ed1-1a71-44f2-b94c-5503c86e52ec.gif)
> 
> Vectorized and multicore n-body simulator(s) written and extensively optimised in C++ for scalability to millions of particles/planets to be run on a single node of a supercomputer.

## Contents:
* step-1.cpp contains the basic non-vectorized code for the n-body solver using the first-order Euler Method numerical time-stepping scheme (18/20).
* step-2.cpp contains a vectorized version of step 1 for more efficient computation of for loops (20/20).
* step-3.cpp upgrades the numerical scheme from steps 1 & 2 to a second order RK2 (Runge-Kutta 2) scheme. Step 3 is also then vectorized (16/20).
* step-4.cpp takes Step 3 and parallelizes the code using OpenMP to run for-loops on multiple cores efficiently (16/20).
* report.pdf contains an indepth analysis of the Strong Scalability and Convergence of the code and algorithms produced (20/20).

## Collision Rules:
![Rules](rules.png?raw=true "Rules")

## Analysis of Strong Scalability and Convergence of the code and algorithms produced (from report.pdf):
![Report](report.png?raw=true "Report")
