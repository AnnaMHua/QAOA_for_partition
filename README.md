# QAOA solver for number partitioning problem
Author: Meng Hua

## Introduction
This project implements the Quantum Approximate Optimization Algorithm(QAOA) using ```qiskit``` and solves the [number partitioning problem](https://en.wikipedia.org/wiki/Partition_problem).

## Number Partitioning
Number partitioning - also referred to as the number bi-partitioning problem, or the 2-partition problem - is one of Karp’s original NP-hard problems, and can be stated as follows. We have a set S of n numbers, and the goal of the optimisation problem is separate them into two disjoint subsets S1 and S2 such that the difference of the sums of the two subsets is minimised.

For example, if ```S = {4, 5, 6, 7, 8}```, then one possible partition would be ```S_1 = {4, 8}``` and ```S2 = {5, 6, 7}```. 
The sums of these two subsets are 12 and 18 respectively, giving a difference of 6. The optimal solution is ```S1 = {4, 5, 6}```, ```S2 = {7, 8}```, since the sum of both of these subsets is 15, and hence the difference is zero.

## Code structure
* ```QAOA.py```: contains a class for solving number partitioning problem. Includes utilities for constructing QAOA quantum circuit and  classical optimizer.
* ```NumberPartition.ipynb```: example and tutorial notebook with result analysis.
* ```sample.partition```: input file of a set of numbers. Used by ```NumberPartition.ipynb```

## Package requirement
* ```qiskit```
* ```scipy```
* ```numpy```
* ```matplotlib```



