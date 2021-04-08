# ROADEF/EURO Challenge 2020: Maintenance Planning Problem

> André L. Maravilha<sup>1, 2</sup>  
> <sup>1</sup> *Dept. of Informatics, Management and Design - Centro Fed. de Edu. Tecnológica de Minas Gerais ([url](https://www.cefetmg.br/))*  
> <sup>2</sup> *Operations Research and Complex Systems Lab. - Universidade Federal de Minas Gerais ([url](https://orcslab.github.io/))*

This repository keeps the source code (and submissions) of my proposed strategy for solving the Maintenance Planning Problem.

Interested people can find more details about the problem and the competition on the event's official page through this [url](https://www.roadef.org/challenge/2020/en/index.php).


## 1. Proposed strategy

The proposed strategy is a heuristic based on Benders' Decomposition technique [1], in which the master problem defines the starting time of each intervention while ensuring the constraints below are met:

- Non-preemptive scheduling;
- Interventions are scheduled once;
- No workflow left;
- Resource constraints;
- Disjunctive constraints.

The subproblems are responsible for determining the expected excess of the planning. The constraints in the master problem guarantee the existence of feasible solutions to the subproblem. Then, in each iteration, new cuts guide the search for improved solutions.

Besides the advantage of producing only cuts that guide the search through improved solutions, the proposed strategy generates cuts without solving the dual version of the subproblems through a solver. The cuts are generated only by inspecting the solution's values obtained in the master problem and the parameters of the subproblem, which speeds up the search for improved solutions.

However, the advantages obtained with the cuts proposed in this strategy bring some drawbacks. The cuts are not valid for the original problem since they may remove valid solutions (including the optimal one). Then, the final solution returned by the proposed method is not guaranteed to be the optimal one. Also, the lower bounds are not valid to the problem too.


## 2. How to build the project

#### 2.1. Important comments before building the project

To compile this project, you need to have Gurobi (version 9.1.+), CMake (version 3.13 or later), and a compatible compiler installed on your computer. The code has not been tested on versions earlier than the ones specified.

#### 2.2. Building the project

Inside the root directory of the project, run the following commands:
```
mkdir build
cd build
cmake -DCMAKE_BUILD_TYPE=Release ..
make
```


## References

1. Benders, J.F. "Partitioning procedure for solving mixed-variables programming problems". Numerische Mathematik, 4, 238-252, 1962.
