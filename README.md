# PhylogeneticTrees.jl

PhylogeneticTrees is a package implementing miqoGraph, a fast algorithm to fit admixture trees using mixed-integer quadratic optimization.  It was developed using the [Julia](http://julialang.org/) language and the [Gurobi](http://www.gurobi.com/) optimization solver.

----
## Installation

PhylogeneticTrees is implemented using the Gurobi optimization solver, which has a [free academic license](http://www.gurobi.com/registration/academic-license-reg).  The instructions for installation can be found [here](http://www.gurobi.com/documentation/).

The Julia language can be downloaded [here](https://julialang.org/downloads/).

PhylogeneticTrees can be installed through the Julia package manager:

```
(v1.0) pkg> add "https://github.com/juliayyan/PhylogeneticTrees.jl" 
```

To check that the installation is working properly, run the following within the package manager:

```
(v1.0) pkg> test PhylogeneticTrees
```

**NOTE**.  The dependencies for PhylogeneticTrees include JuMP (v0.18.5), Gurobi, CSV (v0.4.3), DataFrames (v0.17.1), and LightGraphs.  Newer or older versions of these packages may cause errors.  To install a specific version of these packages, you can use the following commands:

```
(v1.0) pkg> add JuMP@0.18.5
```

```
(v1.0) pkg> add CSV@0.4.3
```

```
(v1.0) pkg> add DataFrames@0.17.1
```

----
## Usage
Example code and datasets can be found in the test/ directory.

----
## Citing PhylogeneticTrees
