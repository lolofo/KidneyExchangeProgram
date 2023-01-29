- [Kidney exchange program](#kidney-exchange-program)
  * [Environment](#environment)
  * [Organization](#organization)

# Kidney exchange program

Kidney exchange is nowaday a big issue. Several people suffer from kidney problem and need transplants.
We propose here an approach to this problem through the framework of stochastic optimization.
This repository has been made within our "Optimization under uncertainty" course given by Jeremy Omer at INSA Rennes.

We chose to treat the problem as following : our goal is to make groups of people (cluster) which are linked together. Within these groups the crossed tests will be proceeded and then after these tests we will create cycles to proceed the real exchange.

## Environment

To execute the different files from this repository you will need several packages.
If you execute the following command in your terminal, at the root of this rep :

```{command line}
julia requirements.jl
```

All the packages will be installed, and then you will be able te execute all the scripts.


## Organization


This repository contains multiple folders :

- `src\sotchastic_framework` : this folder contains all the methods implemented to solve the stochastic problem
- `src\deterministic_framework` : this folder contains some implementation for the deterministic version of the problem.
- `notebook` : this folder contains the notebook for our final evaluation
- `utils` : this folder contains all the function for the secondary tasks such as read and process the data, plot the graphs, ...