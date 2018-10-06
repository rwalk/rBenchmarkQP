# rBenchmarkQP
Benchmark Problems and Analysis for Quadratic Programming in R

This is a collection of simple scripts for profiling several quadratic programming libraries in R. Details are available in this [blog post](https://rwalk.xyz/sparse-quadratic-programming-with-osqp/).

## Prerequisites
To run the profiling, you'll need to install the quadratic programming libraries `quadprog`, `kernlab`, `ipoptr` and `osqp`.  With the exception of `ipoptr`, these packages can be installed with `install.package`. Details on install `ipoptr` can be found [here](https://rwalk.xyz/sparse-quadratic-programming-with-ipoptr/).

## Running
You can run the scripts from the command line with:
```
git clone https://github.com/rwalk/rBenchmarkQP
cd rBenchmarkQP
```
Then from within `R`, execute:
```
source("profile.R")
```
