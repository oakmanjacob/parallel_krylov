# Parallel Krylov Method Analysis
By Jacob Oakman  
For CSE 398: Iterative Methods for Large Sparse Linear Systems

## Overview
This represents a semester long project for a seminar at Lehigh University. The goal of this project is to implement and examine both sequential and parallel versions of several Krylov methods, starting with GMRES. These will be compared for use in solving generated Convection-Diffusion problems.

### Notes
After running gperf, the lines which have the most impact on performace is one that is combined vector addition and vector-scalar product an one which computes the dot product of two vectors. These may be somewhat parallelizable but first it will be worth looking at using AVX SIMD instructions to optimize these.

### Implementation Plan
- [x] Set up a repository to store the project files and create a docker image which can ensure a stable build environment.
- [x] Implement infrastructure for the project including writing classes for storing sparse matrices using the Compressed Sparse Row and Compressed Diagonal formats. Part of this project will be examining how these and other datastructures can affect performance.
- [x] Test that the infrastructure is working correctly using dummy data.
- [x] Implement a sequential version of GMRES based on the existing code from Dr. Saad. This may require adding and testing additional functionality for the sparse matrix implementations.
- [x] Implement logging using spdlog library
- [x] Implement a generator for basic diagonal and tridiagonal matrixes to make sure gmres is working.
- [x] Integrate a generator for toy Convection-Diffusion systems and test that they work correctly.
- [x] Set up the ability to import matrix systems via csv.
- [x] Implement a matvec operation for a matrix stored in the sparse CSR format and a full vector
- [x] Implement a verbatum version of mgmres to ascertain how much the optimizations helped
- [x] Split vector functions into a seperate utility file
- [x] Implement GMRES versions which don't use complex numbers
- [ ] Complete the function level comments to better explain what is going on in the code
- [x] Use AVX instructions to optimize functions for non-complex implementations
- [ ] Generate a sequence of Convection-Diffusion systems to use for testing
- [ ] Test and plot the performance of the GMRES implementations -> Iterations & Execution Time vs Size of System
- [ ] 
- [ ] Implement a parallel version of GMRES using optimization found in research. This may require creating a new data structure for storing the values in a way which prevents data races.


### Stretch Goals
- [ ] Compare the speed of these implementations to the default MATLAB implementation.
- [ ] Test scaling on a machine with a large number of cpu cores.
- [ ] Implement a basic ILU or Jacobi preconditioner to improve convergence and evaluate how this preconditioner affects data 
- [ ] Look for novel optimizations to the parallel implementation of GMRES and attempt improvements using other novel data structures.
- [ ] Implement and test sequential and parallel versions of other iterative methods like CG, MINRES or even non-Krylov methods such as Guass-Siedel and Jacobi.

### Deliverables
- [x] Project Proposal
- [x] Checkpoint
- [ ] Final Paper
- [ ] This codebase

### Ideas
0. What if we tried to use restarting to refine our guess with multiple instances working in parallel and then being combined?
1. What if we implemented a parallel version of the matvec?


## Setup
This section outlines everything whch needs to be done in order to build and run the codebase. I've provided a docker container in order to avoid issues with configuration. More about how to install docker can be found [here](https://docs.docker.com/get-docker/).  
  
This code can also be built and run on ubuntu or in the ubuntu version of WSL.

### Docker Build Image
```bash
$> docker build -t jco222krylov .
```

### Run Docker Container
```bash
$> docker run --privileged --rm -v $(pwd):/root --name jco222krylovenv -it jco222krylov
```

### Build the codebase
The make command will automatically create a build directory and move all built files to it
```bash
$> make
```

### Run
The -h flag will produce a help which explains the other flags.
```bash
$> ./build/krylov -h
```

## Testing

There are currently no test yet run

## Results

There are currently no results to display

## References
### [1] Eric G. Parker and James F. O'Brien. "Real-Time Deformation and Fracture in a Game Environment". In Proceedings of the ACM SIGGRAPH/Eurographics Symposium on Computer Animation, pages 156–166, August 2009.

This paper goes in depth about how a discretization and linear solve was used to efficiently simulate deformation and fractures in a video game environment. It serves as an example of a real world application of these methods. They touch on multiple different problems which needed to be solved 

### [2] E.I. Ioannidis, N. Cheimarios, A.N. Spyropoulos, A.G. Boudouvis. “On the performance of various parallel GMRES implementations on CPU and GPU clusters”. arXiv 1906.04051 (2019)

This paper looks at several different implementations of GMRES on CPU and GPUs. I specifically want to look at the CPU implementations as a jumping off point for my project. A straightforward version of this project could be to replicate the testing done in this paper and verify some of the conclusions.

### [3] David’s Palitta and Valeria Simoncini. “Matrix-Equation-Based Strategies for Convection-Diffusion Equations” arXiv https://arxiv.org/abs/1501.02920 (2015)

This paper looks at how Convection-Diffusion equations can be solved using linear equations. It should be a good resource in understanding how this project relates to its application and may provide some insight into the underlying properties of the matrices involved in these problems.
