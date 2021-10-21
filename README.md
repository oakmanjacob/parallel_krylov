# Parallel Krylov Method Analysis
By Jacob Oakman  
For CSE 398: Iterative Methods for Large Sparse Linear Systems

## Overview
This represents a semester long project for a seminar at Lehigh University. The goal of this project is to implement and examine both sequential and parallel versions of several Krylov methods, starting with GMRES. These will be compared for use in solving generated Convection-Diffusion problems.

### Implementation Plan
0. Set up a repository to store the project files and create a docker image which can ensure a stable build environment.
1. Implement infrastructure for the project including writing classes for storing sparse matrices using the Compressed Sparse Row and Compressed Diagonal formats. Part of this project will be examining how these and other datastructures can affect performance.
2. Test that the infrastructure is working correctly using dummy data.
3. Implement a generator for toy Convection-Diffusion systems and test that they work correctly.
3. Implement a sequential version of GMRES based on the existing code from Dr. Saad. This may require adding and testing additional functionality for the sparse matrix implementations.
4. Test and plot the performance of the GMRES implementation
5. Implement a parallel version of GMRES using optimization found in research. This may require creating a new data structure for storing the values in a way which prevents data races.


### Stretch Goals
0. Compare the speed of these implementations to the default MATLAB implementation.
1. Test scaling on a machine with a large number of cpu cores.
2. Look for novel optimizations to the parallel implementation of GMRES and attempt improvements using other novel data structures.
3. Implement and test sequential and parallel versions of other iterative methods like CG, MINRES or even non-Krylov methods such as Guass-Siedel and Jacobi.

### Deliverables
0. Project Proposal
1. Checkpoint
3. Final Paper
4. This codebase

## Setup
This section outlines everything whch needs to be done in order to build and run the codebase. I've provided a docker container in order to avoid issues with configuration. More about how to install docker can be found [here](https://docs.docker.com/get-docker/).  
  
This code can also be built and run on ubuntu or in the ubuntu version of WSL.

### Docker Build Image
```bash
$> docker build -t jco222krylov .
```

### Run Docker Container
```bash
$> docker run --privileged --rm -v ~/path/to/repo:/root --name jco222krylovenv -it jco222krylov
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
