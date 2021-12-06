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
- [x] Implement a matvec operation for a matrix stored in the sparse CSR format and a full vector.
- [x] Implement a verbatum version of mgmres to ascertain how much the optimizations helped.
- [x] Split vector functions into a seperate utility file.
- [x] Implement GMRES versions which don't use complex numbers.
- [x] Complete the function level comments to better explain what is going on in the code.
- [x] Use AVX instructions to optimize functions for non-complex implementations.
- [x] Change so result information is always printed and there is a more clean way to print the values.
- [x] Generate a sequence of Convection-Diffusion systems to use for testing.
- [x] Test and plot the performance of the GMRES implementations -> Iterations & Execution Time vs Size of System.
- [x] Implement a restarted GMRES method which runs multiple versions in parallel and uses them to coordinate better initial guesses using a linear combination of the results.
- [x] Test the effects on iterations and computation time vs the normal optimized restared method.
- [x] Plot different versions of restarted GMRES vs parallel restarted GMRES


### Stretch Goals
- [x] Compare the speed of these implementations to the default MATLAB implementation.
- [x] Generate plots comparing MATLAB to this implementation
- [ ] Test scaling on a machine with a large number of cpu cores.
- [ ] Implement the ability for GMRES to handle preconditioned systems.
- [ ] Implement a parallel version of GMRES using optimization found in research. This may require creating a new data structure for storing the values in a way which prevents data races.
- [ ] Implement a basic ILU or Jacobi preconditioner to improve convergence and evaluate how this preconditioner affects data.
- [ ] Look for novel optimizations to the parallel implementation of GMRES and attempt improvements using other novel data structures.
- [ ] Implement and test sequential and parallel versions of other iterative methods like CG, MINRES or even non-Krylov methods such as Guass-Siedel and Jacobi.

### Additional TODOs
- [x] Set up a testing script so we don't need to run individual tests from the command line.
- [x] Fix the makefile so we can use make debug to compile with -O0 and debug flags, make gprof to compile with -pg and make to compile for production.
- [ ] Expand this README with more information about the testing and project.

### Deliverables
- [x] Project Proposal
- [x] Checkpoint
- [x] Final Paper
- [x] This codebase

### Ideas
0. What if we implemented a parallel version of the matvec?  
    Answer: We don't actually spend enough time inside the matvec funtion to warrent threaded parallelism.
1. What if we tried to use restarting to refine our guess with multiple instances working in parallel and then being combined?


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

Here is an example of how to run the tests for 4096 x 4096 systems 1..10
```bash
$> bash ./test/run_tests.sh 4096
```

The following scripts are provided for builk testing
- ./test/run_tests_parallel.sh
- ./test/run_tests_sequential.sh
- ./test/run_tests_weird.sh

There is also python script provided to extract an array of just the execution times from the files generated from these tests.
```bash
$> python3 ./test/extract_times.py {RESULT_FILE_NAME_HERE}.txt
```

## Results

### Sequential Optimizations
![Execution Times of many 1024 x 1024 C-D Systems](https://github.com/oakmanjacob/parallel_krylov/blob/main/results/n1024/plot_1024.png)
![Execution Times of many 4096 x 4096 C-D Systems](https://github.com/oakmanjacob/parallel_krylov/blob/main/results/n4096/plot_4096.png)
![Execution Times of many 16384 x 16384 C-D Systems](https://github.com/oakmanjacob/parallel_krylov/blob/main/results/n16384/plot_16384.png)

### Parallel Implementation
![Graph Showing Bad Scaling on Parallel Implementation](https://github.com/oakmanjacob/parallel_krylov/blob/main/results/n16384/plot_parallel_r75_p10.png)

### Weird Restarted GMRES
![Weird and SIMD GMRES on 16384 x 16384 C-D Iteration 4 at Varying Restart Values](https://github.com/oakmanjacob/parallel_krylov/blob/main/results/n16384/plot_wgmres_k4.png)
![Weird and SIMD GMRES on 16384 x 16384 C-D Iteration 10 at Varying Restart Values](https://github.com/oakmanjacob/parallel_krylov/blob/main/results/n16384/plot_wgmres_k10.png)


## References
### [1] Eric G. Parker and James F. O'Brien. "Real-Time Deformation and Fracture in a Game Environment". In Proceedings of the ACM SIGGRAPH/Eurographics Symposium on Computer Animation, pages 156–166, August 2009.

This paper goes in depth about how a discretization and linear solve was used to efficiently simulate deformation and fractures in a video game environment. It serves as an example of a real world application of these methods. They touch on multiple different problems which needed to be solved 

### [2] E.I. Ioannidis, N. Cheimarios, A.N. Spyropoulos, A.G. Boudouvis. “On the performance of various parallel GMRES implementations on CPU and GPU clusters”. arXiv 1906.04051 (2019)

This paper looks at several different implementations of GMRES on CPU and GPUs. I specifically want to look at the CPU implementations as a jumping off point for my project. A straightforward version of this project could be to replicate the testing done in this paper and verify some of the conclusions.

### [3] David’s Palitta and Valeria Simoncini. “Matrix-Equation-Based Strategies for Convection-Diffusion Equations” arXiv https://arxiv.org/abs/1501.02920 (2015)

This paper looks at how Convection-Diffusion equations can be solved using linear equations. It should be a good resource in understanding how this project relates to its application and may provide some insight into the underlying properties of the matrices involved in these problems.

### [4] Y. Saad and M. H. Schultz, "GMRES: A generalized minimal residual algorithm for solving nonsymmetric linear systems", SIAM J. Sci. Stat. Comput., vol. 7, no. 3, pp. 856-869, 1986.

### [5] Arielle Carr, “Recycling Techniques for Sequences of Linear Systems and Eigenproblems”, Dissertation submitted to the Faculty of the Virginia Polytechnic Institute and State University, http://hdl.handle.net/10919/104143

### [6] Ronald B. Morgan, “A Restarted GMRES Method Augmented with Eigenvectors” SIAM J. Matrix Anal. Appl., 16(4), 1154–1171. https://epubs.siam.org/doi/10.1137/S0895479893253975

### [7] Chow, E., Frommer, A. & Szyld, D.B. Asynchronous Richardson iterations: theory and practice. Numer Algor 87, 1635–1651 (2021). https://doi.org/10.1007/s11075-020-01023-3

### [8] Xijian Wang, “Comparison of Fixed Point Methods and Krylov Subspace Methods Solving Convection-Diffusion Equations.” American Journal of Computational Mathematics, (2015) 5, 113-126

