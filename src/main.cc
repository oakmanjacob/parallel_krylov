#include <vector>
#include <iostream>
#include <complex>
#include <cstring>
#include <tuple>
#include <unistd.h>

#include <parallel_krylov/gmres.h>

/**
 * @brief Print a helpful message about how to use this executable.
 * 
 * @param progname The name of the program from argv.
 */
void print_help(char *progname) {
    printf("%s: This executable is meant to demonstrate the performance of " 
           "parallel and sequential iterative solvers.\n", basename(progname));
    printf("  -h          Print help (this message)\n");
}

/**
 * @brief This struct represents the configuration options specified
 * in the arguments passed to the executable.
 */
struct arg_t {

    // Whether the user has requested the help message
    bool help = false;
};

/**
 * @brief Parse the commandline arguments and use them to populate the
 * provided args object
 * 
 * @param argc The number of command-line arguments passed to the program
 * @param argv The list of command-line arguments
 * @param args The struct into which the parsed args should go
 */
void parse_args(int argc, char **argv, arg_t &args) {
    long opt;
    while ((opt = getopt(argc, argv, "h")) != -1) {
        switch (opt) {
            case 'h':
                args.help = true;
                break;
        }
    }
}

/**
 * @brief Main entry point to the executable
 * 
 * @param argc The number of command-line arguments passed to the program
 * @param argv The list of command-line arguments
 * @return int success or error
 */
int main(int argc, char** argv) {
    arg_t args;
    parse_args(argc, argv, args);
    if (args.help) {
        print_help(argv[0]);
        exit(0);
    }

    auto elements = std::vector<std::tuple<size_t,size_t,complex<double>>>{
        std::make_tuple(0,0,2),
        std::make_tuple(1,1,2),
        std::make_tuple(2,2,2),
        std::make_tuple(3,3,2),
        std::make_tuple(4,4,2)
    };

    MatrixCSR<complex<double>> A = MatrixCSR<complex<double>>(5, 5, elements);

    vector<complex<double>> b(5, 1);
    vector<complex<double>> x0(5);
    double tol = 0.000001;
    size_t max_it = 10;
    size_t restart = 10;
    vector<complex<double>> x;
    vector<complex<double>> r;
    vector<double> r_nrm;
    size_t iter;
    bool converged;

    gmres(A,b,x0,tol,max_it,restart,x,r,r_nrm,iter,converged);

    std::cout << "Howdy Partner!" << std::endl;
    for (auto & value : x) {
        std::cout << value << std::endl;
    }

    return 0;
}