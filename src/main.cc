#include <vector>
#include <iostream>
#include <complex>
#include <cstring>
#include <string>
#include <tuple>
#include <unistd.h>

#include <spdlog/spdlog.h>
#include <parallel_krylov/gmres.h>

/**
 * @brief Print a helpful message about how to use this executable.
 * 
 * @param progname The name of the program from argv.
 */
void print_help(char *progname) {
    printf("%s: This executable is meant to demonstrate the performance of " 
           "parallel and sequential iterative solvers.\n", basename(progname));
    printf("  -f          The file name containing the matrix information to import\n");
    printf("  -t          The type of matrix to generate for A (diag, tridiag)\n");
    printf("  -n          The size of the square matrix to generate for A\n");
    printf("  -i          The max number of iterations for gmres to complete\n");
    printf("  -r          The number of iterations to do before restarting\n");
    printf("  -v          Enable verbose logs\n");
    printf("  -h          Print help (this message)\n");
}

/**
 * @brief This struct represents the configuration options specified
 * in the arguments passed to the executable.
 */
struct arg_t {
    std::string filename;
    std::string matrix_type = "diag";
    size_t matrix_size = 5;
    size_t iterations = 100;
    size_t restart = 100;
    double tolerance = 0.000001;
    spdlog::level::level_enum log_level = spdlog::level::warn;
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
    while ((opt = getopt(argc, argv, "f:t:n:i:r:v::h:")) != -1) {
        switch (opt) {
            case 'f':
                args.filename = std::string(optarg);
                break;
            case 't':
                args.matrix_type = std::string(optarg);
                break;
            case 'n':
                args.matrix_size = atoi(optarg);
                break;
            case 'i':
                args.iterations = atoi(optarg);
                break;
            case 'r':
                args.restart = atoi(optarg);
                break;
            case 'v':
                args.log_level = spdlog::level::info;
                break;
            case 'h':
                args.help = true;
                break;
        }
    }
}

void gen_diag(std::vector<std::tuple<size_t,size_t,complex<double>>> &elements, size_t n, complex<double> value) {
    elements.clear();
    elements.reserve(n);

    for (size_t i = 0; i < n; i++) {
        elements.emplace_back(std::make_tuple(i,i,value));
    }
}

void gen_tridiag(std::vector<std::tuple<size_t,size_t,complex<double>>> &elements, size_t n) {
    elements.clear();
    elements.reserve(3*(n-1) + 1);

    for (size_t i = 0; i < n; i++) {
        if (i > 0)
            elements.emplace_back(std::make_tuple(i,i-1,-1));

        elements.emplace_back(std::make_tuple(i,i,2));

        if (i + 1 < n)
            elements.emplace_back(std::make_tuple(i,i+1,-1));
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

    spdlog::set_level(args.log_level);
    spdlog::info("Parallel Krylov Program Started");

    std::vector<std::tuple<size_t,size_t,complex<double>>> elements;

    if (args.matrix_type == "diag") {
        spdlog::info("Generating size {} diagonal matrix with 2 along the diagonal", args.matrix_size);
        gen_diag(elements, args.matrix_size, 2);
    }
    else if (args.matrix_type == "tridiag") {
        spdlog::info("Generating size {} tridiagonal matrix with 2 along the diagonal and -1 along the subdiagonal", args.matrix_size);
        gen_tridiag(elements, args.matrix_size);
    }
    else {
        spdlog::error("Invalid matrix type");
        exit(1);
    }

    spdlog::info("Sparse matrix A:");
    for (auto & element : elements) {
        spdlog::info("row: {}, col: {}, value: {}", std::get<0>(element), std::get<1>(element), real(std::get<2>(element)));
    }

    GMRES_In input{
        MatrixCSR<complex<double>>(args.matrix_size, args.matrix_size, elements),
        vector<complex<double>>(args.matrix_size, 1),
        vector<complex<double>>(args.matrix_size),
        args.tolerance,
        args.iterations,
        args.restart
    };
    
    GMRES_Out output{
        vector<complex<double>>(),
        vector<complex<double>>(),
        vector<double>(),
        0,
        false
    };

    spdlog::info("Starting GMRES Algorithm");
    gmres(input, output);
    spdlog::info("GMRES Algorithm finished in {} iterations with convergence flag set to: {}", output.iter, output.converged);

    spdlog::info("The following is the best guess for vector x with a residual of magnitude: {}", norm(output.r));
    for (auto & value : output.x) {
        std::cout << value << std::endl;
    }

    return 0;
}