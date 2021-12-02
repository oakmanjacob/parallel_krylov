#include <vector>
#include <iostream>
#include <complex>
#include <cstring>
#include <string>
#include <tuple>
#include <unistd.h>
#include <chrono>

#include <spdlog/spdlog.h>
#include <parallel_krylov/gmres.h>
#include <parallel_krylov/matrix_import.h>

/**
 * @brief Print a helpful message about how to use this executable.
 * 
 * @param progname The name of the program from argv.
 */
void print_help(char *progname) {
    printf("%s: This executable is meant to demonstrate the performance of " 
           "parallel and sequential iterative solvers.\n", basename(progname));
    printf("  -k          The number corresponding to which iteration of the convection diffusion problem to import, ie A_$k.csv\n");
    printf("  -m          The type of matrix to generate for A (diag, tridiag, import)\n");
    printf("  -t          The number of threads to use\n");
    printf("  -n          The size of the square matrix to generate for A\n");
    printf("  -i          The max number of iterations for gmres to complete\n");
    printf("  -r          The number of iterations to do before restarting\n");
    printf("  -o          Use optimized GMRES implementation\n");
    printf("  -v          Enable verbose logs\n");
    printf("  -h          Print help (this message)\n");
}

/**
 * @brief This struct represents the configuration options specified
 * in the arguments passed to the executable.
 */
struct arg_t {
    size_t filenum = 1;
    std::string matrix_type = "diag";
    size_t thread_count = 1;
    size_t matrix_size = 5;
    size_t iterations = 100;
    size_t restart = 100;
    double tolerance = 0.000001;
    spdlog::level::level_enum log_level = spdlog::level::warn;
    bool help = false;
    double priority = 0.5;
    size_t optim_level = 1;
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
    while ((opt = getopt(argc, argv, "k:m:n:t:i:r:v::o:h:p:")) != -1) {
        switch (opt) {
            case 'k':
                args.filenum = atoi(optarg);
                break;
            case 'm':
                args.matrix_type = std::string(optarg);
                break;
            case 'n':
                args.matrix_size = atoi(optarg);
                break;
            case 't':
                args.thread_count = atoi(optarg);
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
            case 'o':
                args.optim_level = atoi(optarg);
                break;
            case 'p':
                args.priority = 1 / atoi(optarg);
                break;
            case 'h':
                args.help = true;
                break;
        }
    }
}

/**
 * @brief Convert a complex double to a string
 * 
 * @param value the complex double
 * @return string an a + bi representation of the complex value
 */
string to_string(complex<double> &value) {
    return to_string(real(value)) + " + " + to_string(imag(value)) + "i";
}

/**
 * @brief Convert a tuple representing an element in a sparse matrix to a string
 * 
 * @tparam V the element value type
 * @param value the element tuple of the form <y,x,value>
 * @return string a pretty string containing the element information
 */
template <typename V> string to_string(tuple<size_t,size_t,V> &value) {
    return "(y:" + to_string(std::get<0>(value)) + ", x:" + to_string(std::get<1>(value)) + ", value:" + to_string(std::get<2>(value)) + ")";
}

/**
 * @brief Convert a vector to a pretty string for printing
 * 
 * @tparam V the element type to print
 * @param vec the vector to convert to a string
 * @return string a pretty string containing the elements of the vector
 */
template <typename V> string to_string(vector<V> &vec) {
    string x_str = "\n{";
    for (auto & value : vec) {
        x_str += "\n" + to_string(value) + ",";
    }
    x_str += "\n}";
    return x_str;
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

    std::vector<std::tuple<size_t,size_t,double>> elements;
    std::vector<double> b;

    if (args.matrix_type == "diag") {
        spdlog::info("Generating size {} diagonal matrix with 2 along the diagonal", args.matrix_size);
        gen_diag(elements, args.matrix_size, 2.);
        b.resize(args.matrix_size, 1);
    }
    else if (args.matrix_type == "tridiag") {
        spdlog::info("Generating size {} tridiagonal matrix with 2 along the diagonal and -1 along the subdiagonal", args.matrix_size);
        gen_tridiag(elements, args.matrix_size);
        b.resize(args.matrix_size, 1);
    }
    else if (args.matrix_type == "import") {

        string path = "./data/n" + to_string(args.matrix_size);
        FILE* matrix_file = fopen((path + "/A_" + to_string(args.filenum) + ".csv").c_str(), "r");
        if (!importMatrixElements(matrix_file, elements, args.matrix_size)) {
            spdlog::error("Could not import file!");
        }

        if (matrix_file != nullptr) {
            fclose(matrix_file);
        }

        FILE* vector_file = fopen((path + "/b_" + to_string(args.filenum) + ".csv").c_str(), "r");
        if (!importVector(vector_file, b)) {
            spdlog::error("Could not import file!");
        }
        
        if (vector_file != nullptr) {
            fclose(vector_file);
        }
    }
    else {
        spdlog::error("Invalid matrix type");
        exit(1);
    }

    spdlog::debug("Sparse matrix A of size {} x {}: {}", args.matrix_size, args.matrix_size, to_string(elements));
    spdlog::debug("Result vector b of size {}: {}", b.size(), to_string(b));

    GMRES_In<double> input{
        MatrixCSR<double>(args.matrix_size, args.matrix_size, elements),
        b,
        vector<double>(args.matrix_size),
        args.tolerance,
        args.iterations,
        args.restart
    };
    
    GMRES_Out<double> output;

    spdlog::info("Starting GMRES Algorithm");
    
    auto t1 = std::chrono::high_resolution_clock::now();
    switch (args.optim_level) {
        case 4:
            wgmres(input, output);
            break;
        case 3:
            pgmres_sync(input, output, args.thread_count, args.priority);
            break;
        case 2:
            ogmres_simd(input, output);
            break;
        case 1:
            ogmres(input, output);
            break;
        case 0:
            gmres(input, output);
            break;
        default: {
            spdlog::error("Invalid optimization level {}", args.optim_level);
        }
    }
    auto t2 = std::chrono::high_resolution_clock::now();

    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1);

    printf("GMRES Algorithm finished in %zd iterations with convergence flag set to %d and a residual of magnitude: %e\n", output.iter, output.converged, norm(output.r));
    printf("Algorithm took %ld milliseconds\n", duration.count());

    spdlog::info("The following is the best guess for vector x {}", to_string(output.x));

    vector<double> b0;
    vector<double> r0;
    matvec(input.A, output.x, b0);
    vecsub(b0, input.b, r0);
    spdlog::info("The following is what we got for r0 {}", to_string(r0));

    return 0;
}