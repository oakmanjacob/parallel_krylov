#include <vector>
#include <tuple>
#include <complex>
#include <cstdio>

bool importMatrixElements(FILE* matrixFile, std::vector<std::tuple<size_t,size_t,std::complex<double>>> &elements, size_t &n) {
    if (matrixFile == nullptr) {
        return false;
    }
    n = 0;
    elements.clear();
    size_t row_id;
    size_t col_id;
    double value;
    while (fscanf(matrixFile, "%li,%li,%lf\n", &row_id, &col_id, &value) == 3) {
        elements.emplace_back(std::make_tuple(row_id - 1, col_id - 1, std::complex<double>(value)));
        n = max(max(n, row_id), col_id);
    }

    return true;
}

bool importVector(FILE* vectorFile, std::vector<std::complex<double>> &b) {
    if (vectorFile == nullptr) {
        return false;
    }
    
    b.clear();
    double value;
    while (fscanf(vectorFile, "%lf\n", &value) == 1) {
        b.emplace_back(std::complex<double>(value));
    }

    return true;
}