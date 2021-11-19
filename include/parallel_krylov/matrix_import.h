#include <vector>
#include <tuple>
#include <complex>
#include <cstdio>

void gen_diag(std::vector<std::tuple<size_t,size_t,complex<double>>> &elements, size_t n, complex<double> value) {
    elements.clear();
    elements.reserve(n);

    for (size_t i = 0; i < n; i++) {
        elements.emplace_back(std::make_tuple(i,i,value));
    }
}

void gen_diag(std::vector<std::tuple<size_t,size_t,double>> &elements, size_t n, double value) {
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

void gen_tridiag(std::vector<std::tuple<size_t,size_t,double>> &elements, size_t n) {
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

bool importMatrixElements(FILE* matrixFile, std::vector<std::tuple<size_t,size_t,double>> &elements, size_t &n) {
    if (matrixFile == nullptr) {
        return false;
    }
    n = 0;
    elements.clear();
    size_t row_id;
    size_t col_id;
    double value;
    while (fscanf(matrixFile, "%li,%li,%lf\n", &row_id, &col_id, &value) == 3) {
        elements.emplace_back(std::make_tuple(row_id - 1, col_id - 1, value));
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

bool importVector(FILE* vectorFile, std::vector<double> &b) {
    if (vectorFile == nullptr) {
        return false;
    }
    
    b.clear();
    double value;
    while (fscanf(vectorFile, "%lf\n", &value) == 1) {
        b.emplace_back(value);
    }

    return true;
}