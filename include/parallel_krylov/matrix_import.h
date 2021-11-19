#include <vector>
#include <tuple>
#include <complex>

/**
 * @brief Generate the elements for a sparse diagonal matrix with value across the diagonal
 * 
 * @tparam V the type to use for the values in the matrix
 * @param elements the return vector of matrix elements of the form <y,x,value>
 * @param n the size of the matrix to generate
 * @param value the value to place along the diagonal
 */
template <typename V> void gen_diag(std::vector<std::tuple<size_t,size_t,V>> &elements, size_t n, V value) {
    elements.clear();
    elements.reserve(n);

    for (size_t i = 0; i < n; i++) {
        elements.emplace_back(std::make_tuple(i,i,value));
    }
}

/**
 * @brief Generate the elements for a sparce tridiagonal matrix with 2 across the diagonal and
 * -1 across the subdiagonals
 * 
 * @tparam V the type to use for the values in the matrix
 * @param elements the return vector of matrix elements of the form <y,x,value>
 * @param n the size of the matrix to generate
 */
template <typename V> void gen_tridiag(std::vector<std::tuple<size_t,size_t,V>> &elements, size_t n) {
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
 * @brief Import the elements of a sparse matrix from a csv file
 * 
 * @tparam V the type to use for the values of this matrix
 * @param matrixFile a pointer to the open file to import from
 * @param elements the return vector of matrix elements of the form <y,x,value>
 * @param n the return value of the size of the matrix imported
 * @return true on successful import
 * @return false if the matrixFile pointer was null
 */
template <typename V> bool importMatrixElements(FILE* matrixFile, std::vector<std::tuple<size_t,size_t,V>> &elements, size_t &n) {
    if (matrixFile == nullptr) {
        return false;
    }
    n = 0;
    elements.clear();
    size_t row_id;
    size_t col_id;
    double value;
    while (fscanf(matrixFile, "%li,%li,%lf\n", &row_id, &col_id, &value) == 3) {
        elements.emplace_back(std::make_tuple(row_id - 1, col_id - 1, V(value)));
        n = max(max(n, row_id), col_id);
    }

    return true;
}

/**
 * @brief Import a dense vector from a csv file
 * 
 * @tparam V the type to use for the values of this vector
 * @param vectorFile a pointer to the open file to import from
 * @param b the return vector
 * @return true on successful import
 * @return false if the vectorFile pointer was null
 */
template <typename V> bool importVector(FILE* vectorFile, std::vector<V> &b) {
    if (vectorFile == nullptr) {
        return false;
    }
    
    b.clear();
    double value;
    while (fscanf(vectorFile, "%lf\n", &value) == 1) {
        b.emplace_back(V(value));
    }

    return true;
}