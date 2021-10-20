#include <cassert>
#include <vector>

template <typename V> class Matrix {
    V _data[][];
    size_t _row_count;
    size_t _col_count;

public:
    Matrix(size_t row_count, size_t col_count) {
        _row_count = row_count;
        _col_count = col_count;
        _data = new V[row_count][col_count]
    }

    ~Matrix() {
        for (size_t i = 0; i < _row_count; i++) {
            delete[] _data[i];
        }

        delete[] _data;
    }

    size_t get_row_count() {
        return _row_count;
    }

    size_t get_col_count() {
        return _col_count;
    }
};

/**
 * @brief Define Matvec operation
 * 
 * @tparam V the type of values for both 
 * @param mat 
 * @param vec 
 * @return std::vector<V> 
 */
template <typename V> std::vector<V> operator* (const Matrix<V> mat, std::vector<V> vec) {
    assert(vec.size() == mat._col_count);

    std:vector<V> result(mat._row_count);
    for (size_t i = 0; i < mat._row_count; i++) {
        for (size_t j = 0; j < mat._col_count; j++) {
            result.at(i) += mat[i][j] * vec.at(j);
        }
    }

    return result;
}