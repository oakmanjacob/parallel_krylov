#include <cassert>
#include <vector>
#include <iterator>

template <typename V> class Matrix {
public:
    std::vector<std::vector<V>> _data;
    size_t _row_count;
    size_t _col_count;

    Matrix<V>(size_t row_count, size_t col_count) {
        _row_count = row_count;
        _col_count = col_count;

        _data = new std::vector<std::vector<V>>(_row_count);
        for (size_t i = 0; i < _row_count; i++) {
            _data[i].resize(_col_count);
        }
    }

    Matrix<V>(size_t row_count, size_t col_count, V **data) {
        _row_count = row_count;
        _col_count = col_count;
        _data = new std::vector<std::vector<V>>(_row_count);
        for (size_t i = 0; i < _row_count; i++) {
            _data[i].resize(_col_count);
        }

        for (int i = 0; i < _row_count; i++) {
            for (int j = 0; j < _col_count; j++) {
                _data[i][j] = data[i][j];
            }
        }
    }

    Matrix<V>(const Matrix<V> &other) {
        _row_count = other._row_count();
        _col_count = other._col_count();
        _data = new std::vector<std::vector<V>>(_row_count);
        for (size_t i = 0; i < _row_count; i++) {
            _data[i].resize(_col_count);
        }

        for (int i = 0; i < _row_count; i++) {
            for (int j = 0; j < _col_count; j++) {
                _data[i][j] = other._data[i][j];
            }
        }
    }

    ~Matrix<V>() {}
};

// /**
//  * @brief Define Matvec operation Ax
//  * 
//  * @tparam V the type of values for both 
//  * @param mat a dense matrix
//  * @param vec a dense vector
//  * @return std::vector<V> the result of Ax
//  */
// template <typename V> std::vector<V>& operator* (const Matrix<V> &mat, const std::vector<V> &vec) {
//     assert(vec.size() == mat._col_count);

//     std::vector<V> result(mat._row_count);
//     for (size_t i = 0; i < mat._row_count; i++) {
//         for (size_t j = 0; j < mat._col_count; j++) {
//             result[i] += mat[i][j] * vec[j];
//         }
//     }

//     return result;
// }