#include <cassert>
#include <vector>

template <typename V> class MatrixCSR {
    std::vector<V> _values;
    std::vector<V> _col_indexes;
    std::vector<size_t> _row_pointers;
    size_t _col_count;

public:
    MatrixCSR(size_t row_count, size_t col_count) {
        _col_count = col_count;
        _row_pointers.resize(row_count);
    }

    ~MatrixCSR() {}

    size_t get_row_count() {
        return _row_pointers.size();
    }

    size_t get_col_count() {
        return _col_count;
    }
};

/**
 * @brief Define Matvec operation Ax
 * 
 * @tparam V the type of values for both 
 * @param mat a dense matrix
 * @param vec a dense vector
 * @return std::vector<V> the result of Ax
 */
template <typename V> std::vector<V>& operator* (MatrixCSR<V> &mat, std::vector<V> &vec) {
    assert(vec.size() == mat._col_count);

    std:vector<V> result(mat.get_row_count());
    size_t row_index = 0;
    for (size_t i = 0; i < mat._row_count; i++) {
        while (row_index + 1 < mat._row_pointers.size() && mat._row_pointers[row_index + 1] <= i) {
            row_index++;
        }

        result[row_index] += mat._values[i] * vec[mat._col_indexes[i]];
    }

    return result;
}