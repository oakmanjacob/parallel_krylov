#include <cassert>
#include <vector>
#include <algorithm>
#include <tuple>

#include "parallel_krylov/matrix.h"

template <typename V> class MatrixCSR {
public:
    std::vector<V> _values;
    std::vector<size_t> _col_indexes;
    std::vector<size_t> _row_pointers;
    size_t _col_count;

    MatrixCSR<V>(size_t row_count, size_t col_count) {
        _col_count = col_count;
        _row_pointers.resize(row_count);
    }

    MatrixCSR<V>(Matrix<V> other) {
        _col_count = other.get_col_count();
        _row_pointers.resize(other.get_row_count());
    }

    ~MatrixCSR<V>() {}

    size_t get_row_count() {
        return _row_pointers.size();
    }

    size_t get_col_count() {
        return _col_count;
    }

    void set(size_t x, size_t y, V value) {
        for (size_t i = x + 1; i < _row_pointers.size(); i++) {
            _row_pointers[i]++;
        }
        std::vector<size_t>::iterator location;
        if (x + 1 < _row_pointers.size()) {
            location = lower_bound(_col_indexes.begin() + _row_pointers[x], _col_indexes.begin() + _row_pointers[x + 1], y);
        }
        else {
            location = lower_bound(_col_indexes.begin() + _row_pointers[x], _col_indexes.end(), y); 
        }

        if (location != _col_indexes.end() && *location != y) {
            _col_indexes.insert(location, y);
            _values.insert(_values.begin() + location - _col_indexes.begin(), value);
        }
        else {
            _values[location - _col_indexes.begin()] = value;
        }
    }

    void replace(std::vector<std::tuple<size_t,size_t,V>> &elements) {
        // Remove all old values from the matrix;
        _values.clear();
        _values.reserve(elements.size());
        _col_indexes.clear();
        _col_indexes.reserve(elements.size());


        // Sort the elements in ascending order based on x and then y
        std::sort(elements.begin(), elements.end(), [](std::tuple<size_t,size_t,V> a, std::tuple<size_t,size_t,V> b){
            return std::get<0>(a) < std::get<0>(b) ||
                (std::get<0>(a) == std::get<0>(b) && std::get<1>(a) < std::get<1>(b));
        });

        _row_pointers[0] = 0;

        size_t row_index = 0;
        size_t prev_x = -1;
        size_t prev_y = -1;

        for (auto & element : elements) {
            if (prev_x == std::get<0>(element) && prev_y == std::get<0>(element)) {
                continue;
            }

            if (std::get<0>(element) >= _row_pointers.size() || std::get<1>(element) >= _col_count) {
                continue;
            }

            while (row_index < std::get<0>(element) && row_index + 1 < _row_pointers.size()) {
                row_index++;
                _row_pointers[row_index] = _values.size();
            }

            _col_indexes.push_back(std::get<1>(element));
            _values.push_back(std::get<2>(element));

            prev_x = std::get<0>(element);
            prev_y = std::get<1>(element);
        }
    }
};