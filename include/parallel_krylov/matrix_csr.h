#include <cassert>
#include <vector>
#include <algorithm>
#include <tuple>

#include <spdlog/spdlog.h>

template <typename V> class MatrixCSR {
public:
    std::vector<V> _values;
    std::vector<size_t> _col_indexes;
    std::vector<size_t> _row_pointers;
    size_t _col_count;

    /**
     * @brief Construct a new Matrix C S R< V> object
     * 
     * @param row_count 
     * @param col_count 
     */
    MatrixCSR<V>(size_t row_count, size_t col_count) {
        _col_count = col_count;
        _row_pointers.resize(row_count);
    }

    /**
     * @brief Construct a new Matrix C S R< V> object
     * 
     * @param row_count 
     * @param col_count 
     * @param elements 
     */
    MatrixCSR<V>(size_t row_count, size_t col_count, std::vector<std::tuple<size_t,size_t,V>> &elements) {
        _col_count = col_count;
        _row_pointers.resize(row_count);

        this->replace(elements);
    }

    /**
     * @brief Destroy the Matrix C S R< V> object
     * 
     */
    ~MatrixCSR<V>() {}

    /**
     * @brief Get the row count object
     * 
     * @return size_t 
     */
    size_t get_row_count() {
        return _row_pointers.size();
    }

    /**
     * @brief Get the col count object
     * 
     * @return size_t 
     */
    size_t get_col_count() {
        return _col_count;
    }

    /**
     * @brief 
     * 
     * @param y 
     * @param x 
     * @param value 
     */
    void set(size_t y, size_t x, V value) {
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

    /**
     * @brief 
     * 
     * @param elements 
     */
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
        size_t prev_row = -1;
        size_t prev_col = -1;

        for (auto & element : elements) {
            if (prev_row == std::get<0>(element) && prev_col == std::get<1>(element)) {
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

            prev_row = std::get<0>(element);
            prev_col = std::get<1>(element);
        }
    }
};