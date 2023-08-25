//
// Created by joaovictor on 20/08/23.
//

#ifndef GERVLIB_MATRIX_H
#define GERVLIB_MATRIX_H

#include "Serialize.h"
#include "Configure.h"
#include <format>
#include <iostream>

namespace gervLib::matrix
{

    template <typename O, typename T>
    class Matrix : serialize::Serialize
    {
    private:
        O rows{}, cols{};
        std::unique_ptr<T[]> data;

    public:
        Matrix() = default;

        Matrix(O rows, O cols) : rows(rows), cols(cols)
        {
            data = std::make_unique<T[]>(rows * cols);
        }

        virtual ~Matrix()
        {
            if (data != nullptr)
                data.reset();
        }

        void setValue(O _row, O _col, T _value)
        {

            if (data == nullptr)
                throw std::invalid_argument("Matrix::setValue: Matrix not initialized");

            if (_row >= this->rows || _col >= this->cols || _row < 0 || _col < 0)
                throw std::invalid_argument("Matrix::setValue: Invalid row or column");

            data[_row * cols + _col] = _value;
        }

        T getValue(O _row, O _col)
        {
            if (data == nullptr)
                throw std::invalid_argument("Matrix::getValue: Matrix not initialized");

            if (_row >= this->rows || _col >= this->cols || _row < 0 || _col < 0)
                throw std::invalid_argument("Matrix::getValue: Invalid row or column");

            return data[_row * cols + _col];
        }

        bool isEqual(std::unique_ptr<Matrix<O, T>>& other)
        {

            if ((data == nullptr && other != nullptr) || (data != nullptr && other == nullptr))
                return false;

            if (this->rows != other->rows || this->cols != other->cols)
                return false;

            if (data == nullptr && other == nullptr)
                return true;
            else
            {
                for (O i = 0; i < rows; i++)
                {
                    for (O j = 0; j < cols; j++)
                    {
                        if (data[i * cols + j] != other->data[i * cols + j])
                            return false;
                    }
                }
            }

            return true;
        }

        bool operator==(std::unique_ptr<Matrix<O, T>>& other)
        {
            return isEqual(std::move(other));
        }

        bool operator!=(std::unique_ptr<Matrix<O, T>>& other)
        {
            return !isEqual(std::move(other));
        }

        friend std::ostream &operator<<(std::ostream &_os, Matrix &_matrix) {

            for (O i = 0; i < _matrix.rows; i++) {
                for (O j = 0; j < _matrix.cols; j++) {
                    _os << std::format(" {1:^{0}} ", configure::double_size, _matrix.getValue(i, j));
                }
                _os << std::endl;
            }

            return _os;

        }

        std::unique_ptr<u_char[]> serialize() override
        {

            std::unique_ptr<u_char[]> _data = std::make_unique<u_char[]>(getSerializedSize());

            memcpy(_data.get(), &this->rows, sizeof(O));
            memcpy(_data.get() + sizeof(O), &this->cols, sizeof(O));
            memcpy(_data.get() + sizeof(O) * 2, this->data.get(), sizeof(T) * this->rows * this->cols);

            return _data;

        }

        void deserialize(std::unique_ptr<u_char[]> _data) override
        {

            memcpy(&this->rows, _data.get(), sizeof(O));
            memcpy(&this->cols, _data.get() + sizeof(O), sizeof(O));

            clear();
            data = std::make_unique<T[]>(rows * cols);

            memcpy(this->data.get(), _data.get() + sizeof(O) * 2, sizeof(T) * this->rows * this->cols);

            _data.reset();

        }

        size_t getSerializedSize() override
        {
            return sizeof(O) * 2 + sizeof(T) * this->rows * this->cols;
        }

        void clear()
        {
            if (data != nullptr)
                data.reset();
        }


    };

}

#endif //GERVLIB_MATRIX_H
