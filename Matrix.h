/*
 * Matrix.h
 *
 *Tyler Linderoth
 *Jan 27, 2015
 */

#ifndef MATRIX_H_
#define MATRIX_H_

#include <cstdio>

template <class T>
class Matrix
{
public:
        class Row
        {
                friend class Matrix;
        public:
                T& operator[] (size_t col);
                T operator[] (size_t col) const;
        private:
                Row(Matrix<T>& parent, size_t row);
                Matrix<T>& _parent;
                size_t _row;
        };
        Matrix (size_t rown = 0, size_t coln = 0, T val = 0);
        Matrix(const Matrix<T>& oldmat);
        ~Matrix();
        void allocate (size_t row, size_t col, T val = 0); // reserves memory for matrix
        T max (); // returns maximum value in matrix
        T min (); // returns minimum value in matrix
        size_t rown (); // number of rows
        size_t coln (); // number of columns
        size_t elemn (); // number of elements
        Row operator[] (size_t row);
        void setSubRow (size_t n);
        void setSubCol (size_t n);
        size_t getSubRow ();
        size_t getSubCol ();
private:
        size_t ncol;
        size_t nrow;
        size_t nelem;
        T** data;
        size_t _subncol;
        size_t _subnrow;
};
template<class T> Matrix<T>::Matrix (size_t rown, size_t coln, T val)
		: ncol(0),
		  nrow(0),
		  nelem(0),
		  data(0),
		  _subncol(0),
		  _subnrow (0)
{
	if (rown > 0 && coln > 0)
	{
		allocate(rown, coln, val);
	}
	else if (ncol == 0 && nrow == 0)
		nelem = 0;
	else
	{
		if (rown <= 0)
			fprintf(stderr, "Row number not a positive integer in matrix object...\n");
		if (coln <= 0)
			fprintf(stderr, "Column number not a positive integer in matrix object...\n");
		fprintf(stderr, "Could not allocate memory for matrix object\n");
		nelem = 0;
	}
}

template<class T> Matrix<T>::~Matrix ()
{
	for (size_t i = 0; i < nrow; i++)
		delete [] data[i];
	delete [] data;
	data = 0;
}

template <class T> T& Matrix<T>::Row::operator[] (size_t col)
{
	return _parent.data[_row][col];
}

template <class T> T Matrix<T>::Row::operator[] (size_t col) const
{
	return _parent.data[_row][col];
}

template <class T> Matrix<T>::Row::Row(Matrix<T>& parent, size_t row)
                        : _parent(parent),
                          _row(row)
{}


template<class T> typename Matrix<T>::Row Matrix<T>::operator[] (size_t row)
{
	return Row(*this, row);
}

template<class T> void Matrix<T>::allocate (size_t row, size_t col, T val)
{
	data = new T* [row];
	size_t i = 0;
	size_t j = 0;

	for (i = 0; i < row; ++i)
		data[i] = new T [col];

	// initialize values
	for (i = 0; i < row; ++i)
		for (j = 0; j < col; ++j)
			data[i][j] = val;
	nrow = row;
	ncol = col;
	nelem = row * col;
}

template <class T> Matrix<T>::Matrix(const Matrix<T>& oldmat)
	: ncol(oldmat.ncol),
	  nrow(oldmat.nrow),
	  nelem(oldmat.nelem)
{
	size_t i = 0;
	size_t j = 0;
	for(i = 0; i < nrow; ++i)
		data[i] = new T[ncol];
	for(i = 0; i < nrow; ++i)
		for (j = 0; j < ncol; ++j)
			data[i][j] = oldmat.data[i][j];
}


template<class T> T Matrix<T>::max ()
{
	T max = data[0][0];

	for (size_t i = 0; i < nrow; i++)
		for (size_t j = 0; j < ncol; j++)
		{
			if (data[i][j] > max)
				max = data[i][j];
		}

	return max;
}

template<class T> T Matrix<T>::min ()
{
	T min = data[0][0];

	for (size_t i = 0; i < nrow; i++)
		for (size_t j = 0; j < ncol; j++)
		{
			if (data[i][j] < min)
				min = data[i][j];
		}

	return min;
}

template<class T> size_t Matrix<T>::rown ()
{
	return nrow;
}

template<class T> size_t Matrix<T>::coln ()
{
	return ncol;
}

template<class T> size_t Matrix<T>::elemn ()
{
	return nelem;
}

template<class T> void Matrix<T>::setSubRow (size_t n)
{
	if (n > nrow)
	{
		fprintf(stderr, "row subset size exceeds matrix row capacity\n");
		_subnrow = 0;
	}
	else
		_subnrow = n;
}

template<class T> void Matrix<T>::setSubCol (size_t n)
{
	if (n > ncol)
	{
		fprintf(stderr, "column subset size exceeds matrix column capacity\n");
		_subncol = 0;
	}
	else
		_subncol = n;
}

template<class T> size_t Matrix<T>::getSubRow ()
{
	return _subnrow;
}

template<class T> size_t Matrix<T>::getSubCol ()
{
	return _subncol;
}

#endif /* MATRIX_H_ */
