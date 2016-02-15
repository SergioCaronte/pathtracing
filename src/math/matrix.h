/*
* Copyright (C) 2015 Sergio Nunes da Silva Junior 
*
* Free depedency c++ path tracing algorithm
* Assignment of Advanced Computer Graphic Course - 2/2015
*
* This program is free software; you can redistribute it and/or modify it
* under the terms of the GNU General Public License as published by the Free
* Software Foundation; either version 2 of the License.
*
* This program is distributed in the hope that it will be useful, but WITHOUT
* ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
* FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for
* more details.
*
* author: Sergio Nunes da Silva Junior
* contact: sergio.nunes@dcc.ufmg.com.br
* Universidade Federal de Minas Gerais (UFMG) - Brazil
*/

#ifndef MATRIX_H
#define MATRIX_H

#include "math.h"

#include <cstring>
#include <new>
#include <stdexcept>
#include <cfloat>
#include <string>
#include <sstream>

class Matrix4x4
{
    float m[4][4];

public:

	inline void identity()
	{
		 for(size_t i = 0; i < 4; ++i)
            for(size_t j = 0; j < 4; ++j)
			{
				if(i == j)
					m[i][j] = 1;
				else
					m[i][j] = 0;
            }
	}

	inline void fill(float v)
	{
		 for(size_t i = 0; i < 4; ++i)
            for(size_t j = 0; j < 4; ++j)
                m[i][j] = v;
	}

    inline float *operator[](int index)
	{
        return m[index];
    }


    inline const float *operator[](int index) const
	{
        return m[index];
    }


    inline float &at(int row, int column)
	{
        if(row >= 4 || column >= 4)
            throw std::out_of_range("Trying to access invalid matrix position.");
        return m[row][column];
    }

    inline const float &at(int row, int column) const
	{
        if(row >= 4 || column >= 4)
            throw std::out_of_range("Trying to access invalid matrix position.");
        return m[row][column];
    }

    Matrix4x4 &operator+=(const Matrix4x4 &right)
	{
        for(int i = 0; i < 4; ++i)
            for(int j = 0; j < 4; ++j)
                m[i][j] += right[i][j];
        return *this;
    }


    Matrix4x4 &operator-=(const Matrix4x4 &right)
	{
        for(int i = 0; i < 4; ++i)
            for(int j = 0; j < 4; ++j)
                m[i][j] -= right[i][j];
        return *this;
    }

    Matrix4x4 &operator*=(float right)
	{
        for(int i = 0; i < 4; ++i)
            for(int j = 0; j < 4; ++j)
                m[i][j] *= right;
        return *this;
    }


    Matrix4x4 &operator/=(float right)
	{
        for(int i = 0; i < 4; ++i)
            for(int j = 0; j < 4; ++j)
                m[i][j] /= right;
        return *this;
    }


    Matrix4x4 &operator*=(const Matrix4x4 &right)
	{
        for(size_t i = 0; i < 4; ++i)
            for(size_t j = 0; j < 4; ++j)
                for(size_t k = 0; k < 4; ++k)
                    m[i][k] *= right[k][j];
        return *this;
    }

	std::string print()
	{
		std::stringstream stream;
		for(size_t i = 0; i < 4; ++i)
		{
			for(size_t j = 0; j < 4; ++j)
			{
				stream << m[i][j] << "\t";
			}
			stream << "\n";
		}
		return stream.str();
	}
};

inline Matrix4x4 operator*( Matrix4x4 &left,  Matrix4x4 &right)
{
    Matrix4x4 m;
    m.fill(0);

    for(int i = 0; i < 4; ++i)
        for(int j = 0; j < 4; ++j)
            for(int k = 0; k < 4; ++k)
                m[i][j] += left[i][k] * right[k][j];
    return m;
}

inline Matrix4x4 operator*(const Matrix4x4 &left,const  Matrix4x4 &right)
{
    Matrix4x4 m;
    m.fill(0);

    for(int i = 0; i < 4; ++i)
        for(int j = 0; j < 4; ++j)
            for(int k = 0; k < 4; ++k)
                m[i][j] += left[i][k] * right[k][j];
    return m;
}


#endif
