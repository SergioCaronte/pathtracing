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

#ifndef POINT_H
#define POINT_H

#include <iostream>
#include "vector.h"


class Point
{

public:
    float x, y, z, w;

    Point(float x = 0.0f, float y = 0.0f, float z = 0.0f, float w = 1.0f)
        : x(x), y(y), z(z), w(w)
	{}

    static double distance(const Point &p1, const Point &p2)
	{
        return std::sqrt(std::pow(p1.x - p2.x, 2)
                + std::pow(p1.y - p2.y, 2)
                + std::pow(p1.z - p2.z, 2));
    }

    Point &operator+=(const Point &rhv)
	{
        this->x += rhv.x;
        this->y += rhv.y;
        this->z += rhv.z;
        this->w += rhv.w;
        return *this;
    }

    Point &operator+=(const Vector &rhv)
	{
        this->x += rhv.x;
        this->y += rhv.y;
        this->z += rhv.z;
        this->w += rhv.w;
        return *this;
    }

    Point &operator-=(const Vector &rhv)
	{
        this->x -= rhv.x;
        this->y -= rhv.y;
        this->z -= rhv.z;
        this->w -= rhv.w;
        return *this;
    }

    Point &operator+=(const float &rhv)
	{
        this->x += rhv;
        this->y += rhv;
        this->z += rhv;
        this->w += rhv;
        return *this;
    }

    Point &operator*=(const float &rhv)
	{
        this->x *= rhv;
        this->y *= rhv;
        this->z *= rhv;
        this->w *= rhv;
        return *this;
    }
};


/// + operator for point-vector addition.
inline Point operator+(Point lhv, const Vector &rhv)
{
    lhv += rhv;
    return lhv;
}

inline Point operator-(Point lhv, const Vector &rhv)
{
    lhv -= rhv;
    return lhv;
}

inline Point operator*(Point lhv, const float &rhv)
{
    lhv.x *= rhv;
    lhv.y *= rhv;
    lhv.z *= rhv;
    return lhv;
}

inline Vector operator-(const Point &lhv, const Point &rhv)
{
    return Vector(lhv.x - rhv.x, lhv.y - rhv.y, lhv.z - rhv.z, lhv.w - rhv.w);
}

inline bool operator==(const Point &lhv, const Point &rhv)
{
    return(lhv.x == rhv.x && lhv.y == rhv.y && lhv.z == rhv.z && lhv.w == rhv.w);
}

inline bool operator!=(const Point &lhv, const Point &rhv)
{
    return !operator==(lhv, rhv);
}


inline std::ostream &operator<<(std::ostream &stream, const Point &pt)
{
    stream << "(" << pt.x << ", " << pt.y << ", " << pt.z << ", " << pt.w << ")";
    return stream;
}

#endif
