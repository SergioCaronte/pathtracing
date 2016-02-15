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
#ifndef COLOR_H
#define COLOR_H

#include <iostream>
#include <math.h>

class Color
{
public:
	float r, g, b;

    Color(float _r = 0.0f, float _g = 0.0f, float _b = 0.0f)
        : r(_r), g(_g), b(_b)
	{}

    void clamp()
	{
		r = std::max(0.0f,std::min(r ,1.0f));
		g = std::max(0.0f,std::min(g ,1.0f));
		b = std::max(0.0f,std::min(b ,1.0f));
    }

    inline Color &invScale(float factor)
    {
        r /= factor;
        g /= factor;
        b /= factor;
        return *this;
    }


    Color &operator+=(float other)
	{
        r += other;
        g += other;
        b += other;
        return *this;
    }

    Color &operator+=(const Color &other)
	{
        r += other.r;
        g += other.g;
        b += other.b;
        return *this;
    }

    Color &operator*=(float other)
	{
        r *= other;
        g *= other;
        b *= other;
        return *this;
    }

    Color &operator*=(const Color &other)
	{
        r *= other.r;
        g *= other.g;
        b *= other.b;
        return *this;
    }

    Color &operator/=(const float &factor)
	{
        this->invScale(factor);
        return *this;
    }
};

inline Color operator+(Color lhv, float rhv)
{
    lhv += rhv;
    return lhv;
}

inline Color operator+(Color lhv, const Color &rhv)
{
    lhv += rhv;
    return lhv;
}

inline Color operator*(Color lhv, float rhv)
{
    lhv *= rhv;
    return lhv;
}

inline Color operator*(Color lhv, const Color &rhv)
{
    lhv *= rhv;
    return lhv;
}

inline Color operator/(Color lhv, float rhv)
{
    lhv /= rhv;
    return lhv;
}

inline std::ostream &operator<<(std::ostream &stream, const Color &color)
{
    stream << (int) (color.r * 255) << " " << (int) (color.g * 255) << " "<< (int) (color.b * 255) << std::endl;
    return stream;
}



#endif
