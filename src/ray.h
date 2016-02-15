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

#ifndef RAY_H
#define RAY_H

#include "math/point.h"
#include "math/vector.h"

class Ray
{
public:
	Ray(){}
	Ray(Point o, Vector d);
	~Ray(){}

	Point origin;
	Vector direction;
};

inline std::ostream &operator<<(std::ostream &stream, const Ray &ray) 
{
	stream << "origin: (" <<ray.origin.x << ", " << ray.origin.y << ", " << ray.origin.z << ") direction: ("
		<< ray.direction.x << ", " << ray.direction.y << ", " << ray.direction.z << ")";
    return stream;
}

#endif 