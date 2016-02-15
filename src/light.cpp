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

#include "light.h"
#include "multijittered.h"
#include <string>

Light::Light()
	: type(LightAreaType::LAT_NoArea), num_samples(0), sampler(NULL), area_size(0)
{}

Light::~Light()
{}

void Light::build_area(std::string stype)
{
	half_area_size = area_size/2;


	sampler = new MultiJittered(num_samples);
	if( stype == "noarea")
	{
		num_samples = 1;
		type = LAT_NoArea;
	}
	else if(stype == "square")
	{
		type = LAT_Square;
	}
	else if(stype == "sphere")
	{
		type = LAT_Sphere;
		sampler->map_samples_to_sphere();
	}
	else if(stype == "disk")
	{
		type = LAT_Disk;
		sampler->map_samples_to_unit_disk();
	}
	else if(stype == "hemisphere")
	{
		type = LAT_Hemisphere;
		sampler->map_samples_to_hemisphere(1);
	}
	else
	{
		printf("invalid light type (%s)\n", stype.c_str());
	}
}

Point Light::get_point()
{
	if(num_samples == 1 || type == LAT_NoArea)
		return pos;

	Point dp;
	switch (type)
	{
	case LAT_Square:
		dp = sampler->sample_unit_square();
		break;
	case LAT_Disk:
		dp = sampler->sample_unit_disk();
		break;
	case LAT_Sphere:
		dp = sampler->sample_sphere();
		break;
	case LAT_Hemisphere:
		dp = sampler->sample_hemisphere();
		break;
    case LAT_NoArea:
    default:
        break;
	}
	Point lp = (dp * area_size);
	Point fpos = Point(pos.x + lp.x, pos.y + lp.y, pos.z + lp.z);
	return fpos;
}
