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

#ifndef PPM_H
#define PPM_H
#include "color.h"
#include <string>
#include <vector>

class PPMImage 
{
public:
    PPMImage() : height(0), width(0) {}

	void create(int width, int height);
	void set_color(const int &w, const int &h, Color c);
	bool save(const std::string filename);

    inline bool loaded() {  return height == 0; }
    void load(const std::string &file);

	std::vector< std::vector<Color> > data;
    int height;
    int width;
};

#endif // !PPM_HPP
