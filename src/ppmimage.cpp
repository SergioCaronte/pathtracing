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

#include "ppmimage.h"
#include <algorithm>
#include <fstream>
#include <sstream>
#include <stdexcept>

void PPMImage::create(int w, int h)
{
	width = w;
	height = h;
	data.resize(height);
	for(size_t i = 0; i < data.size(); ++i)
		data[i].resize(width);
}

void PPMImage::set_color(const int &w, const int &h, Color c)
{
	data[h][w] = c;
}

bool PPMImage::save(const std::string filename)
{
	// Create image PPM file.
    std::ofstream file(filename.c_str());
    if(!file.is_open())
	{
        std::cerr << "Failed to create output file: " << filename << std::endl;
        return false;
    }

    // PPM header.
    file << "P3" << std::endl;
    file << width << " " << height << std::endl;
    file << "255" << std::endl;

	for(int h = 0; h < height; ++h)
        for(int w = 0; w < width; ++w)
			 file << data[h][w];
	file.close();
	return true;
}

void PPMImage::load(const std::string &file)
{
    std::ifstream in(file.c_str());
    if(!in)
    {
        printf("Failed to open ppm file.\n");
        exit(1);
    }

    std::string line;
    std::stringstream ss;
    int max_color;

    getline(in, line);
    if(line != "P6")
    {
        printf("Invalid ppm file (%s).\n", line.c_str());
        exit(1);
    }

    do
    {
        getline(in, line);
    } while(line[0] == '#');

    ss << line;
    getline(in, line);
    ss << " " << line;
    ss >> width >> height >> max_color;

    if(width <= 0 || height <= 0 || max_color <= 0)
    {
        printf("Ppm sizes are invalid, dimesions(%dx%d) max color(%d).\n", width, height, max_color);
        exit(1);
    }
    if(max_color > 255)
    {
        printf("Invalid max color(%d), max must be 255.\n", max_color);
        exit(1);
    }

    data.resize(height);

    // Store all the bytes.
    unsigned char r, g, b;
    Color c;
    for(int i = 0; i < height; ++i)
	{
        for(int j = 0; j < width; ++j)
		{
            in >> r >> g >> b;
            if(!in.good())
                printf("Unformatted ppm.");

            c.r = (float) r / max_color;
            c.g = (float) g / max_color;
            c.b = (float) b / max_color;
            data.at(i).push_back(c);
        }
    }
    in.close();
}
