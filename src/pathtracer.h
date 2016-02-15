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

#ifndef PATHTRACER_H
#define PATHTRACER_H

#include <random>

#include "color.h"
#include "ray.h"
#include "math/point.h"
#include "math/vector.h"
#include "structs.h"
#include "scene.h"
#include "intersection.h"
#include "multijittered.h"

class Pathtracer
{
public:
	Pathtracer(int max_ray_travel = 4);
	~Pathtracer(){}

	// compute raytracing. trace a ray for every pixel
	void compute(Scene &scene) ;
	// trace the ray path, raytracing core
	Color trace(Scene &scene, const Ray &ray, size_t depth, Object *excluded_obj = NULL);

    bool intersection(Scene &scene, const Ray &ray, Intersection &interc,
		const Point *max_pos = NULL, const Object *excluded_obj = NULL);

private:
	MultiJittered sampler;
	MultiJittered lightsampler;     // sampler for light position

	void compute_regular(Scene &scene);
	void compute_sampled(Scene &scene);
	inline Vector get_ray_direction(Scene &scene, int w, int h, const Point &lp, const Point &sp);
	inline Vector get_ray_direction(Scene &scene, int w, int h);
	inline Vector get_ray_direction(Scene &scene, int w, int h, const Point &sp);

    // path tracer methods
    void compute_path(Scene &scene);
    // smallpt
    Color radiance(Scene &scene, const Ray &r, int depth, unsigned short *Xi, Object *excluded_obj = NULL);

    Color trace_path(Scene &scene, const Ray &r, int depth, Object *excluded_obj = NULL);
    Point get_light_pos(Scene &scene, const Point &c, Color &emit, double &cos_a_max, int & light_id);
    // used onto random uniform unit generation
    unsigned short Xi[3];

    // return a unit direction vector around an hemisphere of n.
    // @param n, surface normal
    // @param rd, ray direction
    Vector randomHemisphereVector(const Vector &n, const Vector &rd);
    Vector randomHemisphereVectorTowardLight(const SphereObject* light, const Intersection &isect, double & cos_a_max);

    float total_emission;

};
#endif
