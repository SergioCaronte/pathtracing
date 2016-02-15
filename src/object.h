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

#ifndef OBJECT_H
#define OBJECT_H

#include "structs.h"
#include "ray.h"
#include "math/matrix.h"

// Types of objects.
enum ObjectType
{
    Sphere,
    Polyhedron,
	Torus,
	Cylinder
};

// Base Object.
class Object
{
public:
	Object(){}

	virtual void calculate_matrices(float dt = 0) {}
	virtual float hit_test(const Ray &ray, Vector &normal, const Point *endPos = NULL, bool *inside = NULL){return -1.0;}
	Vector mul_vec( Matrix4x4 &mat, const Vector &p);
	Point mul_point( Matrix4x4 &mat, const Point &p);
    Color get_color(const Point &p);
    Color BRDF(const Vector &wo, Vector &wi, const Vector &n, const Point &p, unsigned short *Xi, double &pdf);
    Color sample_blinnphong(const Vector &wo, Vector &wi, const Vector &n, const Point &p, unsigned short *Xi,double &pdf);
    Color sample_lamb(const Vector &wo, Vector &wi, const Vector &n, const Point &p, unsigned short *Xi, double &pdf);
    Color sample_perfect_specular(const Vector &wo, Vector &wi, const Vector &n, const Point &p, unsigned short *Xi, double &pdf);
    Color sample_transmission(const Vector &wo, Vector &wi, const Vector &n, const Point &p, unsigned short *Xi, double &pdf);
    // Torrance-Sparrow
    Color sample_torrancesparrow(const Vector &wo, Vector &wi, const Vector &n, const Point &p, unsigned short *Xi, double &pdf);
    // Fresnel factor for dieletric
    double FDiel(double cos);
    // Fresnel factor for conductors
    double FCond(double cos, double k);
    // Geometric factor
    double G(const Vector &wo, const Vector &wi, const Vector &wh, const Vector &n);
    // Blinn distribuction
    double D(const Vector &wh, const Vector &n);

    Vector randomHemisphereVector(const Vector &n, const Vector &rd, unsigned short * Xi);
    Vector randomImportanceSampling(const Vector &n, const Vector &rd, unsigned short * Xi);
    Vector randomImportanceSamplingHalfway(const Vector &n, unsigned short * Xi);
    inline float AbsCosTheta(const Vector &w) { return fabsf(w.z); }

	size_t id;
    ObjectType type;
	Texture texture;
	Material material;

	Color emission;

	Vector acceleration;
	Point pos;
	// affines transforms
	Point original_pos;
	Point original_rot;
    Point original_scale;
	Matrix4x4 inv_trans;
};


// Sphere object.
class SphereObject : public Object
{
public:
	SphereObject();

	void calculate_matrices(float dt = 0) override;

	float hit_test(const Ray &ray, Vector &normal, const Point *max_pos = NULL, bool *inside = NULL) override;
	// sphere radius.
    float radius;
    //
    float max_emission;
};

// Polyhedron object.
class PolyhedronObject : public Object
{
public:
	PolyhedronObject();

	void calculate_matrices(float dt = 0) override;

	float hit_test(const Ray &ray, Vector &normal, const Point *max_pos = NULL, bool *inside = NULL) override;
    // number of faces.
    int numFaces;
    // faces
    std::vector<Plane> planes;
};

class TorusObject : public Object
{
public:
	TorusObject();

	void calculate_matrices(float dt = 0) override;

	float hit_test(const Ray &ray, Vector &normal, const Point *max_pos = NULL, bool *inside = NULL) override;
	// torus radius
	double radius;
	// torus thickness
	double thickness;

private:
	Vector compute_normal(const Point& p);
};

class CylinderObject : public Object
{
public:
	CylinderObject();

	void calculate_matrices(float dt = 0) override;

	float hit_test(const Ray &ray, Vector &normal, const Point *max_pos = NULL, bool *inside = NULL) override;
	// y position of base
	double bottom;
	// y position of top
	double top;
	// cylinder radius
	double radius;

private:
	Point bottom_pos;
	Point up_pos;
	double inv_radius;

	float hit_cylinder(const Ray &ray, Vector &normal, bool &inside);
	float hit_disk(const Ray &ray, Vector &normal, Vector disk_normal, Point disk_pos);

};

#endif
