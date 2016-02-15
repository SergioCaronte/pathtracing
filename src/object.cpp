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

#include "object.h"
#include "scene.h"
#include "math/math.h"
#include "math/matrix.h"
#include "pathtracer.h"
#include <random>

#define SATURATE(a) std::max(a, 0.0f)
#define	IsZero(x)	((x) > -EQN_EPS && (x) < EQN_EPS)

//////////////////////////////////////////////////////////
/// Object class
//////////////////////////////////////////////////////////
Color Object::get_color(const Point &p)
{
	Color clr;
	switch(texture.type)
	{
        case SolidTexture:
            clr = texture.solid.color;
		break;
        case CheckerTexture:
			int pattern;
            pattern =	floor(p.x / texture.checker.scale) +
						floor(p.y / texture.checker.scale) +
						floor(p.z / texture.checker.scale);
            clr = (pattern % 2) ? texture.checker.color2 : texture.checker.color1;
		break;
        case MapTexture:
			double s, r;
			int u, v;
            //texture interpolation
			s =	texture.map.p0.x * p.x +
                texture.map.p0.y * p.y +
                texture.map.p0.z * p.z +
                texture.map.p0.w;
            r =	texture.map.p1.x * p.x +
                texture.map.p1.y * p.y +
                texture.map.p1.z * p.z +
                texture.map.p1.w;
            u = (int)(r * texture.map.image.height) % texture.map.image.height;
            v = (int)(s * texture.map.image.width) % texture.map.image.width;
            if(u < 0) u += texture.map.image.height;
            if(v < 0) v += texture.map.image.width;
            clr = texture.map.image.data.at(u).at(v);
		break;
    }
	return clr;
}

Color Object::BRDF(const Vector &wo, Vector &wi, const Vector &n, const Point &p, unsigned short *Xi, double &pdf)
{
    double u = erand48(Xi);
    double ka = material.kA;
    double kd = ka + material.kD;
    double ks = kd + material.kS;
    double kr = ks + material.kR;
    double kt = kr + material.kT;

    if(u < ka)
        return sample_torrancesparrow(wo, wi, n, p, Xi, pdf);
    else if(ka <= u && u < kd)
        return sample_lamb(wo, wi, n, p, Xi, pdf);
    else if(kd <= u && u < ks)
        return sample_blinnphong(wo, wi, n, p, Xi, pdf);
    else if(ks <= u && u < kr)
        return sample_perfect_specular(wo, wi, n, p, Xi, pdf);
    else if(kr <= u && u < kt)
        return sample_transmission(wo, wi, n, p, Xi, pdf);
    else
    {   pdf = 1; return Color(); }
}

Color Object::sample_lamb(const Vector &wo, Vector &wi, const Vector &n, const Point &p, unsigned short *Xi, double &pdf)
{
    if(wi.length() > 1.1)
        wi = randomHemisphereVector(n, wo, Xi);

    double cos_i = SATURATE(Vector::dot(wi, n));
    pdf = cos_i * M_1_PI;
    if(pdf <= 0.0)
    {
        wi = randomHemisphereVector(n, wo, Xi);
        return sample_lamb(wo, wi, n, p, Xi, pdf);
    }
    // albedo * cdiff/PI * cos_i
    return get_color(p) * material.kD * M_1_PI * cos_i;
}

Color Object::sample_blinnphong(const Vector &wo, Vector &wi, const Vector &n, const Point &p, unsigned short *Xi, double &pdf)
{
    if(wi.length() > 1.1)
        wi = randomImportanceSampling(n, wo, Xi);
    else // bug fix for direct illumination
        return sample_lamb(wo, wi, n, p, Xi, pdf);
    // halfway vector
    Vector eye = (wo * (-1)).normalize();
    Vector wh = (wi + eye).normalize();
    double cos_h = SATURATE(Vector::dot(wh, n));
    // pdf
    pdf = ((material.alpha + 2.0) / (2.0 * M_PI)) * pow(cos_h, material.alpha);
    // if pdf 0, resampling
    if(pdf <= 0.0)  return sample_blinnphong(wo, wi, n, p, Xi, pdf);
    // albedo * m+8/8PI * cspec * (h . n)^m
    return get_color(p) * (((material.alpha+8)/(8*M_PI)) * material.kS * pow(cos_h, material.alpha));
}

Color Object::sample_perfect_specular(const Vector &wo, Vector &wi, const Vector &n, const Point &p, unsigned short *Xi, double &pdf)
{
    pdf = 1;
    // do not calculante mirror on direct illumination
    if(wi.length() < 1.1)   return Color();

    wi = (wo - n * 2 * Vector::dot(n, wo)).normalize();
    return get_color(p) * material.kR;
}

Color Object::sample_transmission(const Vector &wo, Vector &wi, const Vector &n, const Point &p, unsigned short *Xi, double &pdf)
{
    pdf = 1;
    // do not calculate refraction on direct illumination
    if(wi.length() < 1.1)   return Color();

    Vector nl = Vector::dot(n, wo) < 0 ? n : n * (-1);
    // ray from outside going in
    bool into = Vector::dot(n, nl) > 0;
    // ratio transmition factor
    double eta = into ? 1.0/material.ior : material.ior;
    double cos_o = Vector::dot(wo, nl);
    double cos2t = 1 - (eta * eta) * (1 - cos_o * cos_o);

    // total internal reflection
    if( cos2t < 0 )
        wi = (wo - nl * 2 * Vector::dot(nl, wo)).normalize();
    else
        wi = (wo / eta - ( sqrt(cos2t) - cos_o / eta ) * nl ).normalize();
    return get_color(p) * material.kT;
}

Color Object::sample_torrancesparrow(const Vector &wo, Vector &wi, const Vector &n, const Point &p, unsigned short *Xi, double &pdf)
{
    // inverting wo
    Vector eye = (wo * (-1)).normalize();
    Vector wh;
    // if not direct light
    if(wi.length() > 1.1)
    {
        // importance sampling
        wh = randomImportanceSamplingHalfway(n, Xi);
        // if wh is not in the same hemisphere, turn it;
        if(Vector::dot(eye, wh) <= 0.0) wh = wh * (-1);
        // reflecting wo based on wh
        wi = wo + 2 * Vector::dot(wh, eye) * wh;
    }
    else
    {
        // do not calculate torrancesparrow on direct illumination
        pdf = 1;
        return Color();
    }

    // cos and plus cos
    double cos_o = SATURATE(Vector::dot(n, eye));
    double cos_i = SATURATE(Vector::dot(n, wi));
    double cos_alpha = SATURATE(Vector::dot(eye, wh));
    double cos_h = SATURATE(Vector::dot(n, wh));

    if(cos_i == 0.f || cos_o == 0.f)
    {
        pdf = 1;
        return Color();
    }

    pdf = (((material.alpha+1.0) * pow(cos_h, material.alpha)) / (2.0 * M_PI * 4.f * cos_alpha));
    // if pdf 0, resampling
    if(pdf <= 0.0)  return sample_torrancesparrow(wo, wi, n, p, Xi, pdf);

    Color F;
    F.r = FCond(cos_alpha, material.absortion.r);
    F.g = FCond(cos_alpha, material.absortion.g);
    F.b = FCond(cos_alpha, material.absortion.b);
    return (get_color(p) * (F * ((D(wh, n) * G(eye, wi, wh, n) / (4.f * cos_i * cos_o))))) * material.kA;
}

double Object::FDiel(double cos)
{
    Material m = material;
    double a = m.ior - 1;
    double b = m.ior + 1;
    // fresnel reflectance factor
    double R0 = a * a/ (b * b);
    double c = 1 - cos;
    // Schlick approximation
    return R0 + (1-R0) * c *c * c * c * c;
}

double Object::FCond(double cos, double k)
{
    // Extract from Physically Based Rendering
    double eta = material.ior;

    double tmp_f = eta*eta + k*k;
    double tmp = tmp_f * cos * cos;
    double Rparl2 = (tmp - (2.0 * eta * cos) + 1 ) /
                    (tmp + (2.0 * eta * cos) + 1 );

    double Rperp2 = (tmp_f - (2.0 * eta * cos) + cos * cos) /
                    (tmp_f + (2.0 * eta * cos) + cos * cos);

    return (Rparl2 + Rperp2) / 2.0;
}

double Object::G(const Vector &wo, const Vector &wi, const Vector &wh, const Vector &n)
{
    double NdotWH = SATURATE(Vector::dot(n, wh));
    double NdotWO = SATURATE(Vector::dot(n, wo));
    double NdotWI = SATURATE(Vector::dot(n, wi));
    double WOdotWh = SATURATE(Vector::dot(wo, wh));

    return std::min(1.0, std::min((2.0 * NdotWH * NdotWO / WOdotWh), (2.0 * NdotWH * NdotWI / WOdotWh)));
}

double Object::D(const Vector &wh, const Vector &n)
{
    // expoent
    double e = material.alpha;
    // cos of normal and half vector
    double cos_h = Vector::dot(n, wh);
    // Normalized Blinn Dist
    // (exp+2 / 2PI) * cos_h^exp
    return ((e+2.0)/(2.0*M_PI)) * pow(cos_h, e);
}

Vector Object::randomHemisphereVector(const Vector &n, const Vector &rd, unsigned short * Xi)
{
    Vector nl = Vector::dot(n, rd) < 0 ? n : n * (-1);

    double phi = 2 * M_PI * erand48(Xi);
    double theta = acos(sqrt(erand48(Xi)));
    double x = sin(theta) * cos(phi);
    double y = sin(theta) * sin(phi);
    double z = cos(theta);
    // create orthonormal coordinate frame from normal surface
    Vector bw = nl.normalize();
    Vector bu = Vector::cross((fabs(bw.x) > 0.1 ? Vector(0,1) : Vector(1)), bw).normalize();
    Vector bv = Vector::cross(bw, bu).normalize();
    // random direction inside hemisphere
    return (bu * x + bv * y + bw * z).normalize();
}

Vector Object::randomImportanceSampling(const Vector &n, const Vector &rd, unsigned short * Xi)
{
    double cos_t = pow(erand48(Xi), 1.0/ (material.alpha + 1));
    double sin_t = sqrt(std::max(0.0, 1.0 - cos_t*cos_t));
    double phi = 2 * M_PI * erand48(Xi);
    double x = sin_t * cos(phi);
    double y = sin_t * sin(phi);
    double z = cos_t;

    // create orthonormal coordinate frame from reflected ray
    Vector bw = (rd - n * 2 * Vector::dot(n, rd)).normalize();
    Vector bu = Vector::cross((fabs(bw.x) > 0.1 ? Vector(0,1) : Vector(1)), bw).normalize();
    Vector bv = Vector::cross(bw, bu);
    // random direction inside hemisphere
    Vector dir = (bu * x + bv * y + bw * z).normalize();
    // if inside surface, sample again
    if( Vector::dot(n, dir) <= FLT_EPSILON)
        return randomImportanceSampling(n, rd, Xi);
    return dir;
}

Vector Object::randomImportanceSamplingHalfway(const Vector &n, unsigned short * Xi)
{
    double cos_t = pow(erand48(Xi), 1.0/ (material.alpha + 1));
    double sin_t = sqrt(std::max(0.0, 1.0 - cos_t*cos_t));
    double phi = 2 * M_PI * erand48(Xi);
    double x = sin_t * cos(phi);
    double y = sin_t * sin(phi);
    double z = cos_t;

    // create orthonormal coordinate frame from normal
    Vector bw = n;
    Vector bu = Vector::cross((fabs(bw.x) > 0.1 ? Vector(0,1) : Vector(1)), bw).normalize();
    Vector bv = Vector::cross(bw, bu);

    return (bu * x + bv * y + bw * z).normalize();
}

Vector Object::mul_vec( Matrix4x4 &mat, const Vector &v)
{
	return Vector(	mat[0][0] * v.x + mat[0][1] * v.y + mat[0][2] * v.z,
					mat[1][0] * v.x + mat[1][1] * v.y + mat[1][2] * v.z,
					mat[2][0] * v.x + mat[2][1] * v.y + mat[2][2] * v.z);
}

Point Object::mul_point( Matrix4x4 &mat, const Point &p)
{
	return Point(	mat[0][0] * p.x + mat[0][1] * p.y + mat[0][2] * p.z + mat[0][3],
					mat[1][0] * p.x + mat[1][1] * p.y + mat[1][2] * p.z + mat[1][3],
					mat[2][0] * p.x + mat[2][1] * p.y + mat[2][2] * p.z + mat[2][3]);
}

//////////////////////////////////////////////////////////
/// Sphere class
//////////////////////////////////////////////////////////

SphereObject::SphereObject()
{
	type = Sphere;
}

void SphereObject::calculate_matrices(float dt)
{
	// acceleration at time dt
	Vector c_acc = acceleration * dt;
	// current pos
	pos.x = original_pos.x + c_acc.x;
	pos.y = original_pos.y + c_acc.y;
	pos.z = original_pos.z + c_acc.z;

    Matrix4x4 i_t;
	i_t.identity();
	i_t[0][3] = -pos.x;	i_t[1][3] = -pos.y;	i_t[2][3] = -pos.z;

	//inverse scale matrix
	Matrix4x4 i_s;
	i_s.identity();
	i_s[0][0] = 1/original_scale.x;	i_s[1][1] = 1/original_scale.y;	i_s[2][2] = 1/original_scale.z;

	//object inverse transform matrix
	inv_trans = i_s;
}

float SphereObject::hit_test(const Ray &ray, Vector &normal, const Point *max_pos, bool *inside)
{
	Ray inv_ray;
	inv_ray.origin = mul_point(inv_trans, ray.origin);
	inv_ray.direction = mul_vec(inv_trans, ray.direction);
	inv_ray.direction.normalize();

	Vector e = (inv_ray.origin - pos);
    double a = Vector::dot(inv_ray.direction, inv_ray.direction);
    double b = 2.0f * Vector::dot(inv_ray.direction, e);
    double c = Vector::dot(e, e) - std::pow(radius, 2);
    double delta = std::pow(b, 2) - 4*a*c;

	if(delta < FLT_EPSILON)
        return -1.0f;

    // Resolving the equation.
    delta = std::sqrt(delta);

    // Test both solutions.
    double t1 = (-b - delta) / (2.0f * a);
    double t2 = (-b + delta) / (2.0f * a);

    // Calculate the maximum t.
    double max_t;
    if(max_pos)
        max_t = (max_pos->x - inv_ray.origin.x) / inv_ray.direction.x;
    else
        max_t = FLT_MAX;

    // t1 is always smaller than t2 because it uses -b -discriminant.
    if(t1 > FLT_EPSILON && t1 < max_t)
	{
        if(inside)	*inside = false;
		Point intersection = inv_ray.origin + t1 * inv_ray.direction;
		normal = (intersection - pos).normalize();
		normal = mul_vec(inv_trans, normal);
		normal.normalize();
        return t1;
    }
    else if(t2 > FLT_EPSILON && t2 < max_t)
	{
        if(inside)	*inside = true;
		Point intersection = inv_ray.origin + t2* inv_ray.direction;
		normal = (intersection - pos).normalize();
		normal = mul_vec(inv_trans, normal);
		normal.normalize();
        return t2;
    }
    else
        return -1.0f;
}

//////////////////////////////////////////////////////////
/// Polyhedron class
//////////////////////////////////////////////////////////

PolyhedronObject::PolyhedronObject()
{
	type = Polyhedron;
}

void PolyhedronObject::calculate_matrices(float dt)
{}

float PolyhedronObject::hit_test(const Ray &ray, Vector &normal, const Point *max_pos, bool *inside)
{
    float t0 = 0.0f, t1 = FLT_MAX;
    Vector n_t0, n_t1;

	// run for every face
    for(size_t i = 0; i < planes.size(); ++i)
	{
        Vector p0(ray.origin.x, ray.origin.y, ray.origin.z);
        Vector n = planes[i].normal();
        float dn = Vector::dot(ray.direction, n);
        float vd = Vector::dot(p0, n) + planes[i].d;

        if((dn <= FLT_EPSILON && dn >= -FLT_EPSILON) && vd > FLT_EPSILON)
		{
			t1 = -1.0f;
        }
        else if(dn > FLT_EPSILON && t1 > -vd / dn)
		{
			 t1 = -vd / dn;
             n_t1 = n;
        }
        else if(dn < -FLT_EPSILON && t0 < -vd / dn)
		{
            t0 = -vd / dn;
            n_t0 = n;
        }
    }

    // Calculate the maximum t.
    double max_t;
    if(max_pos)
        max_t = (max_pos->x - ray.origin.x) / ray.direction.x;
    else
        max_t = FLT_MAX;

    if(t1 < t0)
	{
		return -1.0f;
	}
    else if(std::abs(t0) <= FLT_EPSILON && (t1 >= t0) && t1 < max_t)
	{
        normal = (n_t1 * (-1)).normalize();
        return t1;
    }
    else if(t0 > FLT_EPSILON && t1 >= t0 && t0 < max_t)
	{
        normal = n_t0.normalize();
        return t0;
    }
    return -1.0f;
}

/////////////////////////////////////////
/// Torus
/////////////////////////////////////////
TorusObject::TorusObject()
{
	type = Torus;
}

void TorusObject::calculate_matrices(float dt)
{
	// acceleration at time dt
	Vector c_acc = acceleration * dt;
	// current pos
	pos.x = original_pos.x + c_acc.x;
	pos.y = original_pos.y + c_acc.y;
	pos.z = original_pos.z + c_acc.z;

	Point rot;
	rot.x = toRads(original_rot.x);
	rot.y = toRads(original_rot.y);
	rot.z = toRads(original_rot.z);

	//inverse translate matrix
	Matrix4x4 i_t;
	i_t.identity();
	i_t[0][3] = -pos.x;	i_t[1][3] = -pos.y;	i_t[2][3] = -pos.z;

	//inverse rotation matrix
	Matrix4x4 i_rx, i_ry, i_rz;
	i_rx.identity();
	i_rx[1][1] = cos(rot.x);		i_rx[1][2] = sin(rot.x);
	i_rx[2][1] = -sin(rot.x);		i_rx[2][2] = cos(rot.x);

	i_ry.identity();
	i_ry[0][0] = cos(rot.y);		i_ry[0][2] = -sin(rot.y);
	i_ry[2][0] = sin(rot.y);		i_ry[2][2] = cos(rot.y);

	i_rz.identity();
	i_rz[0][0] = cos(rot.z);		i_rz[0][1] = sin(rot.z);
	i_rz[1][0] = -sin(rot.z);		i_rz[1][1] = cos(rot.z);

	//inverse scale matrix
	Matrix4x4 i_s;
	i_s.identity();
	i_s[0][0] = 1/original_scale.x;	i_s[1][1] = 1/original_scale.y;	i_s[2][2] = 1/original_scale.z;

	//object inverse transform matrix
	inv_trans = i_s * i_ry * i_rx * i_rz * i_t;
}

float TorusObject::hit_test(const Ray &ray, Vector &normal, const Point *max_pos, bool *inside)
{
	Ray inv_ray;
	inv_ray.origin = mul_point(inv_trans, ray.origin);
	inv_ray.direction = mul_vec(inv_trans, ray.direction);
	inv_ray.direction.normalize();

	double x1 = inv_ray.origin.x; double y1 = inv_ray.origin.y; double z1 = inv_ray.origin.z;
	double d1 = inv_ray.direction.x; double d2 = inv_ray.direction.y; double d3 = inv_ray.direction.z;

	double coeffs[5];	// coefficient array
	double roots[4];	// solution array

	//define the coefficients
	double sum_d_sqrd = d1*d1 + d2*d2 + d3*d3;
	double e = x1*x1 + y1*y1 + z1*z1 - radius*radius - thickness*thickness;
	double f = x1*d1 + y1*d2 + z1*d3;
	double four_a_sqrd = 4.0 * radius*radius;

	coeffs[0] = e*e - four_a_sqrd * (thickness*thickness-y1*y1);	// constante term
	coeffs[1] = 4.0 * f * e + 2.0 * four_a_sqrd * y1 *d2;
	coeffs[2] = 2.0 * sum_d_sqrd * e + 4.0 * f * f + four_a_sqrd * d2 * d2;
	coeffs[3] = 4.0 * sum_d_sqrd * f;
	coeffs[4] = sum_d_sqrd * sum_d_sqrd;	// coefficient of t^4

	//fin the roots
	int num_real_roots = solveQuartic(coeffs, roots);

	bool intersected = false;
	double t = FLT_MAX;

	if ( num_real_roots == 0)	// ray misses the torus
		return -1.0;

	//find the smallest root greater than FLT_EPSILON, if any
	for( int j = 0; j < num_real_roots; j++)
	{
		if(roots[j] > FLT_EPSILON)
		{
			intersected = true;
			if(roots[j] < t)
				t = roots[j];
		}
	}

	double max_t;
    if(max_pos)
        max_t = (max_pos->x - ray.origin.x) / ray.direction.x;
    else
        max_t = FLT_MAX;

	if( t > max_t || !intersected)
		return -1.0;


	Point hit = inv_ray.origin + t*inv_ray.direction;
	normal = compute_normal(hit);
	normal = mul_vec(inv_trans, normal);
	normal.normalize();
	return t;
}

Vector TorusObject::compute_normal(const Point& p)
{
	Vector normal;
	double param_squared = radius * radius + thickness * thickness;

	double x = p.x;
	double y = p.y;
	double z = p.z;
	double sum_squared = x * x + y * y + z * z;

	normal.x = 4.0 * x * (sum_squared - param_squared);
	normal.y = 4.0 * y * (sum_squared - param_squared + 2.0 * radius * radius);
	normal.z = 4.0 * z * (sum_squared - param_squared);
	normal.normalize();

	return (normal);
}

/////////////////////////////////////////////////////////
/// Cylinder
/////////////////////////////////////////////////////////
CylinderObject::CylinderObject()
{
	type = Cylinder;
}

void CylinderObject::calculate_matrices(float dt)
{
	inv_radius = 1.0/radius;

	// acceleration at time dt
	Vector c_acc = acceleration * dt;
	// current pos
	pos.x = original_pos.x + c_acc.x;
	pos.y = original_pos.y + c_acc.y;
	pos.z = original_pos.z + c_acc.z;

	Point rot;
	rot.x = toRads(original_rot.x);
	rot.y = toRads(original_rot.y);
	rot.z = toRads(original_rot.z);

	//inverse translate matrix
	Matrix4x4  i_t;
	i_t.identity();
	i_t[0][3] = -pos.x;	i_t[1][3] = -pos.y;	i_t[2][3] = -pos.z;

	//inverse rotation matrix
	Matrix4x4 i_rx, i_ry, i_rz;
	i_rx.identity();
	i_rx[1][1] = cos(rot.x);		i_rx[1][2] = sin(rot.x);
	i_rx[2][1] = -sin(rot.x);		i_rx[2][2] = cos(rot.x);

	i_ry.identity();
	i_ry[0][0] = cos(rot.y);		i_ry[0][2] = -sin(rot.y);
	i_ry[2][0] = sin(rot.y);		i_ry[2][2] = cos(rot.y);

	i_rz.identity();
	i_rz[0][0] = cos(rot.z);		i_rz[0][1] = sin(rot.z);
	i_rz[1][0] = -sin(rot.z);		i_rz[1][1] = cos(rot.z);

	//inverse scale matrix
	Matrix4x4 i_s;
	i_s.identity();
	i_s[0][0] = 1/original_scale.x;	i_s[1][1] = 1/original_scale.y;	i_s[2][2] = 1/original_scale.z;

	//object inverse transform matrix
	inv_trans = i_s * i_ry * i_rx * i_rz * i_t;
}

float CylinderObject::hit_test(const Ray &ray, Vector &normal, const Point *max_pos, bool *inside)
{
	Ray inv_ray;
	inv_ray.origin = mul_point(inv_trans, ray.origin);
	inv_ray.direction = mul_vec(inv_trans, ray.direction);
	inv_ray.direction.normalize();

	Vector c_normal;
	bool c_inside;
	double t = FLT_MAX;

	// check collision with cylinder body
	double c_t = hit_cylinder(inv_ray, c_normal, c_inside);
	if(c_t > FLT_EPSILON && c_t < t)
	{
		t = c_t;
		normal = c_normal;
		if(inside)
			*inside = c_inside;
	}

	// check collision with bottom disk
	c_t = hit_disk(inv_ray, c_normal, Vector(0, -1, 0), Point(0, bottom, 0));
	if(c_t > FLT_EPSILON && c_t < t)
	{
		t = c_t;
		normal = c_normal;
		if(inside)
			*inside = false;
	}

	// check collision with top disk
	c_t = hit_disk(inv_ray, c_normal, Vector(0, 1, 0), Point(0, top, 0));
	if(c_t > FLT_EPSILON && c_t < t)
	{
		t = c_t;
		normal = c_normal;
		if(inside)
			*inside = false;
	}

	double max_t;
    if(max_pos)
        max_t = (max_pos->x - ray.origin.x) / ray.direction.x;
    else
        max_t = FLT_MAX;

	if(t < max_t)
	{
		normal = mul_vec(inv_trans, normal);
		normal.normalize();
		return t;
	}
	else
		return -1.0;
}

float CylinderObject::hit_cylinder(const Ray &ray, Vector &normal, bool &inside)
{
	double t;
	double ox = ray.origin.x;	double oy = ray.origin.y;	double oz = ray.origin.z;
	double dx = ray.direction.x;	double dy = ray.direction.y;	double dz = ray.direction.z;

	double a = dx * dx + dz * dz;
	double b = 2.0 * (ox * dx + oz * dz);
	double c = ox * ox + oz * oz - radius * radius;
	double disc = b * b - 4.0 * a * c ;

	Vector inv_rdir = ray.direction * (-1);

	if (disc < 0.0)
		return -1.0;
	else
	{
		double e = sqrt(disc);
		double denom = 2.0 * a;
		t = (-b - e) / denom;    // smaller root

		if (t > FLT_EPSILON)
		{
			double yhit = oy + t * dy;

			if (yhit > bottom && yhit < top)
			{
				normal = Vector((ox + t * dx) * inv_radius, 0.0, (oz + t * dz) * inv_radius);
				inside = false;
				// test for hitting from inside
				if ( Vector::dot(inv_rdir, normal) < 0.0)
				{
					normal = normal * (-1);
					inside = true;
				}
				return t;
			}
		}

		t = (-b + e) / denom;    // larger root
		if (t > FLT_EPSILON)
		{
			double yhit = oy + t * dy;

			if (yhit > bottom && yhit < top)
			{
				normal = Vector((ox + t * dx) * inv_radius, 0.0, (oz + t * dz) * inv_radius);
				inside = false;
				// test for hitting inside surface
				if ( Vector::dot(inv_rdir, normal) < 0.0)
				{
					normal = normal * (-1);
					inside = true;
				}
				return t;
			}
		}
	}
	return -1.0;
}

float CylinderObject::hit_disk(const Ray &ray, Vector &normal, Vector disk_normal, Point disk_pos)
{
	double r_squared = radius * radius;

	double t = Vector::dot((disk_pos - ray.origin), disk_normal) / Vector::dot(ray.direction, disk_normal);

	if (t <= FLT_EPSILON)
		return -1.0;

	Point p = ray.origin + t * ray.direction;
	double d_squared =  pow(disk_pos.x - p.x, 2) +
						pow(disk_pos.y - p.y, 2) +
						pow(disk_pos.z - p.z, 2);

	if (d_squared < r_squared)
	{
		normal = disk_normal;
		return t;
	}
	return -1.0;
}


