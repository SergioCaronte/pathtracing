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

#include <chrono>
#include <fstream>
#include <omp.h>
#include "pathtracer.h"
#include "ppmimage.h"

#define SATURATE(a) std::max(a, 0.0f)

Pathtracer::Pathtracer(int N)
{
    lightsampler = MultiJittered(N);
    lightsampler.map_samples_to_sphere();
}

void Pathtracer::compute(Scene &scene)
{
   compute_path(scene);
}

Vector Pathtracer::get_ray_direction(Scene &scene, int w, int h, const Point &sp)
{
	//pixel position onto view plane
	float ay = -1 * ((float)h -(float)scene.screen.height_px/2 + sp.y);
	float ax = 1 *  ((float)w - (float)scene.screen.width_px/2 + sp.x);
	Vector dir = (ax * scene.camera.x) + (ay * scene.camera.y) - (scene.screen.d * scene.camera.z);
	dir.normalize();
    return dir;
}

Vector Pathtracer::get_ray_direction(Scene &scene, int w, int h, const Point &lp, const Point &sp)
{
    //pixel position onto view plane
	float ay = -1 * ((float)h -(float)scene.screen.height_px/2 + sp.y);
	float ax = 1 *  ((float)w - (float)scene.screen.width_px/2 + sp.x);
	float az = 1;

	az = scene.camera.focal_dist;
	ay = (ay * scene.camera.focal_dist / scene.screen.d);
	ax = (ax * scene.camera.focal_dist / scene.screen.d);

	//Point pt = ori + (ax * scene.camera.x) + (ay * scene.camera.y) - (az * scene.camera.z);
	//Vector dir = pt - ori;
	Vector dir = ((ax - lp.x) * scene.camera.x) + ((ay - lp.y) * scene.camera.y) - (az * scene.camera.z);
	dir.normalize();
	return dir;
}

void Pathtracer::compute_path(Scene &scene)
{
    total_emission = 0;
    for(auto it = scene.luminaries.begin(); it != scene.luminaries.end(); ++it)
    {
        Color e = (*it)->emission;
        total_emission += std::max(e.r, std::max(e.g, e.b));
        (*it)->max_emission = total_emission;
    }
    Screen sc = scene.screen;
    sampler = MultiJittered(sc.samples * scene.N);

    double inv_N = 1.0/scene.N;
    double inv_S = 1.0/sc.samples;

	//output file
	PPMImage output;
	output.create(sc.width_px, sc.height_px);

    // for each pixel of the virtual screen, calculate the color.
    auto start = std::chrono::system_clock::now();
    #pragma omp parallel for schedule(dynamic, 1)
    for(int h = 0; h < sc.height_px; ++h)
	{
		fprintf(stderr,"\rcalculating row %i of %i", h+1, sc.height_px);

		Xi[0] = Xi[1] = 0; Xi[2] = h*h*h;
        for(int w = 0; w < sc.width_px; ++w)
		{
            Color clr;
            // super sampling
            for (int ss = 0; ss < sc.samples; ss++)
            {
                Color p;
                // monte carlo sampling
                for(int n = 0; n < scene.N; n++)
                {
                    Point sp = sampler.sample_unit_square();

                    Point dp = scene.camera.sampler->sample_unit_disk();
                    Point lp = dp * scene.camera.lens_radius;

                    Ray ray;
                    ray.origin = scene.camera.pos + lp.x * scene.camera.x + lp.y * scene.camera.y;
                    // Get a vector width the distance between the camera and the pixel.
                    ray.direction = get_ray_direction(scene, w, h, lp, sp);

                    //Ray ray(scene.camera.pos, get_ray_direction(scene, w, h, sp));
                    p += trace_path(scene, ray, 0) * inv_N;
                }
                p.clamp();
                clr += p * inv_S;
            }
            output.set_color(w, h, clr);
		}
    }
    std::string outname = std::string(scene.output);

    auto duration = std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::system_clock::now() - start);
    long double elapsed = duration.count()/1000000000.0;
    std::cout << "\n" << outname.substr(0,outname.length()-4) << " - elapsed time: " << elapsed << " sec\n";
    std::cout.flush();

    std::stringstream outname_stream;

    if(!output.save(scene.output))
		std::cerr << "Failed to save output file: " << scene.output << std::endl;

}

Color Pathtracer::trace_path(Scene &scene, const Ray &r, int depth, Object *excluded_obj)
{
    Intersection isect;
    if(!intersection(scene, r, isect, NULL, excluded_obj))
        return Color();

    // if it has collided with light, we just return emission
    Color emit = isect.object.emission;
    if(emit.r > 0 || emit.g > 0 || emit.b > 0)  return isect.object.get_color(isect.contact) * emit;

    // Radiance, surface color
    Color L, f;
    // probability distribution function
    double pdf;

    // Direct illumination
    double cos_a_max;
    int light_id;
    Color light_emission;
    Point lp = get_light_pos(scene, isect.contact, light_emission, cos_a_max, light_id);
    Vector shadow_dir = (lp - isect.contact).normalize();
    Intersection direct;
    // Shadow ray
    if(Vector::dot(shadow_dir, isect.normal) > 0 && intersection(scene, Ray(isect.contact, shadow_dir), direct, &lp, &isect.object) && direct.object.id == light_id)
    {
        f = isect.object.BRDF(r.direction, shadow_dir, isect.normal, isect.contact, Xi, pdf);
        //light attenuation
        double omega = 2 * M_PI * ( 1 - cos_a_max);
        double cos_i = SATURATE(Vector::dot(shadow_dir, isect.normal));
        L += (f * (light_emission*cos_i*omega))/pdf;
    }

    // Indirect light
    // not unitary vector, keep it!
    Vector wi = Vector(1,1,1);
    // BRDF returns color as well as gives wi direction and pdf of wi direction
    f = isect.object.BRDF(r.direction, wi, isect.normal, isect.contact, Xi, pdf);
    // Russian roullete based on max(f.intensity) after depth 1
    double c_rr = depth < 1 ? 1 : std::max(f.r, std::max(f.g, f.b));
    if(erand48(Xi) < c_rr)
    {
        // scattering ray
        Ray scat_ray = Ray(isect.contact, wi);
        L += (f * trace_path(scene, scat_ray, depth+1, &isect.object) / (pdf * c_rr) );
    }

    return L;
}

Point Pathtracer::get_light_pos(Scene &scene, const Point &c, Color &emit, double &cos_a_max, int &light_id)
{
    SphereObject* l = scene.luminaries[0];

    float u = erand48(Xi) * total_emission;
    // importance lighting
    for(auto it = scene.luminaries.begin(); it != scene.luminaries.end(); ++it)
        if(u < (*it)->max_emission)
        {
            l = *it;
            break;
        }

    emit = l->emission;
    light_id = l->id;

    cos_a_max = sqrt(1 - l->radius * l->radius / Vector::dot(c-l->original_pos, c-l->original_pos));

    Point p = lightsampler.sample_sphere();
    p *= l->radius/2.0;
    p += l->original_pos;
    return p;
}

bool Pathtracer::intersection(Scene &scene, const Ray &ray,
		Intersection &intersection, const Point *max_pos, const Object *excluded_obj)
{
    // temporarily store the intersections.
    Object *closest_obj = NULL;
    float closest_t = FLT_MAX;
    bool closest_inside = false;
	Vector closest_normal;

    // loop through all objects looking where it intersects.
    for(auto it = scene.objects.begin(); it != scene.objects.end(); ++it)
	{
        // if object is excluded, we are not taking into account
        if(excluded_obj && excluded_obj->id == (*it)->id)
			continue;

		float t;
		Vector n;
		bool insided;
		t = (*it)->hit_test(ray, n, max_pos, &insided);

		//if found closest until now, save info
		if(t > FLT_EPSILON && t < closest_t)
		{
            closest_obj = (*it);
            closest_t = t;
            closest_inside = insided;
			closest_normal = n;
        }
    }

    // if object was found.
    if(closest_obj)
	{
		intersection.contact = ray.origin + closest_t * ray.direction;
		intersection.object = *closest_obj;
		intersection.normal = closest_normal;

		// check if camera is inside hit sphere.
		if(intersection.object.type == ObjectType::Sphere)
		{
			intersection.inside = false;
            if(closest_inside)
			{
				intersection.normal = intersection.normal * (-1);
				intersection.inside = true;
            }
        }
        return true;
    }
    // no intersections.
    return false;
}

Vector Pathtracer::randomHemisphereVectorTowardLight(const SphereObject* light, const Intersection &isect, double & cos_a_max)
{
    Point center = light->original_pos;
    float radius = light->radius;

    Vector sw = center - isect.contact;
    Vector su = Vector::cross((fabs(sw.x) > 0.1 ? Vector(0,1) : Vector(1)), sw).normalize();
    Vector sv = Vector::cross(sw, su);
    // max angle
    cos_a_max = sqrt(1 - radius * radius / Vector::dot(isect.contact - center, isect.contact - center));
    double eps1 = erand48(Xi), eps2 = erand48(Xi);
    double cos_a = 1 - eps1 + eps1 * cos_a_max;
    double sin_a = sqrt(1-cos_a*cos_a);
    double phi = 2 * M_PI * eps2;
    Vector ldir = su * cos(phi) * sin_a + sv * sin(phi) * sin_a + sw * cos_a;
    ldir.normalize();
    return ldir;
}

Vector Pathtracer::randomHemisphereVector(const Vector &n, const Vector &rd)
{
    Vector nl = Vector::dot(n, rd) < 0 ? n : n * (-1);
    // angle around sphere
    double phi = 2 * M_PI * erand48(Xi);
    // distance from center
    double r = erand48(Xi);
    double rs = sqrt(r);
    double x = rs * cos(phi);
    double y = rs * sin(phi);
    double z = sqrt(1 - x*x - y*y);
    // create orthonormal coordinate frame from normal surface
    Vector bw = nl;
    Vector bu = Vector::cross((fabs(bw.x) > 0.1 ? Vector(0,1) : Vector(1)), bw).normalize();
    Vector bv = Vector::cross(bw, bu);
    // random direction inside hemisphere
    return (bu * x + bv * y + bw * z).normalize();
}
