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

#include "scene.h"
#include "pathtracer.h"
#include <fstream>
#include <iostream>
#include <cstdlib>

void Scene::load_file(int argc, char **argv)
{
	if(argc < 4)
	{
		std::cerr << "Cmd line usage: " << argv[0] << " input output number_ray_sample [width] [height]\n";
        exit(1);
	}

	// input file name
	input = argv[1];
	// output file name
    output = argv[2];
    // number of ray sample per pixel
    N = atoi(argv[3]);

	// check if resolution parameters were set
    if(argc == 6)
	{
        screen.width_px = strtoul(argv[4], NULL, 10);
        screen.height_px = strtoul(argv[5], NULL, 10);
    }
    else
	{
        screen.width_px = 800;
        screen.height_px = 600;
    }

    // open the input file.
    std::ifstream f_input(input);
    if(!f_input.is_open())
    {
        printf("Failed to open input file (%s).\n", input);
        exit(1);
    }

    // parse the input file.
    parse_camera(f_input);
    parse_light(f_input);
    parse_texture(f_input);
    parse_material(f_input);
    parse_object(f_input);

    // close the input file.
    f_input.close();
}

void Scene::calculate_cam_base()
{
	Vector view = camera.lookat - camera.pos;
	view.normalize();
	camera.z = view * (-1);
	camera.x = Vector::cross(view, camera.up);
	camera.x.normalize();
	camera.y = Vector::cross(camera.z, camera.x);
	camera.y.normalize();
}

void Scene::parse_camera(std::ifstream &in)
{
    in >> camera.pos.x >> camera.pos.y >> camera.pos.z;

    in >> camera.lookat.x >> camera.lookat.y >> camera.lookat.z;
    camera.dir = camera.lookat - camera.pos;
    camera.dir.normalize();

    in >> camera.up.x >> camera.up.y >> camera.up.z;
    camera.up.normalize();

	in >> camera.fovy >> screen.samples >> camera.lens_radius >> camera.focal_dist;
	in >> camera.shutter_time >> camera.exposure;

	calculate_cam_base();

	// screen features
	// calculate viewplane distance
	screen.d = 0.5f * screen.height_px  / tan((M_PI/180.0)*0.5f*camera.fovy);

    screen.theta = M_PI * camera.fovy/180.0;
	screen.aspect = float(screen.width_px)/float(screen.height_px);

	// add viewplane distance to focal_distance
	camera.focal_dist += Point::distance(camera.pos, camera.lookat);

	// start camera sampler
	camera.sampler = new MultiJittered(screen.samples * N);
	camera.sampler->map_samples_to_unit_disk();
}

void Scene::parse_light(std::ifstream &in)
{
    Light light;

    in >> num_lights;
    if(in.fail())
    {
        printf("Failed to read the input file in light description.\n");
        exit(1);
    }
    // ambient light.
    in >> ambient.pos.x >> ambient.pos.y >> ambient.pos.z;
    in >> ambient.color.r >> ambient.color.g >> ambient.color.b;
    in >> ambient.att.a >> ambient.att.b >> ambient.att.c;

    // other lights.
	std::string type;
    for(int i = 1; i < num_lights; ++i)
	{
        in >> light.pos.x >> light.pos.y >> light.pos.z;
        in >> light.color.r >> light.color.g >> light.color.b;
        in >> light.att.a >> light.att.b >> light.att.c;
		in >> type;
		in >> light.num_samples >> light.area_size;

		light.build_area(type);
        lights.push_back(light);
    }
}

void Scene::parse_texture(std::ifstream &in)
{
    std::string type;
    Texture tex;

    in >> num_textures;
    if(in.fail())
    {
        printf("Failed to read the input file in texture description.\n");
        exit(1);
    }
    for(int i = 0; i < num_textures; ++i)
	{
        tex.id = i;
        in >> type;

        if(type == "solid")
		{
            tex.type = SolidTexture;
            in >> tex.solid.color.r >> tex.solid.color.g >> tex.solid.color.b;
        }
        else if(type == "checker")
		{
            tex.type = CheckerTexture;
            in >> tex.checker.color1.r >> tex.checker.color1.g
                >> tex.checker.color1.b;
            in >> tex.checker.color2.r >> tex.checker.color2.g
                >> tex.checker.color2.b;
            in >> tex.checker.scale;

            tex.checker.map[0][0] = tex.checker.color1;
            tex.checker.map[0][1] = tex.checker.color2;
            tex.checker.map[1][0] = tex.checker.color1;
            tex.checker.map[1][1] = tex.checker.color2;
        }
        else if(type == "texmap")
		{
            tex.type = MapTexture;
            in >> tex.map.filename;
            in >> tex.map.p0.x >> tex.map.p0.y >> tex.map.p0.z >> tex.map.p0.w;
            in >> tex.map.p1.x >> tex.map.p1.y >> tex.map.p1.z >> tex.map.p1.w;
            tex.map.image.load(tex.map.filename);
        }
        else
		{
            printf("invalid texture type (%s).\n", type.c_str());
        }
        textures.push_back(tex);
    }
}

void Scene::parse_material(std::istream &in)
{
    Material material;

    in >> num_materials;
    if(in.fail())
    {
        printf("Failed to read the input file in material description.\n");
        exit(1);
    }
    for(int i = 0; i < num_materials; ++i)
	{
        material.id = i;
        in >> material.kA >> material.kD >> material.kS >> material.alpha
           >> material.kR >> material.kT >> material.ior;

        double norm = material.kA + material.kD + material.kS + material.kR + material.kT;
        material.kA /= norm;
        material.kD /= norm;
        material.kS /= norm;
        material.kR /= norm;
        material.kT /= norm;

        materials.push_back(material);
    }
}

void Scene::parse_object(std::istream &in)
{
    std::string type;

    in >> num_objects;
    if(in.fail())
    {
        printf("Failed to read the input file in object description.\n");
        exit(1);
    }
    for(int i = 0; i < num_objects; ++i)
	{
        int texId, matId;
		in >> texId >> matId;
        in >> type;

        if(type == "sphere")
		{
			SphereObject* object = new SphereObject();
			object->id = i;
			object->texture = textures[texId];
			object->material = materials[matId];

			in >> object->original_pos.x >> object->original_pos.y >> object->original_pos.z;
            in >> object->radius;
			in >> object->original_scale.x >> object->original_scale.y >> object->original_scale.z;
			in >> object->acceleration.x >> object->acceleration.y >> object->acceleration.z;
			in >> object->emission.r >> object->emission.g >> object->emission.b;
            object->calculate_matrices();
			if(object->emission.r > 0 || object->emission.g > 0 || object->emission.b > 0)
                luminaries.push_back(object);

            if(object->material.kA > 0){
                in >> object->material.absortion.r;
                in >> object->material.absortion.g;
                in >> object->material.absortion.b;
            }


			objects.push_back(object);
        }
        else if(type == "polyhedron")
		{
			PolyhedronObject* object = new PolyhedronObject();
			object->id = i;
			object->texture = textures[texId];
			object->material = materials[matId];

            in >> object->numFaces;
            if(in.fail())
                printf("Failed to read the input file in polyhedron faces description.\n");

            Plane plane;
            for(int j = 0; j < object->numFaces; ++j)
			{
                in >> plane.a >> plane.b >> plane.c >> plane.d;
                object->planes.push_back(plane);
            }
			objects.push_back(object);
        }
		else if(type == "torus")
		{
			TorusObject* object = new TorusObject();
			object->id = i;
			object->texture = textures[texId];
			object->material = materials[matId];

			in >> object->radius >> object->thickness;
			in >> object->original_pos.x >> object->original_pos.y >> object->original_pos.z;
			in >> object->original_rot.x >> object->original_rot.y >> object->original_rot.z;
			in >> object->original_scale.x >> object->original_scale.y >> object->original_scale.z;
			in >> object->acceleration.x >> object->acceleration.y >> object->acceleration.z;
			object->calculate_matrices();
			objects.push_back(object);
		}
		else if(type == "cylinder")
		{
			CylinderObject* object = new CylinderObject();
			object->id = i;
			object->texture = textures[texId];
			object->material = materials[matId];

			in >> object->bottom >> object->top >> object->radius;
			in >> object->original_pos.x >> object->original_pos.y >> object->original_pos.z;
			in >> object->original_rot.x >> object->original_rot.y >> object->original_rot.z;
			in >> object->original_scale.x >> object->original_scale.y >> object->original_scale.z;
			in >> object->acceleration.x >> object->acceleration.y >> object->acceleration.z;
			object->calculate_matrices();
			objects.push_back(object);
		}
        else
		{
            printf("invalid object type(%s).\n", type.c_str());
        }
    }
}


