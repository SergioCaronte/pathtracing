A dependency free c++ pathtracing project coded as an assignment of a computer graphics course. 

Coded and tested using Code::Blocks and G++ on Ubuntu 14.04.

### How do I get set up? ###

On Windows:

Open Code::Blocks project (Pathtracing.cbp) and compile.

On Linux:

Open terminal, travel to root folder and run compile.sh (./compile.sh)

### Running ###

It works in command line. Call the program name passing five parameters
input file that describes the scene.
output file it is the pathtraced image (it only saves as ppm file).
number of samples per pixel
[image width]
[image height]

for example:
On Windows:
Pathtracing.exe input.in output.ppm 10 800 600 

On Linux:
./Pathtracing input.in output.ppm 10 800 600

### Features ###

Multi-sampling using Multi-jittering.
Collision with spheres, planes, torus and cylinders.
Support to three kinds of texture (solid, checker and image(.ppm)).
Simulate lens aperture and depth of field.
Support to area lights (plane, spheres, hemispheres and circles).
Rigid transformations
BRDFs (Lambertian, Phong, Torrance-Sparrow, Reflection, Refraction)

### Input file ###

In order to configure a 3D scene, it must provide an text file as follows:

The first four lines describe the camera parameters:

- CameraCenterPosition (3 floats)

- CameraLookAtPosition (3 floats)

- CameraVectorUp (3 floats)

- FOVy (1 float) -NumSamplesPerPixel (1 int) -RadiusLens (1 float) -FocalLength (1 float) -ShutterSpeed (1 int) -ExposureValue (1 int)


A integer value informing number of light sources. 
!This project does not uses this information for light source. It is garbage from an old raytracing project.
!The emission of sphere object in object description is used as light source instead.

- NumberOfLightSources (1 int)

First light is used as ambient

- Position (3 floats) -SpectrumRGB (3 floats) -OpenGLFallOff (3 floats)

Following lights

- Position (3 floats) -SpectrumRGB (3 floats) -OpenGLFallOff (3 floats) -TypeOfLightArea (string) -NumberOfSamplingPerLight (1 int) -LightAreaSize (1 float)

A integer value informing number of textures

- NumberOfTextures (1 int)

Following lines describe textures

- Type(string) following parameters vary depending of the texture type

- solid -SpectrumRGB (3 floats)

- checker -SpectrumRGB1 (3 floats) -SpectrumRGB2 (3 floats) -SquareSize (1 int)

- texmap -FileName(ppm image name) -VectorOfInterpolation1 (4 floats)  -VectorOfInterpolation2 (4 floats)  

A integer value informing number of materials

- NumberOfTextures (1 int)

Following lines describe material coefficients (The total sums of all coefficients must be lesser or equal to 1.0)

- kTorranceSparrow (1 float) -kLambertian (1 float) -kBlinnPhong( 1 float) -kSpecularAlpha (1 int) -kReflection (1 float) -kTransmission (1 float) -kIndexOfRefraction(a.k.a n2) (1 float)  

A integer value informing number of objects

- NumberOfObjects (1 int)

Following lines describe objects

- index of texture (1 int) -index of material (1 int) -typeOfObject (string) -following parameters vary depending of the object type

- sphere -CenterPosition (3 floats) -Radius (1 float) -Scale (3 floats) -Acceleration (3 floats) -Emission ( 3 floats)

- cylinder -BottomZPosition (1 float) -TopZPosition (1 float) -Radius (1 float) -CenterPosition (3 floats) -Rotation (3 floats) - Acceleration (3 floats)

- torus -OutsideRadius (1 float) -InsideRadius (1 float)  -CenterPosition (3 floats) -Rotation (3 floats) - Acceleration (3 floats)

- polyhedron -NumberOfPlanes (1 int) 

For each polyhedron plane -NormalVector (3 floats) -Distance (1 float) 
	

