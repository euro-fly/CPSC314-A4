#include <vector>
#include <cstdlib>
#include <cfloat>
#include <cmath>
#include <map>
#include <vector>
#include <cstdio>
#include <string>
#include <fstream>
#include <algorithm>
#include <iostream>
#include <sstream>
#include "raytracer.hpp"
#include "image.hpp"

//2016 Version

void Raytracer::render(const char *filename, const char *depth_filename, Scene const &scene)
{
    // Allocate the two images that will ultimately be saved.
    Image colorImage(scene.resolution[0], scene.resolution[1]);
    Image depthImage(scene.resolution[0], scene.resolution[1]);
    
    // Create the zBuffer.
    double *zBuffer = new double[scene.resolution[0] * scene.resolution[1]];
    for(int i = 0; i < scene.resolution[0] * scene.resolution[1]; i++) {
        zBuffer[i] = DBL_MAX;
    }

	//////////////////
	// YOUR CODE HERE 
	// calculate camera parameters for rays, refer to the slides for details
	//!!! USEFUL NOTES: tan() takes rad rather than degree, use deg2rad() to transform
	//!!! USEFUL NOTES: view plane can be anywhere, but it will be implemented differently,
	//you can find references from the course slides 22_GlobalIllum.pdf

	Vector camera_position = scene.camera.position;
	Vector camera_center = scene.camera.center;

	Vector w_vector = camera_center - camera_position;
	w_vector.normalize();

	Vector v_vector = scene.camera.up;
	v_vector.normalize();

	Vector u_vector = w_vector.cross(v_vector);
	u_vector.normalize();
	

	double angle = tan(deg2rad(scene.camera.fov / 2));
	//double atangent = atan(deg2rad(scene.camera.fovy)/2);
	// get length of top from centre of image plane
	double top = scene.camera.zNear * angle;
	double right = top * scene.camera.aspect;
	double left = -right;
	double bottom = -top;

	// calculate vector from camera to left top of image plane
	Vector centerVec = camera_position + (scene.camera.zNear * w_vector);
	Vector oVec = centerVec + (left * u_vector) + (bottom * v_vector);
	double deltaU = (right - left) / scene.resolution[0];
	double deltaV = (top - bottom) / scene.resolution[1];

    // Iterate over all the pixels in the image.
    for(int y = 0; y < scene.resolution[1]; y++) {
        for(int x = 0; x < scene.resolution[0]; x++) {

            // Generate the appropriate ray for this pixel
			Ray ray;
			if (scene.objects.empty())
			{
				//no objects in the scene, then we render the default scene:
				//in the default scene, we assume the view plane is at z = 640 with width and height both 640
				ray = Ray(scene.camera.position, (Vector(-320, -320, 640) + Vector(x + 0.5, y + 0.5, 0) - scene.camera.position).normalized());
			}
			else
			{
				//////////////////
				// YOUR CODE HERE
				// set primary ray using the camera parameters
				//!!! USEFUL NOTES: all world coordinate rays need to have a normalized direction

				// calculate the vector at the pixel (x,y) using equation in slides
				Vector changeU = (x + 0.5) * deltaU * u_vector;
				Vector changeY = (y + 0.5) * deltaV * v_vector;
				Vector ray_position = oVec + changeU + changeY;

				Vector ray_direction = ray_position - camera_position;
				ray_direction.normalize();

				ray = Ray(camera_position, ray_direction);
				
			}

            // Initialize recursive ray depth.
            int ray_depth = 0;
           
            // Our recursive raytrace will compute the color and the z-depth
            Vector color;

            // This should be the maximum depth, corresponding to the far plane.
            // NOTE: This assumes the ray direction is unit-length and the
            // ray origin is at the camera position.
            double depth = scene.camera.zFar;

            // Calculate the pixel value by shooting the ray into the scene
            trace(ray, ray_depth, scene, color, depth);

            // Depth test
            if(depth >= scene.camera.zNear && depth <= scene.camera.zFar && 
                depth < zBuffer[x + y*scene.resolution[0]]) {
                zBuffer[x + y*scene.resolution[0]] = depth;

                // Set the image color (and depth)
                colorImage.setPixel(x, y, color);
                depthImage.setPixel(x, y, (depth-scene.camera.zNear) / 
                                        (scene.camera.zFar-scene.camera.zNear));
            }
        }

		//output step information
		if (y % 100 == 0)
		{
			printf("Row %d pixels done.\n", y);
		}
    }

	//save image
    colorImage.writeBMP(filename);
    depthImage.writeBMP(depth_filename);

	printf("Ray tracing terminated and images are saved.\n");

    delete[] zBuffer;
}


bool Raytracer::trace(Ray const &ray, int &ray_depth, Scene const &scene, Vector &rayOutColor, double &depth)
{
    // Increment the ray depth.
	ray_depth++;

    // - iterate over all objects calling Object::intersect.
    // - don't accept intersections not closer than given depth.
    // - call Raytracer::shade with the closest intersection.
    // - return true iff the ray hits an object.
	if (scene.objects.empty())
	{
		// no objects in the scene, then we render the default scene:
		// For default, we assume there's a cube centered on (0, 0, 1280 + 160) with side length 320 facing right towards the camera
		// test intersection:
		double x = 1280 / ray.direction[2] * ray.direction[0] + ray.origin[0];
		double y = 1280 / ray.direction[2] * ray.direction[1] + ray.origin[1];
		if ((x <= 160) && (x >= -160) && (y <= 160) && (y >= -160))
		{
			//if intersected:
			Material m; m.emission = Vector(16.0, 0, 0); m.reflect = 0; //just for default material, you should use the intersected object's material
			Intersection intersection;	//just for default, you should pass the intersection found by calling Object::intersect()
			rayOutColor = shade(ray, ray_depth, intersection, m, scene);
			depth = 1280;	//the depth should be set inside each Object::intersect()
		}
	}
	else
	{
		//////////////////
		// YOUR CODE HERE
		// Note that for Object::intersect(), the parameter hit is the current hit
		// your intersect() should be implemented to exclude intersection far away than hit.depth

		Intersection hit;
		hit.depth = depth;

		for (int j = 0; j < scene.objects.size(); ++j) {
			if (scene.objects[j]->intersect(ray, hit)) {
				if (hit.depth <= depth) {
					depth = hit.depth;
					rayOutColor = shade(ray, ray_depth, hit, scene.objects[j]->material, scene);
				}
			}
		}
	}

    // Decrement the ray depth.
	ray_depth--;

    return false; 
}


Vector Raytracer::shade(Ray const &ray, int &ray_depth, Intersection const &intersection, Material const &material, Scene const &scene)
{

    // - iterate over all lights, calculating ambient/diffuse/specular contribution
    // - use shadow rays to determine shadows
    // - integrate the contributions of each light
    // - include emission of the surface material
    // - call Raytracer::trace for reflection/refraction colors
    // Don't reflect/refract if maximum ray recursion depth has been reached!
	//!!! USEFUL NOTES: attenuate factor = 1.0 / (a0 + a1 * d + a2 * d * d)..., ambient light doesn't attenuate, nor does it affected by shadow
	//!!! USEFUL NOTES: don't accept shadow intersection far away than the light position
	//!!! USEFUL NOTES: for each kind of ray, i.e. shadow ray, reflected ray, and primary ray, the accepted furthest depth are different
// !!!!! edited lines start 
	Vector diffuse(0);
	Vector ambient(0);
	Vector specular(0);		
// !!!!! edited lines end
	double shadow_bias = 0.00001;
	for (auto lightIter = scene.lights.begin(); lightIter != scene.lights.end(); lightIter++)
	{
		//////////////////
		// YOUR CODE HERE 
		// First you can assume all the light sources are directly visible. You should calculate the ambient, diffuse, 
		// and specular terms.You should think of this part in terms of determining the color at the point where the ray 
		// intersects the scene.
		// After you finished, you will be able to get the colored resulting image with local illumination, just like in programming assignment 3.

		Vector Ia = material.ambient * lightIter->ambient;

		Vector n = intersection.normal.normalized();
		Vector l = lightIter->position - intersection.position;
		l.normalize();

		float n_dot_l = n.dot(l) > 0.0 ? n.dot(l) : 0.0; // if n dot l is < 0, make it 0 - total internal reflection

		// Id = kd*Il*(n dot l)
		Vector Id = material.diffuse * lightIter->diffuse * n_dot_l;

		// we want to reflect l about n...
		Vector r = l - (2 * (l.dot(n)) * n);
		r.normalize();

		Vector v = ray.direction.normalized();

		float v_dot_r = v.dot(r) > 0.0 ? v.dot(r) : 0.0;

		Vector Is = material.specular * lightIter->specular * pow(v_dot_r, material.shininess);


		//////////////////
		// YOUR CODE HERE 
		// Emit the shadow ray from a point you're computing direct illumination for to determine which lights 
		// are contributing to the lighting at that point.Be careful to exclude the origin of the ray from the 
		// intersection points, but do remember that the intersection points could be other points on the same 
		// object if the object is not convex(for example, a teapot).
		// For points in the shadow, scale their original lighting color by the factor  (1 - material.shadow)

		Intersection shadow;
		
		Ray shadow_ray = Ray(intersection.position + (l * shadow_bias), l);
		double light_distance = (lightIter->position - intersection.position).length();
		shadow.depth = light_distance;

		bool in_shadow = false;

		for (int j = 0; j < scene.objects.size(); ++j) {
			if (scene.objects[j]->intersect(shadow_ray, shadow)) {
				in_shadow = true;
				break;
			}
		}
		// our attenuation = CONSTANT + LINEAR * DISTANCE + QUADRATIC * DISTANCE^2
		double attenuation = lightIter->attenuation[0] + (lightIter->attenuation[1] * light_distance) + (lightIter->attenuation[2] * pow(light_distance, 2));
		// and luminosity = 1/attenuation
		double luminosity = 1 / attenuation;

		if (in_shadow) {
			Is = Is * (1 - material.shadow);
			Id = Id * (1 - material.shadow);
		}

		Is = Is * luminosity;
		Id = Id * luminosity;

		// update our ambient, specular, and diffuse lights by adding these ones
		ambient += Ia;
		specular += Is;
		diffuse += Id;

		//////////////////
		// YOUR CODE HERE 
		// Use the ray_depth recursion depth variable to stop the recursion process. (The default used in the solution is 10.) 
		// Update the lighting computation at each step to account for the secondary component.
		// You can think of this part as an extended shadow ray calculation, recursively iterating to determine contributing 
		// light(and weighting newly determined light sources into the original pixel).

		// okay tbh i have no idea what goes here

	}

	Vector reflectedLight(0);
	if ((!(ABS_FLOAT(material.reflect) < 1e-6)) && (ray_depth < MAX_RAY_RECURSION))
	{
		//////////////////
		// YOUR CODE HERE 
		// calculate reflected color using trace() recursively

		Vector n = intersection.normal.normalized();

		// get the reflection of the ray about n
		Vector r = ray.direction - (2 * (ray.direction.dot(n)) * n);
		r.normalize();

		Ray reflection = Ray((intersection.position + (r * shadow_bias)), r);
		double depth = DBL_MAX;

		trace(reflection, ray_depth, scene, reflectedLight, depth);
	}


	// !!!! edited line starts
	return material.emission + ambient + diffuse + specular + material.reflect * reflectedLight;  
	// !!!! edited line ends
}