#include "object.hpp"
#include <fstream>
#include <sstream>
#include <map>
#include <vector>
#include <cmath>
#include <cfloat>
#include <iostream>

// 2016 Version

bool Object::intersect(Ray ray, Intersection &hit) const 
{
    // Assert the correct values of the W coords of the origin and direction.
    // You can comment this out if you take great care to construct the rays
    // properly.
    ray.origin[3] = 1;
    ray.direction[3] = 0;

    Ray localRay(i_transform * ray.origin, i_transform * ray.direction);
	//!!! USEFUL NOTES: to calculate depth in localIntersect(), if the intersection happens at
	//ray.origin + ray.direction * t, then t is just the depth
	//!!! USEFUL NOTES: Here direction might be scaled, so you must not renormalize it in
	//localIntersect(), or you will get a depth in local coordinate system,
	//which can't be compared with intersections with other objects
    if (localIntersect(localRay, hit))
	{
        // Assert correct values of W.
        hit.position[3] = 1;
        hit.normal[3] = 0;
        
		// Transform intersection coordinates into global coordinates.
        hit.position = transform * hit.position;
        hit.normal = (n_transform * hit.normal).normalized();
        
		return true;
    }

    return false;
}

//INTERSECTIONS 

bool Sphere::localIntersect(Ray const &ray, Intersection &hit) const 
{
	//////////////////
    // YOUR CODE HERE 
	// For this part you are required to calculate if a ray has intersected your sphere.
	// Be sure to cover all the possible intersection scenarios(zero, one, and two points of intersection).
	// Test your result by comparing the output depth image of your algorithm with the provided example solution's results.

	// Here in local coordinate system, the sphere is centered on (0, 0, 0)
	//with radius 1.0
	//
	// NOTE: hit.depth is the current closest intersection depth, so don't
	// accept any intersection that happens further away than that.

	// The formula for a sphere: 
	// (x - Cx)^2 + (y - Cy)^2 + (Z - Cz)^2 - R^2 = 0
	// But since the sphere is centered at the origin, that simplifies to
	// x^2 + y^2 + z^2 - r^2 = 0

	// The parametric equations for a ray in the direction <d_x, d_y, d_z>
	// X = x0 + d_x*t
	// Y = y0 + d_y*t
	// Z = z0 + d_z*t

	double x_0 = ray.origin[0];
	double y_0 = ray.origin[1];
	double z_0 = ray.origin[2];



	double d_x = ray.direction[0];
	double d_y = ray.direction[1];
	double d_z = ray.direction[2];

	double t;
	double r = this->radius;

	double a = ray.direction.length2();
	double b = 2 * ((x_0 * d_x) + (y_0 * d_y) + (z_0 * d_z));
	double c = ray.origin.length2() - pow(r, 2);

	double d = (pow(b, 2)) - (4 * a * c);
	bool inSphere = false;

	if (d < 0) { // no solution
		return false;
	}
	else if (d > 0) { //two solutions
		double t1 = (-b + sqrt(d)) / (2 * a);
		double t2 = (-b - sqrt(d)) / (2 * a);
		if (t1 >= 0 && t2 >= 0) {

			if (t1 < t2) {
				t = t1;
			}
			else {
				t = t2;
			}
		}
		else if (t1 >= 0 && t2 < 0) {
			t = t1;
		}
		else if (t2 >= 0 && t1 < 0) {
			t = t2;
		}
		else {
			return false;
		}
	}
	else if (d == 0) { // one solution
		double t0 = (-b) / (2 * a);
		if (t0 >= 0) {
			t = t0;
		}
		else {
			return false;
		}
	}
	if (hit.depth > t) {
		hit.depth = t;
		Vector ray_position = Vector((x_0 + d_x * t), (y_0 + d_y * t), (z_0 + d_z * t)); // our parametric ray
		hit.position = ray_position;
		hit.normal = ray_position.normalized();
		return true;
	}
	else {
		return false;
	}
}


bool Plane::localIntersect(Ray const &ray, Intersection &hit) const
{
	//////////////////
	// YOUR CODE HERE 
	// The implementation of this part is similar to the previous part in that you are 
	// calculating if a line has intersected your plane.Test this function in a similar 
	// way to the previous part(note that as you do new objects will appear).

	// Do not accept intersection while ray is inside the plane
	// Here in local coordinate system, the plane is at z = 0
	//
	// NOTE: hit.depth is the current closest intersection depth, so don't
	// accept any intersection that happens further away than that.

	
	double x_0 = ray.origin[0];
	double y_0 = ray.origin[1];
	double z_0 = ray.origin[2];

	double d_x = ray.direction[0];
	double d_y = ray.direction[1];
	double d_z = ray.direction[2];

	// since the plane equation is z = 0, we can derive t ourselves
	// z_0 + t*d_z = 0
	// t  = - z_0/d_z

	double t = -z_0 / d_z;

	if (t < 0) {
		return false;
	}
	if (hit.depth > t) {
		hit.depth = t;
		Vector ray_position = Vector((x_0 + d_x * t), (y_0 + d_y * t), (z_0 + d_z * t)); // our parametric ray
		hit.position = ray_position;
		hit.normal = Vector(0, 0, 1, 0); // since z = 0, the normal points up from z
		return true;
	} 
    return false;
}

bool Mesh::intersectTriangle(Ray const &ray, Triangle const &tri, Intersection &hit) const
{
	// Extract vertex positions from the mesh data.
	Vector const &p0 = positions[tri[0].pi];
	Vector const &p1 = positions[tri[1].pi];
	Vector const &p2 = positions[tri[2].pi];

	//////////////////
	// YOUR CODE HERE 
	// Think back through the course and try to decide what equations might help you decide on which 
	// side of the bounding lines of the triangle the ray intersects.
	// Test this part just like the above two parts; when triangle intersection is working properly, 
	// you should be able to see full meshes appear in your scenes.
	//
	// NOTE: hit.depth is the current closest intersection depth, so don't
	// accept any intersection that happens further away than that.
	//!!! USEFUL NOTES: for the intersection point, its normal should satisfy hit.normal.dot(ray.direction) < 0


	//algo c/o scratchapixel
	double x_0 = ray.origin[0];
	double y_0 = ray.origin[1];
	double z_0 = ray.origin[2];

	double d_x = ray.direction[0];
	double d_y = ray.direction[1];
	double d_z = ray.direction[2];
	
	Vector p0p1 = p1 - p0;
	Vector p0p2 = p2 - p0;
	Vector normal = p0p1.cross(p0p2);
	normal.normalize();

	float normal_dot_dir = normal.dot(ray.direction);
	if (fabs(normal_dot_dir) < 0.0001) {
		return false;
	}

	float d = normal.dot(p0);
	float t = (-normal.dot(ray.origin) + d) / normal.dot(ray.direction);

	if (t < 0) {
		return false;
	}

	Vector ray_position = ray.origin + t * ray.direction;

	Vector edge0 = p1 - p0;
	Vector edge1 = p2 - p1;
	Vector edge2 = p0 - p2;
	Vector ray_p0 = ray_position - p0;
	Vector ray_p1 = ray_position - p1;
	Vector ray_p2 = ray_position - p2;

	if ((normal.dot(edge0.cross(ray_p0)) < 0) || (normal.dot(edge1.cross(ray_p1)) < 0) || (normal.dot(edge2.cross(ray_p2)) < 0)) {
		return false;
	}

	if (hit.depth > t) {
		hit.depth = t;
		hit.position = ray_position;
		hit.normal = normal.normalized();
		return true;
	}
	return false;
	
}

bool Conic::localIntersect(Ray const &ray, Intersection &hit) const {
	//////////////////
	// YOUR CODE HERE (creative license)
    return false;
}


// Intersections!
bool Mesh::localIntersect(Ray const &ray, Intersection &hit) const
{
	// Bounding box check
	double tNear = -DBL_MAX, tFar = DBL_MAX;
	for (int i = 0; i < 3; i++) {
		if (ray.direction[i] == 0.0) {
			if (ray.origin[i] < bboxMin[i] || ray.origin[i] > bboxMax[i]) {
				// Ray parallel to bounding box plane and outside of box!
				return false;
			}
			// Ray parallel to bounding box plane and inside box: continue;
		}
		else {
			double t1 = (bboxMin[i] - ray.origin[i]) / ray.direction[i];
			double t2 = (bboxMax[i] - ray.origin[i]) / ray.direction[i];
			if (t1 > t2) std::swap(t1, t2); // Ensure t1 <= t2

			if (t1 > tNear) tNear = t1; // We want the furthest tNear
			if (t2 < tFar) tFar = t2; // We want the closest tFar

			if (tNear > tFar) return false; // Ray misses the bounding box.
			if (tFar < 0) return false; // Bounding box is behind the ray.
		}
	}
	// If we made it this far, the ray does intersect the bounding box.

	// The ray hits the bounding box, so check each triangle.
	bool isHit = false;
	for (size_t tri_i = 0; tri_i < triangles.size(); tri_i++) {
		Triangle const &tri = triangles[tri_i];

		if (intersectTriangle(ray, tri, hit)) {
			isHit = true;
		}
	}
	return isHit;
}

double Mesh::implicitLineEquation(double p_x, double p_y,
	double e1_x, double e1_y,
	double e2_x, double e2_y) const
{
	return (e2_y - e1_y)*(p_x - e1_x) - (e2_x - e1_x)*(p_y - e1_y);
}


