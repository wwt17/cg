// C++ include
#include <iostream>
#include <string>
#include <vector>
#include <cmath>

// Utilities for the Assignment
#include "utils.h"

// Image writing library
#define STB_IMAGE_WRITE_IMPLEMENTATION // Do not include this line twice in your project!
#include "stb_image_write.h"


// Shortcut to avoid Eigen:: everywhere, DO NOT USE IN .h
using namespace Eigen;


const double eps = 1e-7;
const size_t H = 800, W = 800;

template<class T>
T sqr(const T x) {
	return x * x;
}


enum ProjectionType {orthographic, perspective};
const char *ProjectionTypeString[] = {"orthographic", "perspective"};


class Object {
public:
	virtual bool intersect_with_ray(const Vector3d &ray_origin, const Vector3d &ray_direction, Vector3d &ray_intersection, Vector3d &ray_normal) const = 0;
};


class Sphere: public Object {
public:
	Vector3d center;
	double radius;
	Sphere(const Vector3d &center, const double radius): center(center), radius(radius) {
	}
	bool intersect_with_ray(const Vector3d &ray_origin, const Vector3d &ray_direction, Vector3d &ray_intersection, Vector3d &ray_normal) const {
		Vector3d v = ray_origin - center;
		double a = ray_direction.squaredNorm(), b_over_2 = ray_direction.dot(v), c = v.squaredNorm() - sqr(radius);
		double delta = sqr(b_over_2) - a * c, t;
		if (delta >= 0) {
			t = (-b_over_2 - std::sqrt(delta)) / a;
			ray_intersection = ray_origin + t * ray_direction;
			ray_normal = (ray_intersection - center).normalized();
			return t >= 0;
		}
		return false;
	}
};


class Parallelogram: public Object {
public:
	Vector3d origin, u, v;
	Parallelogram(const Vector3d &origin, const Vector3d &u, const Vector3d &v): origin(origin), u(u), v(v) {
	}
	bool intersect_with_ray(const Vector3d &ray_origin, const Vector3d &ray_direction, Vector3d &ray_intersection, Vector3d &ray_normal) const {
		/* The basic idea is: set up a new basis (u, v, -ray_direction)
		 * with origin at origin and compute the coordinate of ray_origin in
		 * this basis.
		*/

		// The basis degenerates (i.e., volume of its span is 0)
		if (std::abs(u.cross(v).dot(-ray_direction)) <= eps)
			return false;

		// Solve the coordinate (x', y', z') so that u * x' + v * y' + (-ray_direction) * z' = ray_origin - origin
		Matrix3d basis; basis << u, v, -ray_direction;
		Vector3d coord = basis.inverse() * (ray_origin - origin);
		ray_intersection = origin + u * coord(0) + v * coord(1);

		// Compute normal
		ray_normal = u.cross(v).normalized();
		if ((ray_origin - ray_intersection).dot(ray_normal) < 0)
			ray_normal = -ray_normal;

		// Intersect if within the span of u and v and ray_direction is onto the plane
		return coord(0) >= 0 && coord(0) <= 1 && coord(1) >= 0 && coord(1) <= 1 && coord(2) >= 0;
	}
};


void raytrace(
	const Object *object,
	// RBGA Color (in [0, 1])
	const Array4d color = Array4d(1, 1, 1, 1),
	// Single light source
	const Vector3d light_position = Vector3d(-1, 1, 1),
	const double diffuse_weight = 1,
	const double specular_weight = 0,
	const double specular_power = 1000,
	const double ambient = 0,
	// Projection type
	const ProjectionType projection_type = perspective,
	// Perspective position
	const Vector3d perspective_position = Vector3d(0, 0, 2),
	// Ray direction for orthographic projection
	const Vector3d orthographic_direction = Vector3d(0, 0, -1),
	// Image size (H, W)
	const size_t H = 800,
	const size_t W = 800,
	// Output file
	const std::string filename = "raytrace.png"
) {
	MatrixXd
		// brightness of color
		C = MatrixXd::Zero(W, H),
		// alpha mask
		A = MatrixXd::Zero(W, H),
		// diffuse and specular values
		diffuse = MatrixXd::Zero(W, H),
		specular = MatrixXd::Zero(W, H);

	// The camera is perspective/orthographic, pointing in the direction -z and covering the unit square (-1,1) in x and y
	Vector3d origin(-1, 1, 1);
	Vector3d x_displacement(2.0 / W, 0, 0);
	Vector3d y_displacement(0, -2.0 / H, 0);

	for (unsigned i=0; i < W; ++i) {
		for (unsigned j=0; j < H; ++j) {
			// Prepare the ray
			Vector3d ray_origin = origin + double(i)*x_displacement + double(j)*y_displacement;
			Vector3d ray_direction = projection_type == perspective
				? (ray_origin - perspective_position).normalized()
				: orthographic_direction;

			// Check if the ray intersects with the parallelogram
			Vector3d ray_intersection, ray_normal;
			if (object->intersect_with_ray(ray_origin, ray_direction, ray_intersection, ray_normal)) {
				// The ray hit the parallelogram, the exact intersection point is ray_intersection and the normal at the intersection point is ray_normal

				Vector3d l = (light_position - ray_intersection).normalized();

				// diffuse shading
				diffuse(i,j) = std::max(0., l.dot(ray_normal));

				// specular shading
				if (specular_weight) {
					Vector3d v = (perspective_position - ray_intersection).normalized();
					Vector3d h = (v + l).normalized();
					specular(i,j) = pow(std::max(0., h.dot(ray_normal)), specular_power);
				}

				// Simple ambient + diffuse + specular model
				C(i,j) = ambient + diffuse_weight * diffuse(i,j) + specular_weight * specular(i,j);

				// Clamp to zero
				C(i,j) = std::max(C(i,j),0.);

				// Disable the alpha mask for this pixel
				A(i,j) = 1;
			}
		}
	}

	// Save to png
	write_matrix_to_png(C * color(0), C * color(1), C * color(2), A * color(3), filename);
}


int main() {
	Sphere sphere(Vector3d(0, 0, 0), 0.9);
	Parallelogram pgram(
		Vector3d(-0.5, 0.6, -0.1),
		Vector3d(0.6, -0.2, -0.2),
		Vector3d(0.4, -0.5, -0.1));

	// RGBA Color
	Array4d color(87/255., 6/255., 140/255., 1);
	// Light source
	Vector3d light_position(-1, 1, 1);
	// Perspective position
	Vector3d perspective_position(0, 0, 3);
	// Ray direction for orthographic projection
	Vector3d orthographic_direction(0, 0, -1);

	Object *objects[2] = {&sphere, &pgram};
	std::string object_names[2] = {"sphere", "parallelogram"};

	for (ProjectionType projection: {orthographic, perspective}) {
		const char *projection_string = ProjectionTypeString[projection];
		for (int i = 0; i < 2; i++) {
			std::cout << "Simple ray tracer, one " << object_names[i] << " with " << projection_string << " projection" << std::endl;
			raytrace(
				objects[i],
				color,
				light_position,
				1, 0, 0, 0,
				projection,
				perspective_position,
				orthographic_direction,
				H, W,
				object_names[i] + "_" + projection_string + ".png"
			);

			double specular_power = 1000, ambient = 0.1;
			std::cout << "Simple ray tracer, one " << object_names[i] << " with " << projection_string << " and different shading" << std::endl;
			raytrace(
				objects[i],
				color,
				light_position,
				1, 1, specular_power, ambient,
				projection,
				perspective_position,
				orthographic_direction,
				H, W,
				object_names[i] + "_" + projection_string + "_shading.png"
			);
		}
	}

	return 0;
}
