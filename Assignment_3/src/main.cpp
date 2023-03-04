////////////////////////////////////////////////////////////////////////////////
// C++ include
#include <iostream>
#include <string>
#include <vector>
#include <limits>
#include <cmath>

// Utilities for the Assignment
#include "utils.h"

// Image writing library
#define STB_IMAGE_WRITE_IMPLEMENTATION // Do not include this line twice in your project!
#include "stb_image_write.h"

// Shortcut to avoid Eigen:: everywhere, DO NOT USE IN .h
using namespace Eigen;

////////////////////////////////////////////////////////////////////////////////
// Global variables, functions and classes
////////////////////////////////////////////////////////////////////////////////
typedef Vector4d Color; // RBGA Color (in [0, 1])

const double eps = 1e-7, inf = 1./0., pi = acos(-1);

template<class T>
T sqr(const T x) {
	return x * x;
}


enum ProjectionType {orthographic, perspective};
const char *ProjectionTypeString[] = {"orthographic", "perspective"};


class Object {
public:
	virtual double intersect_with_ray(const Vector3d &ray_origin, const Vector3d &ray_direction, Vector3d &ray_intersection, Vector3d &ray_normal) const = 0;
};


class Sphere: public Object {
public:
	Vector3d center;
	double radius;
	Sphere(const Vector3d &center, const double radius): center(center), radius(radius) {
	}
	double intersect_with_ray(const Vector3d &ray_origin, const Vector3d &ray_direction, Vector3d &ray_intersection, Vector3d &ray_normal) const {
		Vector3d v = ray_origin - center;
		double a = ray_direction.squaredNorm(), b_over_2 = ray_direction.dot(v), c = v.squaredNorm() - sqr(radius);
		double delta = sqr(b_over_2) - a * c, t;
		if (delta >= 0) {
			t = (-b_over_2 - std::sqrt(delta)) / a;
			ray_intersection = ray_origin + t * ray_direction;
			ray_normal = (ray_intersection - center).normalized();
			return t;
		}
		return -inf;
	}
};


class Parallelogram: public Object {
public:
	Vector3d origin, u, v;
	Parallelogram(const Vector3d &origin, const Vector3d &u, const Vector3d &v): origin(origin), u(u), v(v) {
	}
	double intersect_with_ray(const Vector3d &ray_origin, const Vector3d &ray_direction, Vector3d &ray_intersection, Vector3d &ray_normal) const {
		/* The basic idea is: set up a new basis (u, v, -ray_direction)
		 * with origin at origin and compute the coordinate of ray_origin in
		 * this basis.
		*/

		// The basis degenerates (i.e., volume of its span is 0)
		if (std::abs(u.cross(v).dot(-ray_direction)) <= eps)
			return -inf;

		// Solve the coordinate (x', y', z') so that u * x' + v * y' + (-ray_direction) * z' = ray_origin - origin
		Matrix3d basis; basis << u, v, -ray_direction;
		Vector3d coord = basis.inverse() * (ray_origin - origin);
		ray_intersection = origin + u * coord(0) + v * coord(1);

		// Compute normal
		ray_normal = u.cross(v).normalized();
		if ((ray_origin - ray_intersection).dot(ray_normal) < 0)
			ray_normal = -ray_normal;

		// Intersect if within the span of u and v and ray_direction is onto the plane
		return coord(0) >= 0 && coord(0) <= 1 && coord(1) >= 0 && coord(1) <= 1 ? coord(2) : -inf;
	}
};


class Light {
public:
	Vector3d position;
	Color color;
	Light(const Vector3d &position, const Color &color): position(position), color(color) {
	}
};


////////////////////////////////////////////////////////////////////////////////
// Scene setup, global variables
////////////////////////////////////////////////////////////////////////////////
const std::string filename("raytrace.png");

//Camera settings
const double focal_length = 10;
const double field_of_view = pi / 4; //45 degrees
const ProjectionType projection_type = perspective;
const Vector3d camera_position(0, 0, 5);

//Maximum number of recursive calls
const int max_bounce = 5;

// Objects
std::vector<Sphere> spheres;
std::vector<Parallelogram> parallelograms;

//Material for the object, same material for all objects
const Color obj_ambient_color(0.5, 0.1, 0.1, 0);
const Color obj_diffuse_color(0.5, 0.5, 0.5, 0);
const Color obj_specular_color(0.2, 0.2, 0.2, 0);
const double obj_specular_exponent = 256.0;
const Color obj_reflection_color(0.7, 0.7, 0.7, 0);
const Color obj_refraction_color(0.7, 0.7, 0.7, 0);

// Precomputed (or otherwise) gradient vectors at each grid node
const int grid_size = 20;
std::vector<std::vector<Vector2d>> grid;

//Lights
std::vector<Light> lights;
//Ambient light
const Color ambient_light(0.2, 0.2, 0.2, 0);

//Fills the different arrays
void setup_scene() {
    grid.resize(grid_size + 1);
    for (int i = 0; i < grid_size + 1; ++i) {
        grid[i].resize(grid_size + 1);
        for (int j = 0; j < grid_size + 1; ++j)
            grid[i][j] = Vector2d::Random().normalized();
    }

    //Spheres
    spheres.emplace_back(Vector3d(10, 0, 1), 1);
    spheres.emplace_back(Vector3d(7, 0.05, -1), 1);
    spheres.emplace_back(Vector3d(4, 0.1, 1), 1);
    spheres.emplace_back(Vector3d(1, 0.2, -1), 1);
    spheres.emplace_back(Vector3d(-2, 0.4, 1), 1);
    spheres.emplace_back(Vector3d(-5, 0.8, -1), 1);
    spheres.emplace_back(Vector3d(-8, 1.6, 1), 1);

    //parallelograms
	Vector3d pgram_origin(-100, -1.25, -100), pgram_A(100, 0, -100), pgram_B(-100, -1.2, 100);
    parallelograms.emplace_back(pgram_origin, pgram_A - pgram_origin, pgram_B - pgram_origin);

    //Lights
	Color light_color(16, 16, 16, 0);
    lights.emplace_back(Vector3d(8, 8, 0), light_color);
    lights.emplace_back(Vector3d(6, -8, 0), light_color);
    lights.emplace_back(Vector3d(4, 8, 0), light_color);
    lights.emplace_back(Vector3d(2, -8, 0), light_color);
    lights.emplace_back(Vector3d(0, 8, 0), light_color);
    lights.emplace_back(Vector3d(-2, -8, 0), light_color);
    lights.emplace_back(Vector3d(-4, 8, 0), light_color);
}

//We need to make this function visible
Color shoot_ray(const Vector3d &ray_origin, const Vector3d &ray_direction, int max_bounce);

////////////////////////////////////////////////////////////////////////////////
// Perlin noise code
////////////////////////////////////////////////////////////////////////////////

// Function to linearly interpolate between a0 and a1
// Weight w should be in the range [0.0, 1.0]
inline double lerp(const double a0, const double a1, const double w) {
    assert(w >= 0);
    assert(w <= 1);
    return a0 + w * (a1 - a0);
}

// Computes the dot product of the distance and gradient vectors.
double dotGridGradient(int ix, int iy, double x, double y) {
    // Compute the distance vector
	Vector2d d(x - ix, y - iy);
    // Compute and return the dot-product
    return d.dot(grid[ix][iy]);
}

// Compute Perlin noise at coordinates x, y
double perlin(double x, double y) {
    // Determine grid cell coordinates x0, y0
    int x0 = int(x);
    int x1 = x0 + 1;
    int y0 = int(y);
    int y1 = y0 + 1;

    // Determine interpolation weights
    double sx = x - x0;
    double sy = y - y0;

    // Interpolate between grid point gradients
    double n0 = dotGridGradient(x0, y0, x, y);
    double n1 = dotGridGradient(x1, y0, x, y);

    double ix0 = lerp(n0, n1, sx);

    n0 = dotGridGradient(x0, y1, x, y);
    n1 = dotGridGradient(x1, y1, x, y);

    double ix1 = lerp(n0, n1, sx);
    double value = lerp(ix0, ix1, sy);

    return value;
}

Color procedural_texture(const double tu, const double tv, const bool use_perlin=true) {
    assert(tu >= 0);
    assert(tv >= 0);

    assert(tu <= 1);
    assert(tv <= 1);

	double color;
	if (use_perlin) {
		color = (perlin(tu * grid_size, tv * grid_size) + 1) / 2;
	}
	else { // Example for checkerboard texture
		color = (int(tu * grid_size) + int(tv * grid_size)) % 2;
	}
	return Color(0, color, 0, 0);
}

////////////////////////////////////////////////////////////////////////////////
// Intersection code
////////////////////////////////////////////////////////////////////////////////

//Finds the closest intersecting object returns its index
//In case of intersection it writes into p and N (intersection point and normals)
std::pair<int, double> find_nearest_object(const Vector3d &ray_origin, const Vector3d &ray_direction, Vector3d &p, Vector3d &N, const double t_threshold=0) {
    // Find the object in the scene that intersects the ray first
    // we store the index and the 'closest_t' to their expected values
    int closest_index = -1;
    double closest_t = std::numeric_limits<double>::max(); //closest t is "+ infinity"

    Vector3d tmp_p, tmp_N;
    for (int i = 0; i < spheres.size(); ++i) {
        //returns t and writes on tmp_p and tmp_N
        const double t = spheres[i].intersect_with_ray(ray_origin, ray_direction, tmp_p, tmp_N);
        //We have intersection
        if (t >= t_threshold) {
            //The point is before our current closest t
            if (t < closest_t) {
                closest_index = i;
                closest_t = t;
                p = tmp_p;
                N = tmp_N;
            }
        }
    }

    for (int i = 0; i < parallelograms.size(); ++i) {
        //returns t and writes on tmp_p and tmp_N
        const double t = parallelograms[i].intersect_with_ray(ray_origin, ray_direction, tmp_p, tmp_N);
        //We have intersection
        if (t >= t_threshold) {
            //The point is before our current closest t
            if (t < closest_t) {
                closest_index = spheres.size() + i;
                closest_t = t;
                p = tmp_p;
                N = tmp_N;
            }
        }
    }

    return std::pair<int, double>(closest_index, closest_t);
}

////////////////////////////////////////////////////////////////////////////////
// Raytracer code
////////////////////////////////////////////////////////////////////////////////

//Checks if the light is visible
bool is_light_visible(const Vector3d &ray_origin, const Vector3d &ray_direction, const Vector3d &light_position) {
    // Determine if the light is visible here
	Vector3d p, N;
	std::pair<int, double> nearest = find_nearest_object(ray_origin, ray_direction, p, N, eps);
	if (nearest.first < 0) // If the ray reach no object
		return true;
	else
		return nearest.second >= (light_position - ray_origin).norm(); // The reached point is farther than the light
}

Color shoot_ray(const Vector3d &ray_origin, const Vector3d &ray_direction, int max_bounce) {
    //Intersection point and normal, these are output of find_nearest_object
    Vector3d p, N;

    const int nearest_object = find_nearest_object(ray_origin, ray_direction, p, N, eps).first;

    if (nearest_object < 0) {
        // Return a transparent color
        return Color(0, 0, 0, 0);
    }

    // Ambient light contribution
    const Color ambient_color = obj_ambient_color.array() * ambient_light.array();

    // Punctual lights contribution (direct lighting)
    Color lights_color(0, 0, 0, 0);
    for (int i = 0; i < lights.size(); ++i) {
		const Light &light = lights[i];
        const Vector3d D = light.position - p;
        const Vector3d Li = D.normalized();

        // Shoot a shadow ray to determine if the light should affect the intersection point
		if (!is_light_visible(p, Li, light.position))
			continue;

        Color diff_color = obj_diffuse_color;

        if (nearest_object == 4) {
            //Compute UV coodinates for the point on the sphere
            const double x = p(0) - spheres[nearest_object].center[0];
            const double y = p(1) - spheres[nearest_object].center[1];
            const double z = p(2) - spheres[nearest_object].center[2];
            const double tu = acos(z / spheres[nearest_object].radius) / pi;
            const double tv = (pi + atan2(y, x)) / (2 * pi);

            diff_color = procedural_texture(tu, tv);
        }

        // shading parameters
		const double diffuse_weight = 1, specular_weight = 1;

        // Diffuse contribution
        const Color diffuse = diff_color * std::max(Li.dot(N), 0.0);

        // Specular contribution
		const Vector3d v = (ray_origin - p).normalized();
		const Vector3d h = (v + Li).normalized();
        const Color specular = pow(std::max(0., h.dot(N)), obj_specular_exponent) * obj_specular_color;

        // Attenuate lights according to the squared distance to the lights
        lights_color += (diffuse_weight * diffuse + specular_weight * specular).cwiseProduct(light.color) / D.squaredNorm();
    }

    Color refl_color = obj_reflection_color;
    if (nearest_object == 4) {
        refl_color = Color(0.5, 0.5, 0.5, 0);
    }
    // Compute the color of the reflected ray and add its contribution to the current point color.
    // use refl_color
    Color reflection_color(0, 0, 0, 0);
	if (max_bounce > 0) {
		reflection_color = shoot_ray(p, ray_direction - 2. * ray_direction.dot(N) * N, max_bounce - 1).cwiseProduct(refl_color);
	}

    // TODO: Compute the color of the refracted ray and add its contribution to the current point color.
    //       Make sure to check for total internal reflection before shooting a new ray.
    Color refraction_color(0, 0, 0, 0);

    // Rendering equation
    Color C = ambient_color + lights_color + reflection_color + refraction_color;

    //Set alpha to 1
    C(3) = 1;

    return C;
}

////////////////////////////////////////////////////////////////////////////////

void raytrace_scene(
	// Projection type
	const ProjectionType projection_type = projection_type,
	// Image size (h, w)
	const size_t w = 800,
	const size_t h = 400
) {
    std::cout << "Simple ray tracer." << std::endl;

    MatrixXd R = MatrixXd::Zero(w, h);
    MatrixXd G = MatrixXd::Zero(w, h);
    MatrixXd B = MatrixXd::Zero(w, h);
    MatrixXd A = MatrixXd::Zero(w, h); // Store the alpha mask

    // The camera always points in the direction -z
    // The sensor grid is at a distance 'focal_length' from the camera center,
    // and covers an viewing angle given by 'field_of_view'.
    double aspect_ratio = double(w) / double(h);
    double image_y = tan(field_of_view / 2) * focal_length; // compute the correct pixels size
    double image_x = aspect_ratio * image_y; // compute the correct pixels size

    // The pixel grid through which we shoot rays is at a distance 'focal_length'
    const Vector3d image_origin(-image_x, image_y, camera_position.z() - focal_length);
    const Vector3d x_displacement(2.0 / w * image_x, 0, 0);
    const Vector3d y_displacement(0, -2.0 / h * image_y, 0);

    for (unsigned i = 0; i < w; ++i)
        for (unsigned j = 0; j < h; ++j) {
            // TODO: Implement depth of field
            const Vector3d pixel_center = image_origin + (i + 0.5) * x_displacement + (j + 0.5) * y_displacement;

            // Prepare the ray
            Vector3d ray_origin;
            Vector3d ray_direction;

            if (projection_type == perspective) {
                // Perspective camera
				ray_origin = camera_position;
				ray_direction = (pixel_center - camera_position).normalized();
            }
            else if (projection_type == orthographic) {
                // Orthographic camera
                ray_origin = camera_position + Vector3d(pixel_center[0], pixel_center[1], 0);
                ray_direction = Vector3d(0, 0, -1);
            }

            const Color C = shoot_ray(ray_origin, ray_direction, max_bounce);
            R(i, j) = C(0);
            G(i, j) = C(1);
            B(i, j) = C(2);
            A(i, j) = C(3);
        }

    // Save to png
    write_matrix_to_png(R, G, B, A, filename);
}

////////////////////////////////////////////////////////////////////////////////

int main(int argc, char *argv[]) {
    setup_scene();

    raytrace_scene();
    return 0;
}
