////////////////////////////////////////////////////////////////////////////////
// C++ include
#include <iostream>
#include <string>
#include <vector>
#include <limits>
#include <cmath>
#include <cstdlib>
#include <random>
#include <getopt.h>

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
			if (t <= eps) {  // if the first intersection is at the back
				// try the second intersection
				t = (-b_over_2 + std::sqrt(delta)) / a;
			}
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


std::random_device rd;
std::mt19937 gen(rd());


class Sampler {
public:
	double stddev;
	Sampler(const double stddev): stddev(stddev) {
	}
	Vector3d operator ()(const Vector3d &camera_position) const {
		if (stddev <= 0)
			return camera_position;

		std::normal_distribution<> d(0., stddev);
		const Vector3d deviation(d(gen), d(gen), 0);
		return camera_position + deviation;
	}
};


////////////////////////////////////////////////////////////////////////////////
// Scene setup, global variables
////////////////////////////////////////////////////////////////////////////////

// Objects
std::vector<Sphere> spheres;
std::vector<Parallelogram> parallelograms;

//Material for the object, same material for all objects
const Color obj_ambient_color(0.5, 0.1, 0.1, 0);
const Color obj_diffuse_color(0.5, 0.5, 0.5, 0);
const Color obj_specular_color(0.2, 0.2, 0.2, 0);
double obj_specular_exponent = 256.0;
// refractive indexes
const double eta_0 = 1;  // air
double eta = 1e18;       // object; very large value for no refraction
double s_fraction = 1;  // fraction of s-polarized light (in s-polarized and p-polarized light)
const Color obj_reflection_color(0.7, 0.7, 0.7, 0);
const Color obj_refraction_color(0.7, 0.7, 0.7, 0);

// Precomputed (or otherwise) gradient vectors at each grid node
int grid_size = 20;
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
	Vector3d pgram_origin(-100, -1.25, -100), pgram_A(-100, -1.2, 100), pgram_B(100, 0, -100);
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

inline double refract(
	const Vector3d &N,   // norm. points out.
	const Vector3d &l,   // the light or view vector; these two cases are the same
	// See the figure in the README
	double &cos_theta_1,
	double &sin_theta_1,
	double &cos_theta_2,
	double &sin_theta_2
) { // Compute refraction, put computed things in the references, and return reflectance
	cos_theta_1 = std::max(-1., std::min(1., l.dot(N)));
	const bool out_light = cos_theta_1 >= 0;
	double eta_1, eta_2;
	if (out_light) {
		eta_1 = eta_0;
		eta_2 = eta;
	}
	else {
		eta_1 = eta;
		eta_2 = eta_0;
	}
	sin_theta_1 = sqrt(1 - sqr(cos_theta_1));
	sin_theta_2 = eta_1 / eta_2 * sin_theta_1;  // Snell-Descartes law
	double reflectance;
	if (sin_theta_2 >= 1) {  // total internal reflection
		reflectance = 1;
	}
	else {
		cos_theta_2 = (out_light ? 1 : -1) * sqrt(1 - sqr(sin_theta_2));  // make it having the same sign as cos_theta_1
		// Fresnel equations
		const double reflectance_s = sqr((eta_1 * cos_theta_1 - eta_2 * cos_theta_2) / (eta_1 * cos_theta_1 + eta_2 * cos_theta_2)),
					 reflectance_p = sqr((eta_1 * cos_theta_2 - eta_2 * cos_theta_1) / (eta_1 * cos_theta_2 + eta_2 * cos_theta_1));
		reflectance = s_fraction * reflectance_s + (1 - s_fraction) * reflectance_p;
		cos_theta_2 *= -1;  // flip the sign
	}
	return reflectance;
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

		double cos_theta_1, sin_theta_1, cos_theta_2, sin_theta_2;
		const double reflectance = refract(N, Li, cos_theta_1, sin_theta_1, cos_theta_2, sin_theta_2);
		const bool out_light = cos_theta_1 >= 0;
		Vector3d Li_refracted_flipped;
		if (sin_theta_2 < 1) {
			Vector3d N_T = Li.cross(N).cross(N);
			if (N_T.norm() != 0)
				N_T.normalize();
			Li_refracted_flipped = -sin_theta_2 * N_T + cos_theta_2 * N;
		}

		const Vector3d v = -ray_direction.normalized();
		const bool out_view = v.dot(N) >= 0;
		const bool is_reflection = out_view == out_light;

		// Diffuse contribution
		const Color diffuse = diff_color * std::abs(cos_theta_1);

		// Specular contribution
		const Vector3d h = (v + (is_reflection ? Li : Li_refracted_flipped)).normalized();
		const Color specular = pow(std::abs(h.dot(N)), obj_specular_exponent) * obj_specular_color;

        // Attenuate lights according to the squared distance to the lights
		// Also ration the amount of light by whether it is reflection/refraction
        lights_color += (diffuse_weight * diffuse + specular_weight * specular).cwiseProduct(light.color) / D.squaredNorm() * (is_reflection ? reflectance : 1 - reflectance);
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

    // Compute the color of the refracted ray and add its contribution to the current point color.
    // Make sure to check for total internal reflection before shooting a new ray.
    Color refraction_color(0, 0, 0, 0);
	if (max_bounce > 0) {
		//refraction_color = 
	}

    // Rendering equation
    Color C = ambient_color + lights_color + reflection_color + refraction_color;

    //Set alpha to 1
    C(3) = 1;

    return C;
}

////////////////////////////////////////////////////////////////////////////////

void raytrace_scene(
	// Camera settings
	Vector3d camera_position,
	ProjectionType projection_type,
	double focal_length,
	double field_of_view,
	// Image size (h, w)
	const size_t w,
	const size_t h,
	//Maximum number of recursive calls
	int max_bounce,
	// Sample the deviation from the camera position
	const Sampler sampler,
	const unsigned n_sample,
	const std::string filename
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
            const Vector3d pixel_center = image_origin + (i + 0.5) * x_displacement + (j + 0.5) * y_displacement;

            // Prepare the ray
            Vector3d ray_origin;

            if (projection_type == perspective) {  // Perspective camera
				ray_origin = camera_position;
            }
            else if (projection_type == orthographic) {  // Orthographic camera
                ray_origin = camera_position + Vector3d(pixel_center[0], pixel_center[1], 0);
            }

            // Implement depth of field
			// Sample multiple ray origins and average the colors
            Color C(0, 0, 0, 0);
			for (unsigned i_sample = 0; i_sample < n_sample; ++i_sample) {
				Vector3d ray_origin_sample = sampler(ray_origin);
				Vector3d ray_direction = (pixel_center - ray_origin_sample).normalized();
				Color sample_color = shoot_ray(ray_origin_sample, ray_direction, max_bounce);
				C += sample_color;
			}
			C /= n_sample;

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
	static struct option long_options[] = {
		{"filename",              required_argument, 0, 0  },
		{"focal_length",          required_argument, 0, 'f'},
		{"field_of_view",         required_argument, 0, 1  },
		{"projection_type",       required_argument, 0, 'p'},
		{"obj_specular_exponent", required_argument, 0, 2  },
		{"eta",                   required_argument, 0, 'e'},
		{"s_fraction",            required_argument, 0, 5  },
		{"max_bounce",            required_argument, 0, 'b'},
		{"grid_size",             required_argument, 0, 'g'},
		{"stddev",                required_argument, 0, 3  },
		{"n_sample",              required_argument, 0, 4  },
		{0,                       0,                 0, 0  }
	};

	//Camera settings
	const Vector3d camera_position(0, 0, 5);
	ProjectionType projection_type = perspective;
	double focal_length = 10;
	double field_of_view = pi / 4; // 45 degrees
	size_t w = 800, h = 400;
	//Maximum number of recursive calls
	int max_bounce = 5;
	//Sample
	Sampler sampler(0);
	unsigned n_sample = 1;
	std::string filename("raytrace.png");

	int opt;
	while ((opt = getopt_long(argc, argv, "w:h:f:p:e:b:g:", long_options, NULL)) != -1) {
		switch (opt) {
		case 0:
			filename = optarg;
			break;
		case 1:
			field_of_view = pi / 180 * atoi(optarg); // field of view in degree
			break;
		case 2:
			obj_specular_exponent = atof(optarg);
			break;
		case 3:
			sampler = Sampler(atof(optarg));
			break;
		case 4:
			n_sample = atoi(optarg);
			break;
		case 5:
			s_fraction = atof(optarg);
			break;
		case 'w':
			w = atoi(optarg);
			break;
		case 'h':
			h = atoi(optarg);
			break;
		case 'f':
			focal_length = atof(optarg);
			break;
		case 'p':
			if (optarg[0] == 'p') {
				projection_type = perspective;
			}
			else if (optarg[0] == 'o') {
				projection_type = orthographic;
			}
			else {
				std::cerr << "Unknown projection type: " << optarg << std::endl;
				exit(EXIT_FAILURE);
			}
			break;
		case 'e':
			eta = atof(optarg);
			break;
		case 'b':
			max_bounce = atoi(optarg);
			break;
		case 'g':
			grid_size = atoi(optarg);
			break;
		default: // '?'
			std::cerr << "Read the code for usage." << std::endl;
			exit(EXIT_FAILURE);
		}
	}

    setup_scene();

    raytrace_scene(
		camera_position,
		projection_type,
		focal_length,
		field_of_view,
		w, h,
		max_bounce,
		sampler, n_sample,
		filename);

    return 0;
}
