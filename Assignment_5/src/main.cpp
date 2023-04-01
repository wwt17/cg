// C++ include
#include <iostream>
#include <string>
#include <vector>
#include <limits>
#include <fstream>
#include <algorithm>
#include <functional>
#include <numeric>
#include <cmath>
#include <getopt.h>

// Utilities for the Assignment
#include "raster.h"
#include <gif.h>

// Image writing library
#define STB_IMAGE_WRITE_IMPLEMENTATION // Do not include this line twice in your project!
#include "stb_image_write.h"

// Shortcut to avoid Eigen:: everywhere, DO NOT USE IN .h
using namespace Eigen;


enum ProjectionType {orthographic, perspective};
const char *ProjectionTypeString[] = {"orthographic", "perspective"};

////////////////////////////////////////////////////////////////////////////////
// Scene setup, global variables
////////////////////////////////////////////////////////////////////////////////
const std::string data_dir = DATA_DIR;
std::string mesh_filename(data_dir + "bunny.off");
bool refit_mesh = true;

// Triangle Mesh
MatrixX3d vertices; // n x 3 matrix (n points)
MatrixX3i facets;   // m x 3 matrix (m triangles)

//Material for the object, same material for all objects
const Color obj_ambient_color(0.0, 0.5, 0.0, 0);
const Color obj_diffuse_color(0.5, 0.5, 0.5, 0);
const Color obj_specular_color(0.2, 0.2, 0.2, 0);
double obj_specular_exponent = 256.0;
const Color obj_reflection_color(0.7, 0.7, 0.7, 0);

//Lights
std::vector<Light> lights;
//Ambient light
const Color ambient_light(0.2, 0.2, 0.2, 0);

double alpha = -1;  // not set

//Fills the different arrays
void setup_scene() {
    //Loads file
    std::ifstream in(mesh_filename);
    std::string token;
    in >> token;
    int nv, nf, ne;
    in >> nv >> nf >> ne;
    vertices.resize(nv, 3);
    facets.resize(nf, 3);
    for (int i = 0; i < nv; ++i) {
        in >> vertices(i, 0) >> vertices(i, 1) >> vertices(i, 2);
    }
    for (int i = 0; i < nf; ++i) {
        int s;
        in >> s >> facets(i, 0) >> facets(i, 1) >> facets(i, 2);
        assert(s == 3);
    }

	//Refit the mesh into [-1, +1]^3
	if (refit_mesh) {
		Vector3d center(0, 0, 0);
		for (size_t i = 0; i < nf; ++i) {
			Vector3d triangle_center(0, 0, 0);
			for (size_t k = 0; k < 3; ++k)
				triangle_center += vertices.row(facets(i, k)).transpose();
			triangle_center /= 3;
			center += triangle_center;
		}
		center /= nf;
		vertices.rowwise() -= center.transpose();  // centering
		double rad = 0;
		for (size_t i = 0; i < nv; ++i)
			rad = std::max(rad, vertices.row(i).norm());
		vertices /= rad;  // rescaling
	}

    //Lights
	lights.emplace_back(Vector3d( 8,  8, 0), Color(16, 16, 16, 0));
	lights.emplace_back(Vector3d( 6, -8, 0), Color(16, 16, 16, 0));
	lights.emplace_back(Vector3d( 4,  8, 0), Color(16, 16, 16, 0));
	lights.emplace_back(Vector3d( 2, -8, 0), Color(16, 16, 16, 0));
	lights.emplace_back(Vector3d( 0,  8, 0), Color(16, 16, 16, 0));
	lights.emplace_back(Vector3d(-2, -8, 0), Color(16, 16, 16, 0));
	lights.emplace_back(Vector3d(-4,  8, 0), Color(16, 16, 16, 0));
}

////////////////////////////////////////////////////////////////////////////////

void render_scene(
	// Shading
	const std::string shading,
	// Transformation
	const bool transformation,
	const int timesteps,
	const int delay,
	// Camera settings
	const Vector3d camera_position,
	const ProjectionType projection_type,
	const double focal_length,
	const double field_of_view,
	// Image size (h, w)
	const size_t w,
	const size_t h,
	const std::string filename
) {
    std::cout << "Simple rendering system for triangulated 3D models based on rasterization." << std::endl;

	// TODO: will only be used in perspective camera
    // The camera always points in the direction -z
    // The sensor grid is at a distance 'focal_length' from the camera center,
    // and covers an viewing angle given by 'field_of_view'.
    double aspect_ratio = double(w) / double(h);
    double image_y = tan(field_of_view / 2) * focal_length;
    double image_x = aspect_ratio * image_y;

	// The Framebuffer storing the image rendered by the rasterizer
	Eigen::Matrix<FrameBufferAttributes,Eigen::Dynamic,Eigen::Dynamic> frameBuffer(w, h);

	UniformAttributes uniform;

	Program program;
	program.VertexShader = [](const VertexAttributes& va, const UniformAttributes& uniform) {
		VertexAttributes new_va = va.transform(uniform.affine, uniform.rotation);
		new_va.color = uniform.color;
		return new_va;
	};
	program.FragmentShader = [](const VertexAttributes& va, const UniformAttributes& uniform) {
		return FragmentAttributes(va.color, va.position);
	};
	program.BlendingShader = [](const FragmentAttributes& fa, const FrameBufferAttributes& previous) {
		FrameBufferAttributes new_fba = previous;
		new_fba.blend_with(fa.color, fa.position[2]);
		return new_fba;
	};

	Program shading_program = program;

	const size_t nv = vertices.rows(), nf = facets.rows();

	// Triangles
	std::vector<VertexAttributes> triangle_vertices;
	for (size_t i = 0; i < nf; ++i) {
		Vector3d points[3];
		for (size_t k = 0; k < 3; ++k)
			points[k] = vertices.row(facets(i, k));
		Vector3d normal = (points[1] - points[0]).cross(points[2] - points[0]).normalized();
		for (size_t k = 0; k < 3; ++k)
			triangle_vertices.push_back(VertexAttributes(position3_to_position(points[k]), normal));
	}

	// Edges
	std::vector<VertexAttributes> edge_vertices;
	for (size_t i = 0; i < nf; ++i) {
		for (size_t k = 0; k < 3; ++k) {
			for (size_t dk = 0; dk < 2; ++dk) {
				const Vector3d position = vertices.row(facets(i, (k + dk) % 3));
				edge_vertices.push_back(VertexAttributes(position3_to_position(position)));
			}
		}
	}

	std::function<void()> rasterize;

	if (shading == "silhouette") {
		if (alpha == -1)
			alpha = 1;
		rasterize = [&]() {
			uniform.color = Color(1, 1, 1, alpha);
			rasterize_triangles(program, uniform, triangle_vertices, frameBuffer);
		};
	}
	else if (shading == "wireframe") {
		if (alpha == -1)
			alpha = 0.5;
		rasterize = [&]() {
			uniform.color = Color(1, 1, 1, alpha);
			rasterize_triangles(program, uniform, triangle_vertices, frameBuffer);
			uniform.color = Color(0, 0, 0, 1);
			rasterize_lines(program, uniform, edge_vertices, 0.5, frameBuffer);
		};
	}
	else {
		if (alpha == -1)
			alpha = 1;
		uniform.obj_ambient_color = obj_ambient_color;
		uniform.obj_diffuse_color = obj_diffuse_color;
		uniform.obj_specular_color = obj_specular_color;
		uniform.obj_specular_exponent = obj_specular_exponent;
		uniform.ambient_light = ambient_light;
		uniform.lights = lights;
		uniform.alpha = alpha;
		shading_program.VertexShader = [](const VertexAttributes& va, const UniformAttributes& uniform) {
			VertexAttributes new_va = va.transform(uniform.affine, uniform.rotation);

			const Vector3d ray_direction(0, 0, 1);
			const Vector3d &p = position_to_position3(new_va.position), &N = new_va.normal;

			// Ambient light contribution
			const Color ambient_color = uniform.obj_ambient_color.array() * uniform.ambient_light.array();

			// Punctual lights contribution (direct lighting)
			Color lights_color(0, 0, 0, 0);
			for (const Light &light: uniform.lights) {
				const Vector3d D = light.position - p;

				// Diffuse contribution
				const Vector3d Li = D.normalized();
				const Color diffuse = uniform.obj_diffuse_color * std::max(Li.dot(N), 0.0);

				// Specular contribution
				const Vector3d Hi = (Li - ray_direction).normalized();
				const Color specular = uniform.obj_specular_color * std::pow(std::max(N.dot(Hi), 0.0), uniform.obj_specular_exponent);

				// Attenuate lights according to the squared distance to the lights
				lights_color += (diffuse + specular).cwiseProduct(light.color) / D.squaredNorm();
			}

			// Rendering equation
			new_va.color = ambient_color + lights_color;

			//Set alpha
			new_va.color(3) = uniform.alpha;

			return new_va;
		};
		if (shading == "per-vertex") {
			// Compute the per-vertex normals by averaging per-face normals
			MatrixX3d sum_normals = MatrixX3d::Zero(nv, 3);
			VectorXi cnt_facets = VectorXi::Zero(nv);
			auto it = triangle_vertices.begin();
			for (size_t i = 0; i < nf; ++i) {
				for (size_t k = 0; k < 3; ++k) {
					size_t vertex_id = facets(i, k);
					cnt_facets(vertex_id) += 1;
					sum_normals.row(vertex_id) += it->normal;
					++it;
				}
			}
			MatrixX3d normals = sum_normals.array().colwise() / cnt_facets.array().cast<double>();
			it = triangle_vertices.begin();
			for (size_t i = 0; i < nf; ++i) {
				for (size_t k = 0; k < 3; ++k) {
					size_t vertex_id = facets(i, k);
					it->normal = normals.row(vertex_id);
					++it;
				}
			}
			rasterize = [&]() {
				rasterize_triangles(shading_program, uniform, triangle_vertices, frameBuffer);
			};
		}
		else {
			rasterize = [&]() {
				rasterize_triangles(shading_program, uniform, triangle_vertices, frameBuffer);
				uniform.color = Color(0, 0, 0, 1);
				rasterize_lines(program, uniform, edge_vertices, 0.5, frameBuffer);
			};
		}
	}

	Matrix4d view = Matrix4d::Identity();
	if (aspect_ratio < 1)
		view(1, 1) = aspect_ratio;
	else
		view(0, 0) = 1 / aspect_ratio;
	uniform.affine = view * uniform.affine;

	std::vector<uint8_t> image;
	if (transformation) {
		GifWriter g;
		GifBegin(&g, filename.c_str(), frameBuffer.rows(), frameBuffer.cols(), delay);

		for (int t = 0; t <= timesteps; t++) {
			// Rotation
			double theta = 2 * pi * t / timesteps;
			uniform.rotation <<
				cos(theta), 0, -sin(theta), 0,
				         0, 1,           0, 0,
				sin(theta), 0,  cos(theta), 0,
				         0, 0,           0, 1;
			// Affine transformation
			uniform.affine = Matrix4d::Identity();
			uniform.affine.block<3, 1>(0, 3) = (double)(t - timesteps) / timesteps * Vector3d(0, 0.7, -0.7);
			uniform.affine = view * uniform.affine;

			frameBuffer.setConstant(FrameBufferAttributes());
			rasterize();
			framebuffer_to_uint8(frameBuffer, image);
			GifWriteFrame(&g, image.data(), frameBuffer.rows(), frameBuffer.cols(), delay);
		}

		GifEnd(&g);
	}
	else {
		rasterize();
		framebuffer_to_uint8(frameBuffer, image);
		stbi_write_png(filename.c_str(), frameBuffer.rows(), frameBuffer.cols(), 4, image.data(), frameBuffer.rows()*4);
	}
}

////////////////////////////////////////////////////////////////////////////////

int main(int argc, char *argv[]) {
	std::ios::sync_with_stdio(false);

	static struct option long_options[] = {
		{"filename",              required_argument, 0, 0  },
		{"mesh_filename",         required_argument, 0, 6  },
		{"no_refit_mesh",         no_argument,       0, 7  },
		{"shading",               required_argument, 0, 8  },
		{"alpha",                 required_argument, 0, 'a'},
		{"transformation",        no_argument,       0, 't'},
		{"timesteps",             required_argument, 0, 9  },
		{"delay",                 required_argument, 0, 10 },
		{"focal_length",          required_argument, 0, 'f'},
		{"field_of_view",         required_argument, 0, 1  },
		{"projection_type",       required_argument, 0, 'p'},
		{"obj_specular_exponent", required_argument, 0, 2  },
		{"grid_size",             required_argument, 0, 'g'},
		{0,                       0,                 0, 0  }
	};

	//Camera settings
	const Vector3d camera_position(0, 0, 2);
	ProjectionType projection_type = perspective;
	double focal_length = 2;
	double field_of_view = pi / 4; // 45 degrees
	size_t w = 500, h = 500;
	std::string filename("triangle.png");
	std::string shading("wireframe");
	bool transformation = false;
	int timesteps = 20;
	int delay = 25;

	int opt;
	while ((opt = getopt_long(argc, argv, "w:h:f:p:a:", long_options, NULL)) != -1) {
		switch (opt) {
		case 0:
			filename = optarg;
			break;
		case 6:
			mesh_filename = optarg;
			break;
		case 7:
			refit_mesh = false;
			break;
		case 8:
			shading = optarg;
			break;
		case 'a':
			alpha = atof(optarg);
			break;
		case 't':
			transformation = true;
			break;
		case 9:
			timesteps = atoi(optarg);
			break;
		case 10:
			delay = atoi(optarg);
			break;
		case 1:
			field_of_view = pi / 180 * atoi(optarg); // field of view in degree
			break;
		case 2:
			obj_specular_exponent = atof(optarg);
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
		default: // '?'
			std::cerr << "Read the code for usage." << std::endl;
			exit(EXIT_FAILURE);
		}
	}

	if (transformation) {
		if (!endsWith(filename, ".gif")) {
			std::cerr << "With object transformation, output filename should ends with .gif." << std::endl;
			exit(EXIT_FAILURE);
		}
	}
	else {
		if (!endsWith(filename, ".png")) {
			std::cerr << "Without object transformation, output filename should ends with .png." << std::endl;
			exit(EXIT_FAILURE);
		}
	}

    setup_scene();

    render_scene(
		shading,
		transformation,
		timesteps,
		delay,
		camera_position,
		projection_type,
		focal_length,
		field_of_view,
		w, h,
		filename);

    return 0;
}
