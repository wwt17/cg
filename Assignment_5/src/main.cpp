// C++ include
#include <iostream>
#include <string>
#include <vector>
#include <limits>
#include <fstream>
#include <algorithm>
#include <numeric>
#include <cmath>
#include <getopt.h>

// Utilities for the Assignment
#include "raster.h"

// Image writing library
#define STB_IMAGE_WRITE_IMPLEMENTATION // Do not include this line twice in your project!
#include "stb_image_write.h"

// Shortcut to avoid Eigen:: everywhere, DO NOT USE IN .h
using namespace Eigen;

////////////////////////////////////////////////////////////////////////////////
// Classes for the scene
////////////////////////////////////////////////////////////////////////////////

template<class T>
T sqr(const T x) {
	return x * x;
}


enum ProjectionType {orthographic, perspective};
const char *ProjectionTypeString[] = {"orthographic", "perspective"};


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
		Matrix<double, 1, 3>
			min_coord = vertices.colwise().minCoeff(),
			max_coord = vertices.colwise().maxCoeff(),
			mid_coord = (min_coord + max_coord) / 2;
		double resize_length = (max_coord - min_coord).maxCoeff();
		vertices.rowwise() -= mid_coord;  // centering
		vertices /= resize_length / 2;  // rescaling
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
	// Camera settings
	Vector3d camera_position,
	ProjectionType projection_type,
	double focal_length,
	double field_of_view,
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
		return VertexAttributes(va.position, uniform.color);
	};
	program.FragmentShader = [](const VertexAttributes& va, const UniformAttributes& uniform) {
		return FragmentAttributes(va.color, va.position);
	};
	program.BlendingShader = [](const FragmentAttributes& fa, const FrameBufferAttributes& previous) {
		FrameBufferAttributes new_fba = previous;
		new_fba.blend_with(fa.color, fa.position[2]);
		return new_fba;
	};

	// Triangles
	std::vector<VertexAttributes> triangle_vertices;
	const size_t nf = facets.rows();
	for (size_t i = 0; i < nf; ++i) {
		for (size_t k = 0; k < 3; ++k) {
			const Vector3d position = vertices.row(facets(i, k));
			triangle_vertices.push_back(VertexAttributes(position3_to_position(position)));
		}
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

	if (shading == "silhouette") {
		uniform.color = Color(1, 1, 1, 1);
		rasterize_triangles(program, uniform, triangle_vertices, frameBuffer);
	}
	else if (shading == "wireframe") {
		uniform.color = Color(1, 1, 1, 0.5);
		rasterize_triangles(program, uniform, triangle_vertices, frameBuffer);
		uniform.color = Color(0, 0, 0, 1);
		rasterize_lines(program, uniform, edge_vertices, 0.5, frameBuffer);
	}

	std::vector<uint8_t> image;
	framebuffer_to_uint8(frameBuffer,image);
	stbi_write_png(filename.c_str(), frameBuffer.rows(), frameBuffer.cols(), 4, image.data(), frameBuffer.rows()*4);
}

////////////////////////////////////////////////////////////////////////////////

int main(int argc, char *argv[]) {
	std::ios::sync_with_stdio(false);

	static struct option long_options[] = {
		{"filename",              required_argument, 0, 0  },
		{"mesh_filename",         required_argument, 0, 6  },
		{"no_refit_mesh",         no_argument,       0, 7  },
		{"shading",               required_argument, 0, 8  },
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

	int opt;
	while ((opt = getopt_long(argc, argv, "w:h:f:p:", long_options, NULL)) != -1) {
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

    setup_scene();

    render_scene(
		shading,
		camera_position,
		projection_type,
		focal_length,
		field_of_view,
		w, h,
		filename);

    return 0;
}
