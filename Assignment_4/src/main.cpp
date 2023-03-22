////////////////////////////////////////////////////////////////////////////////
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
#include "utils.h"

// Image writing library
#define STB_IMAGE_WRITE_IMPLEMENTATION // Do not include this line twice in your project!
#include "stb_image_write.h"

// Shortcut to avoid Eigen:: everywhere, DO NOT USE IN .h
using namespace Eigen;

////////////////////////////////////////////////////////////////////////////////
// Class to store tree
////////////////////////////////////////////////////////////////////////////////
class AABBTree {
public:
    class Node {
    public:
        AlignedBox3d bbox;
        int parent;   // Index of the parent node (-1 for root)
        int left;     // Index of the left child (-1 for a leaf)
        int right;    // Index of the right child (-1 for a leaf)
        int triangle; // Index of the node triangle (-1 for internal nodes)
    };

    std::vector<Node> nodes;
	int n_nodes;
    int root;

    AABBTree() = default;                             // Default empty constructor
    AABBTree(const MatrixX3d &V, const MatrixX3i &F); // Build a BVH from an existing mesh

	int top_down_construct(int* const st, int* const ed, const MatrixX3d &V, const MatrixX3i &F, const MatrixX3d &centroids);
	int bottom_up_construct(const MatrixX3d &V, const MatrixX3i &F, double (* const criterion)(const Node &n1, const Node &n2));
	void find_nearest_object(const Vector3d &ray_origin, const Vector3d &ray_direction, Vector3d &p, Vector3d &N, std::pair<int, double> &nearest, const int cur_node);
};

////////////////////////////////////////////////////////////////////////////////
// Classes for the scene
////////////////////////////////////////////////////////////////////////////////

typedef Vector4d Color; // RBGA Color (in [0, 1])

const double eps = 1e-7, inf = 1./0., pi = acos(-1);

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
std::string mesh_filename(data_dir + "dodeca.off");

// Triangle Mesh
MatrixX3d vertices; // n x 3 matrix (n points)
MatrixX3i facets;   // m x 3 matrix (m triangles)
AABBTree bvh;
bool use_bvh = true;

//Material for the object, same material for all objects
const Color obj_ambient_color(0.0, 0.5, 0.0, 0);
const Color obj_diffuse_color(0.5, 0.5, 0.5, 0);
const Color obj_specular_color(0.2, 0.2, 0.2, 0);
double obj_specular_exponent = 256.0;
const Color obj_reflection_color(0.7, 0.7, 0.7, 0);

// Precomputed (or otherwise) gradient vectors at each grid node
int grid_size = 20;
std::vector<std::vector<Vector2d>> grid;

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

    //setup tree
	if (use_bvh) {
		bvh = AABBTree(vertices, facets);
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
// BVH Code
////////////////////////////////////////////////////////////////////////////////

double squared_centroid_distance(const AABBTree::Node &n1, const AABBTree::Node &n2) {
	return (n1.bbox.center() - n2.bbox.center()).squaredNorm();
}

double volume_increase(const AABBTree::Node &n1, const AABBTree::Node &n2) {
	AlignedBox3d bbox = n1.bbox;
	bbox.extend(n2.bbox);
	return bbox.volume() - n1.bbox.volume() - n2.bbox.volume();
}

double (*criterion)(const AABBTree::Node &n1, const AABBTree::Node &n2) = &squared_centroid_distance;

bool top_down_construction = true;
AABBTree::AABBTree(const MatrixX3d &V, const MatrixX3i &F) {
	n_nodes = 0;
	if (F.rows() == 0) {
		root = -1;
	}
	else {
		nodes.resize(F.rows() * 2 - 1);
		if (top_down_construction) {
			// Compute the centroids of all the triangles in the input mesh
			MatrixX3d centroids(F.rows(), V.cols());
			centroids.setZero();
			for (int i = 0; i < F.rows(); ++i) {
				for (int k = 0; k < F.cols(); ++k) {
					centroids.row(i) += V.row(F(i, k));
				}
				centroids.row(i) /= F.cols();
			}

			int triangle_indexes[F.rows()];
			for (int i = 0; i < F.rows(); ++i)
				triangle_indexes[i] = i;

			root = top_down_construct(triangle_indexes, triangle_indexes + F.rows(), V, F, centroids);
		}
		else {
			root = bottom_up_construct(V, F, criterion);
		}
	}
}

inline void build_leaf_node(AABBTree::Node &node, const int triangle, const MatrixX3d &V, const MatrixX3i &F) {
	node.bbox = AlignedBox3d();
	for (int k = 0; k < 3; ++k) {
		node.bbox.extend(V.row(F(triangle, k)).transpose());
	}
	node.left = -1;
	node.right = -1;
	node.triangle = triangle;
}

inline void build_internal_node(std::vector<AABBTree::Node> &nodes, const int cur_node_id, const int left, const int right) {
	AABBTree::Node &cur_node = nodes[cur_node_id];
	cur_node.bbox = AlignedBox3d();
	cur_node.left = left;
	cur_node.bbox.extend(nodes[left ].bbox);
	nodes[left ].parent = cur_node_id;
	cur_node.right = right;
	cur_node.bbox.extend(nodes[right].bbox);
	nodes[right].parent = cur_node_id;
	cur_node.triangle = -1;
}

bool centroid_sorting_split = false;
int AABBTree::top_down_construct(int* const st, int* const ed, const MatrixX3d &V, const MatrixX3i &F, const MatrixX3d &centroids) {
	const int cur_node_id = n_nodes++;
	nodes[cur_node_id].parent = -1;
	if (ed - st == 1) {  // leaf: only 1
		build_leaf_node(nodes[cur_node_id], *st, V, F);
	}
	else {
		// Split each set of primitives into 2 sets of roughly equal size,
		// based on sorting the centroids along one direction or another.
		if (centroid_sorting_split) {
			AlignedBox3d centroid_bbox;
			for (auto it = st; it < ed; ++it) {
				const int triangle = *it;
				centroid_bbox.extend(centroids.row(triangle).transpose());
			}
			const Vector3d span_size = centroid_bbox.max() - centroid_bbox.min();
			int max_span_dim = 0;
			for (int k = 0; k < 3; ++k)
				if (span_size(k) > span_size(max_span_dim))
					max_span_dim = k;
			std::nth_element(
				st, st + (ed - st) / 2, ed,
				[=](const int a, const int b) {
					return centroids(a, max_span_dim) < centroids(b, max_span_dim);
				}
			);
		}
		auto mid = st + ((ed - st + 1) / 2);
		build_internal_node(
			nodes,
			cur_node_id,
			top_down_construct(st, mid, V, F, centroids),
			top_down_construct(mid, ed, V, F, centroids)
		);
	}
	return cur_node_id;
}

int AABBTree::bottom_up_construct(const MatrixX3d &V, const MatrixX3i &F, double (* const criterion)(const AABBTree::Node &n1, const AABBTree::Node &n2)) {
	int cur_set_size = 0;
	int cur_set[F.rows()];  // containing node indices
	// Initialize cur_set with all leaf nodes created from every triangle
	for (; n_nodes < F.rows(); ++n_nodes) {
		build_leaf_node(nodes[n_nodes], n_nodes, V, F);
		cur_set[cur_set_size++] = n_nodes;
	}

	while (cur_set_size > 1) {
		// find i_to_merge, j_to_merge with min criterion value
		double min_criterion_value = inf;
		int i_to_merge = 0, j_to_merge = 1;
		for (int i = 0; i < cur_set_size; ++i)
			for (int j = i + 1; j < cur_set_size; ++j) {
				double criterion_value = criterion(nodes[cur_set[i]], nodes[cur_set[j]]);
				if (criterion_value < min_criterion_value) {
					min_criterion_value = criterion_value;
					i_to_merge = i;
					j_to_merge = j;
				}
			}

		// merge cur_set[i_to_merge], cur_set[j_to_merge]
		const int cur_node_id = n_nodes++;
		build_internal_node(
			nodes,
			cur_node_id,
			cur_set[i_to_merge],
			cur_set[j_to_merge]
		);

		// update cur_set
		cur_set[i_to_merge] = cur_node_id;
		cur_set[j_to_merge] = cur_set[--cur_set_size];
	}

	const int root = cur_set[0];
	nodes[root].parent = -1;
	return root;
}


////////////////////////////////////////////////////////////////////////////////
// Intersection code
////////////////////////////////////////////////////////////////////////////////

double ray_triangle_intersection(const Vector3d &ray_origin, const Vector3d &ray_direction, const Vector3d &a, const Vector3d &b, const Vector3d &c, Vector3d &p, Vector3d &N) {
    // Compute whether the ray intersects the given triangle.

	Matrix4d M; M << a, b, c, -ray_direction,
	                 1, 1, 1,              0;
	Vector4d target; target << ray_origin, 1;
	Vector4d coord = M.inverse() * target;
	p = coord(0) * a + coord(1) * b + coord(2) * c;

	N = (b - a).cross(c - a).normalized();
	if (N.dot(-ray_direction) < 0)
		N = -N;

    return coord(0) >= 0 && coord(1) >= 0 && coord(2) >= 0 ? coord(3) : -inf;
}

double ray_box_intersection(const Vector3d &ray_origin, const Vector3d &ray_direction, const AlignedBox3d &box) {
    // Compute whether the ray intersects the given box.
    // we are not testing with the real surface here anyway.
	// Intersection of t ranges in all dimensions
	double box_t_l = -inf, box_t_r = +inf;
	for (int j = 0; j < 3; ++j) {
		double t_l, t_r;
		if (ray_direction(j) == 0) {
			t_l = box.min()(j) <= ray_origin(j) ? -inf : +inf;
			t_r = box.max()(j) >= ray_origin(j) ? +inf : -inf;
		}
		else {
			t_l = (box.min()(j) - ray_origin(j)) / ray_direction(j);
			t_r = (box.max()(j) - ray_origin(j)) / ray_direction(j);
			if (t_l > t_r) std::swap(t_l, t_r);
		}
		box_t_l = std::max(box_t_l, t_l);
		box_t_r = std::min(box_t_r, t_r);
	}
	return box_t_l > box_t_r ? -inf : box_t_r >= 0 ? std::max(box_t_l, 0.) : box_t_r;
}

inline void update_nearest_triangle(const Vector3d &ray_origin, const Vector3d &ray_direction, Vector3d &p, Vector3d &N, std::pair<int, double> &nearest, const size_t i) {
	Vector3d tmp_p, tmp_N;
	const double t = ray_triangle_intersection(
		ray_origin, ray_direction,
		vertices.row(facets(i, 0)),
		vertices.row(facets(i, 1)),
		vertices.row(facets(i, 2)),
		tmp_p, tmp_N);
	if (t >= eps && t < nearest.second) {
		nearest.first = i;
		nearest.second = t;
		p = tmp_p;
		N = tmp_N;
	}
}

//Finds the closest intersecting object returns its index
//In case of intersection it writes into p and N (intersection point and normals)
std::pair<int, double> find_nearest_object(const Vector3d &ray_origin, const Vector3d &ray_direction, Vector3d &p, Vector3d &N) {
	std::pair<int, double> nearest(-1, inf);  // index of nearest (object, t)
	if (!use_bvh) {
		// Method (1): Traverse every triangle and return the closest hit.
		const size_t nf = facets.rows();
		for (size_t i = 0; i < nf; ++i) {
			update_nearest_triangle(ray_origin, ray_direction, p, N, nearest, i);
		}
	}
	else {
		// Method (2): Traverse the BVH tree and test the intersection with a
		// triangles at the leaf nodes that intersects the input ray.
		bvh.find_nearest_object(ray_origin, ray_direction, p, N, nearest, bvh.root);
	}
	return nearest;
}

void AABBTree::find_nearest_object(const Vector3d &ray_origin, const Vector3d &ray_direction, Vector3d &p, Vector3d &N, std::pair<int, double> &nearest, const int cur_node) {
	const double box_t = ray_box_intersection(ray_origin, ray_direction, nodes[cur_node].bbox);
	// optimization: terminate if distance to the bbox is farther than the
	// known nearest object
	if (box_t < 0 || box_t >= nearest.second) {
		return;
	}
	if (nodes[cur_node].triangle != -1) {  // leaf
		update_nearest_triangle(ray_origin, ray_direction, p, N, nearest, nodes[cur_node].triangle);
	}
	else {
		find_nearest_object(ray_origin, ray_direction, p, N, nearest, nodes[cur_node].left);
		find_nearest_object(ray_origin, ray_direction, p, N, nearest, nodes[cur_node].right);
	}
}

////////////////////////////////////////////////////////////////////////////////
// Raytracer code
////////////////////////////////////////////////////////////////////////////////

Color shoot_ray(const Vector3d &ray_origin, const Vector3d &ray_direction) {
    //Intersection point and normal, these are output of find_nearest_object
    Vector3d p, N;

    const int nearest_object = find_nearest_object(ray_origin, ray_direction, p, N).first;

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

        Color diff_color = obj_diffuse_color;

        // Diffuse contribution
        const Vector3d Li = (light.position - p).normalized();
        const Color diffuse = diff_color * std::max(Li.dot(N), 0.0);

        // Specular contribution
        const Vector3d Hi = (Li - ray_direction).normalized();
        const Color specular = obj_specular_color * std::pow(std::max(N.dot(Hi), 0.0), obj_specular_exponent);
        // Vector3d specular(0, 0, 0);

        // Attenuate lights according to the squared distance to the lights
        const Vector3d D = light.position - p;
        lights_color += (diffuse + specular).cwiseProduct(light.color) / D.squaredNorm();
    }

    // Rendering equation
    Color C = ambient_color + lights_color;

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
    double image_y = tan(field_of_view / 2) * focal_length;
    double image_x = aspect_ratio * image_y;

    // The pixel grid through which we shoot rays is at a distance 'focal_length'
    const Vector3d image_origin(-image_x, image_y, camera_position[2] - focal_length);
    const Vector3d x_displacement(2.0 / w * image_x, 0, 0);
    const Vector3d y_displacement(0, -2.0 / h * image_y, 0);

    for (size_t i = 0; i < w; ++i) {
        for (size_t j = 0; j < h; ++j) {
            const Vector3d pixel_center = image_origin + (i + 0.5) * x_displacement + (j + 0.5) * y_displacement;

            // Prepare the ray
            Vector3d ray_origin;
            Vector3d ray_direction;

            if (projection_type == perspective) {  // Perspective camera
                ray_origin = camera_position;
                ray_direction = (pixel_center - camera_position).normalized();
            }
            else if (projection_type == orthographic) {  // Orthographic camera
                ray_direction = Vector3d(0, 0, -1);
                ray_origin = pixel_center - ray_direction * focal_length;
            }

            const Color C = shoot_ray(ray_origin, ray_direction);
            R(i, j) = C(0);
            G(i, j) = C(1);
            B(i, j) = C(2);
            A(i, j) = C(3);
        }
    }

    // Save to png
    write_matrix_to_png(R, G, B, A, filename);
}

////////////////////////////////////////////////////////////////////////////////

int main(int argc, char *argv[]) {
	std::ios::sync_with_stdio(false);

	static struct option long_options[] = {
		{"filename",              required_argument, 0, 0  },
		{"mesh_filename",         required_argument, 0, 6  },
		{"brute_force",           no_argument,       0, 7  },
		{"centroid_sorting_split",no_argument,       0, 8  },
		{"construction",          required_argument, 0, 9  },
		{"criterion",             required_argument, 0, 10 },
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
	size_t w = 640, h = 480;
	std::string filename("raytrace.png");

	int opt;
	while ((opt = getopt_long(argc, argv, "w:h:f:p:g:", long_options, NULL)) != -1) {
		switch (opt) {
		case 0:
			filename = optarg;
			break;
		case 6:
			mesh_filename = optarg;
			break;
		case 7:
			use_bvh = false;
			break;
		case 8:
			centroid_sorting_split = true;
			break;
		case 9:
			switch (optarg[0]) {
				case 'T':
				case 't':
					top_down_construction = true;
					break;
				case 'B':
				case 'b':
					top_down_construction = false;
					break;
				default:
					std::cerr << "Unknown construction type: " << optarg << std::endl;
					exit(EXIT_FAILURE);
			}
			break;
		case 10:
			if (optarg == std::string("centroid_distance")) {
				criterion = &squared_centroid_distance;
			}
			else if (optarg == std::string("volume_increase")) {
				criterion = &volume_increase;
			}
			else {
				std::cerr << "Unknown criterion: " << optarg << std::endl;
				exit(EXIT_FAILURE);
			}
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
		filename);

    return 0;
}
