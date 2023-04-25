#include "SDLViewer.h"

#include <Eigen/Core>

#include <functional>
#include <iostream>
#include <map>
#include <getopt.h>

#include "raster.h"

// Image writing library
#define STB_IMAGE_WRITE_IMPLEMENTATION // Do not include this line twice in your project!
#include "stb_image_write.h"


class Triangle {
public:
	VertexAttributes vas[3];
	Triangle() {
	}
	Triangle(const VertexAttributes _vas[3]) {
		for (int i = 0; i < 3; i++) vas[i] = _vas[i];
	}
	void selected() {
		for (int i = 0; i < 3; i++) vas[i].color = Color(0.5, 0.5, 0.5, 1);
	}
	void unselected() {
		for (int i = 0; i < 3; i++) vas[i].color = Color(1, 1, 1, 1);
	}
	void move(const Position4& v) {
		for (int i = 0; i < 3; i++) vas[i].position += v;
	}
};


int main(int argc, char *argv[]) {
	static struct option long_options[] = {
		{"redraw_interval", required_argument, 0, 'i'},
		{"width"          , required_argument, 0, 'w'},
		{"height"         , required_argument, 0, 'h'},
		{"line_thickness" , required_argument, 0, 't'},
		{0                ,                 0, 0,   0},
	};

	int redraw_interval = 30;
    int width = 500, height = 500;
	double line_thickness = 1;

	int opt;
	while ((opt = getopt_long(argc, argv, "i:w:h:t:", long_options, NULL)) != -1) {
		switch (opt) {
		case 'i':
			redraw_interval = atoi(optarg);
			break;
		case 'w':
			width = atoi(optarg);
			break;
		case 'h':
			height = atoi(optarg);
			break;
		case 't':
			line_thickness = atof(optarg);
			break;
		default:
			std::cerr << "Read the code for usage." << std::endl;
			exit(EXIT_FAILURE);
		}
	}

    // The Framebuffer storing the image rendered by the rasterizer
	Eigen::Matrix<FrameBufferAttributes,Eigen::Dynamic,Eigen::Dynamic> frameBuffer(width, height);

	// Global Constants (empty in this example)
	UniformAttributes uniform;

	// Basic rasterization program
	Program program;

	// The vertex shader is the identity
	program.VertexShader = [](const VertexAttributes& va, const UniformAttributes& uniform) {
		return va;
	};

	// The fragment shader uses a fixed color
	program.FragmentShader = [](const VertexAttributes& va, const UniformAttributes& uniform) {
		return FragmentAttributes(va.color, va.obj_id);
	};

	// The blending shader converts colors between 0 and 1 to uint8
	program.BlendingShader = [](const FragmentAttributes& fa, const FrameBufferAttributes& previous) {
		return FrameBufferAttributes(color_to_color8(fa.color), fa.obj_id);
	};

	// triangles
	std::map<int, Triangle> triangles;
	int new_obj_id = 1, selected_obj_id = 0;
	bool moving_selected_obj = false;

	// new triangle
	int n_new_vertices = 0;
	VertexAttributes new_vertices[3];
	for (int i = 0; i < 3; i++) new_vertices[i].color = Color(1, 1, 1, 1);
	// mode
	char mode = 0;

    // Initialize the viewer and the corresponding callbacks
    SDLViewer viewer;
    viewer.init("Viewer Example", width, height);

	auto sdl_position_to_position4 = [&](const int x, const int y) {
		return Position4(
			(double(x)/double(width) * 2) - 1,
			(double(height-1-y)/double(height) * 2) - 1,
			0,
			1
		);
	};
	auto sdl_vector_to_position4 = [&](const int x, const int y) {
		return Position4(
			x/double(width) * 2,
			-y/double(height) * 2,
			0,
			0
		);
	};

    viewer.mouse_move = [&](int x, int y, int xrel, int yrel) {
		//std::cerr << "mouse_move(x=" << x << ", y=" << y << ", xrel=" << xrel << ", yrel=" << yrel << ")" << std::endl;

		Position4 position = sdl_position_to_position4(x, y);

		switch (mode) {
		case 'i':
			if (n_new_vertices) {
				new_vertices[n_new_vertices].position = position;
				viewer.redraw_next = true;
			}
			break;
		case 'o':
			if (selected_obj_id && moving_selected_obj) {
				auto v = sdl_vector_to_position4(xrel, yrel);
				//std::cerr << "move triangle " << selected_obj_id << " by " << v.transpose() << std::endl;
				triangles[selected_obj_id].move(v);
				viewer.redraw_next = true;
			}
			break;
		case 'p':
			break;
		default:
			break;
		}
    };

    viewer.mouse_pressed = [&](int x, int y, bool is_pressed, int button, int clicks) {
		std::cerr << "mouse_pressed(x=" << x << ", y=" << y << ", is_pressed=" << is_pressed << ")" << std::endl;

		if (mode != 'o') {
			if (selected_obj_id) {
				triangles[selected_obj_id].unselected();
				selected_obj_id = 0;
				viewer.redraw_next = true;
			}
		}

		switch (mode) {
		case 'i':
			if (!is_pressed) {
				Position4 position = sdl_position_to_position4(x, y);
				new_vertices[n_new_vertices].position = position;
				n_new_vertices++;
				if (n_new_vertices == 3) {
					for (int i = 0; i < 3; i++) new_vertices[i].obj_id = new_obj_id;
					triangles[new_obj_id++] = Triangle(new_vertices);
					n_new_vertices = 0;
				}
				else {
					new_vertices[n_new_vertices].position = position;
				}
				viewer.redraw_next = true;
			}
			break;
		case 'o':
			if (is_pressed) {
				if (selected_obj_id) {
					triangles[selected_obj_id].unselected();
					viewer.redraw_next = true;
				}
				selected_obj_id = frameBuffer(x, height-1-y).obj_id;
				std::cerr << "selected_obj_id = " << selected_obj_id << std::endl;
				if (selected_obj_id) {
					triangles[selected_obj_id].selected();
					moving_selected_obj = true;
				}
			}
			else {
				moving_selected_obj = false;
			}
			break;
		case 'p':
			if (is_pressed) {
				selected_obj_id = frameBuffer(x, height-1-y).obj_id;
				std::cerr << "selected_obj_id = " << selected_obj_id << std::endl;
				if (selected_obj_id) {
					triangles.erase(selected_obj_id);
					selected_obj_id = 0;
				}
				viewer.redraw_next = true;
			}
			break;
		default:
			break;
		}

		viewer.update();
    };

    viewer.mouse_wheel = [&](int dx, int dy, bool is_direction_normal) {
    };

    viewer.key_pressed = [&](char key, bool is_pressed, int modifier, int repeat) {
		std::cerr << "key_pressed: " << "key = " << key << " is_pressed = " << is_pressed << std::endl;

		switch (key) {
		case 'i':
		case 'o':
		case 'p':
			mode = key;
			break;
		default:
			break;
		}
    };

    viewer.redraw = [&](SDLViewer &viewer) {
        // Clear the framebuffer
        for (size_t i = 0; i < frameBuffer.rows(); i++)
            for (size_t j = 0; j < frameBuffer.cols(); j++)
                frameBuffer(i,j) = FrameBufferAttributes();

		std::cerr << "n_triangles = " << triangles.size() << std::endl;
		std::vector<VertexAttributes> triangle_vertices;
		for (const auto& id_obj: triangles) {
			const Triangle& triangle = id_obj.second;
			for (int i = 0; i < 3; i++) triangle_vertices.push_back(triangle.vas[i]);
		}
		for (auto& va: triangle_vertices)
			std::cerr << va.position.transpose() << std::endl;
		rasterize_triangles(program, uniform, triangle_vertices, frameBuffer);

		std::vector<VertexAttributes> line_vertices;
		// draw a line from every vertex to every other vertex
		for (int i = 0; i <= n_new_vertices; i++)
			for (int j = i + 1; j <= n_new_vertices; j++) {
				line_vertices.push_back(new_vertices[i]);
				line_vertices.push_back(new_vertices[j]);
			}
		std::cerr << "n_lines = " << line_vertices.size() / 2 << std::endl;
		for (auto& va: line_vertices)
			std::cerr << va.position.transpose() << std::endl;
		rasterize_lines(program, uniform, line_vertices, line_thickness, frameBuffer);

        // Buffer for exchanging data between rasterizer and sdl viewer
        Eigen::Matrix<uint8_t, Eigen::Dynamic, Eigen::Dynamic> R(width, height);
        Eigen::Matrix<uint8_t, Eigen::Dynamic, Eigen::Dynamic> G(width, height);
        Eigen::Matrix<uint8_t, Eigen::Dynamic, Eigen::Dynamic> B(width, height);
        Eigen::Matrix<uint8_t, Eigen::Dynamic, Eigen::Dynamic> A(width, height);

        for (size_t i = 0; i < frameBuffer.rows(); i++) {
            for (size_t j = 0; j < frameBuffer.cols(); j++) {
                R(i,frameBuffer.cols()-1-j) = frameBuffer(i,j).color(0);
                G(i,frameBuffer.cols()-1-j) = frameBuffer(i,j).color(1);
                B(i,frameBuffer.cols()-1-j) = frameBuffer(i,j).color(2);
                A(i,frameBuffer.cols()-1-j) = frameBuffer(i,j).color(3);
            }
        }
        viewer.draw_image(R, G, B, A);
    };

    viewer.launch(redraw_interval=redraw_interval);

    return 0;
}
