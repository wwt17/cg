#include "SDLViewer.h"

#include <Eigen/Core>

#include <functional>
#include <iostream>
#include <getopt.h>

#include "raster.h"

// Image writing library
#define STB_IMAGE_WRITE_IMPLEMENTATION // Do not include this line twice in your project!
#include "stb_image_write.h"

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
		return FragmentAttributes(va.color(0),va.color(1),va.color(2));
	};

	// The blending shader converts colors between 0 and 1 to uint8
	program.BlendingShader = [](const FragmentAttributes& fa, const FrameBufferAttributes& previous) {
		return FrameBufferAttributes(fa.color[0]*255,fa.color[1]*255,fa.color[2]*255,fa.color[3]*255);
	};

	// triangle vertices
	std::vector<VertexAttributes> triangle_vertices;
	// new triangle
	int n_new_vertices = 0;
	VertexAttributes new_vertices[3];
	for (int i = 0; i < 3; i++) new_vertices[i].color = Color(1, 1, 1, 1);
	// mode
	char mode = 0;

    // Initialize the viewer and the corresponding callbacks
    SDLViewer viewer;
    viewer.init("Viewer Example", width, height);

	auto sdl_xy_to_position4 = [&](const int x, const int y) {
		return Position4(
			(double(x)/double(width) * 2) - 1,
			(double(height-1-y)/double(height) * 2) - 1,
			0,
			1
		);
	};

    viewer.mouse_move = [&](int x, int y, int xrel, int yrel) {
		//std::cerr << "mouse_move: " << std::endl;

		Position4 position = sdl_xy_to_position4(x, y);

		switch (mode) {
		case 'i':
			if (n_new_vertices) {
				new_vertices[n_new_vertices].position = position;
				viewer.redraw_next = true;
			}
			break;
		case 'o':
			break;
		case 'p':
			break;
		default:
			break;
		}

    };

    viewer.mouse_pressed = [&](int x, int y, bool is_pressed, int button, int clicks) {
		std::cerr << "mouse_pressed: " << "is_pressed = " << is_pressed << std::endl;
		if (is_pressed)
			return;

		Position4 position = sdl_xy_to_position4(x, y);

		switch (mode) {
		case 'i':
			new_vertices[n_new_vertices].position = position;
			n_new_vertices++;
			if (n_new_vertices == 3) {
				for (int i = 0; i < n_new_vertices; i++)
					triangle_vertices.push_back(new_vertices[i]);
				n_new_vertices = 0;
			}
			else {
				new_vertices[n_new_vertices].position = position;
			}
			viewer.redraw_next = true;
			break;
		case 'o':
			break;
		case 'p':
			break;
		default:
			break;
		}
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
                frameBuffer(i,j).color << 0,0,0,1;

		std::cerr << "n_triangle_vertices = " << triangle_vertices.size() << std::endl;
		for (auto& va: triangle_vertices)
			std::cerr << va.color.transpose() << std::endl;
		rasterize_triangles(program, uniform, triangle_vertices, frameBuffer);

		std::vector<VertexAttributes> line_vertices;
		// draw a line from every vertex to every other vertex
		for (int i = 0; i <= n_new_vertices; i++)
			for (int j = i + 1; j <= n_new_vertices; j++) {
				line_vertices.push_back(new_vertices[i]);
				line_vertices.push_back(new_vertices[j]);
			}
		std::cerr << "n_line_vertices = " << line_vertices.size() << std::endl;
		for (auto& va: line_vertices)
			std::cerr << va.color.transpose() << std::endl;
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
