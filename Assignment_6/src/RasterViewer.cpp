#include "SDLViewer.h"

#include <Eigen/Core>

#include <functional>
#include <iostream>
#include <map>
#include <stdexcept>
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
		{"keyframe_seconds",required_argument, 0, 'k'},
		{"n_inter_frames" , required_argument, 0, 'f'},
		{0                ,                 0, 0,   0},
	};

	int redraw_interval = 30;
    int width = 500, height = 500;
	double line_thickness = 1;
	double keyframe_seconds = 1;  // time between two contiguous keyframes in seconds
	int n_inter_frames = 1;       // number of interpolating frames between two contiguous keyframes

	int opt;
	while ((opt = getopt_long(argc, argv, "i:w:h:t:k:f:", long_options, NULL)) != -1) {
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
		case 'k':
			keyframe_seconds = atof(optarg);
			break;
		case 'f':
			n_inter_frames = atoi(optarg);
			break;
		default:
			std::cerr << "Read the code for usage." << std::endl;
			exit(EXIT_FAILURE);
		}
	}

    // The Framebuffer storing the image rendered by the rasterizer
	Eigen::Matrix<FrameBufferAttributes,Eigen::Dynamic,Eigen::Dynamic> frameBuffer(width, height);

	// Global Constants
	UniformAttributes uniform(1, n_inter_frames, 0);

	// Basic rasterization program
	Program program;

	// The vertex shader is the identity
	program.VertexShader = [](const VertexAttributes& va, const UniformAttributes& uniform) {
		VertexAttributes new_va = va;
		try {
			new_va.position = uniform.triangles.at(va.obj_id).transform_at_time(uniform.cur_t()) * new_va.position;
		}
		catch (const std::out_of_range& oor) {
		}
		new_va.position = uniform.view_transform * new_va.position;
		return new_va;
	};

	// The fragment shader uses a fixed color
	program.FragmentShader = [](const VertexAttributes& va, const UniformAttributes& uniform) {
		return FragmentAttributes(va.color, va.obj_id);
	};

	// The blending shader converts colors between 0 and 1 to uint8
	program.BlendingShader = [](const FragmentAttributes& fa, const FrameBufferAttributes& previous) {
		return FrameBufferAttributes(color_to_color8(blend(color8_to_color(previous.color), fa.color)), fa.obj_id);
	};

	// triangles
	std::map<int, Triangle> &triangles = uniform.triangles;
	int new_obj_id = 1, selected_obj_id = 0;
	bool moving_selected_obj = false;
	VertexAttributes *selected_vertex = NULL;

	// new triangle
	int n_new_vertices = 0;
	VertexAttributes new_vertices[3];
	for (int i = 0; i < 3; i++) new_vertices[i].obj_id = new_obj_id;
	for (int i = 0; i < 3; i++) new_vertices[i].color = Color(1, 1, 1, 1);
	// mode
	char mode = 'i';
	// interpolation type
	InterpolationType interpolation;

    // Initialize the viewer and the corresponding callbacks
    SDLViewer viewer;
    viewer.init("Viewer Example", width, height);

	auto sdl_position_to_position4 = [&](const int x, const int y) {
		return Position4(uniform.view_transform.inverse() * Position4(
			(double(x)/double(width) * 2) - 1,
			(double(height-1-y)/double(height) * 2) - 1,
			0,
			1
		));
	};
	auto sdl_vector_to_position4 = [&](const int x, const int y) {
		return Position4(sdl_position_to_position4(x, y) - sdl_position_to_position4(0, 0));
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
				triangles[selected_obj_id].append_transform(translation(v), uniform.cur_keyframe_id());
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

		switch (mode) {
		case 'i':
			if (!is_pressed) {
				Position4 position = sdl_position_to_position4(x, y);
				new_vertices[n_new_vertices].position = position;
				n_new_vertices++;
				if (n_new_vertices == 3) {
					triangles[new_obj_id++] = Triangle(new_vertices);
					for (int i = 0; i < 3; i++) new_vertices[i].obj_id = new_obj_id;
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
					viewer.redraw_next = true;
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
		case 'c':
			if (!is_pressed) {
				Position4 position = sdl_position_to_position4(x, y);
				selected_vertex = NULL;
				double nearest_sqr_distance = inf;
				// find vertex nearest to the mouse
				for (auto& id_obj: triangles) {
					Triangle& triangle = id_obj.second;
					for (int i = 0; i < 3; i++) {
						VertexAttributes& vertex = triangle.vas[i];
						double cur_sqr_distance = sqr_distance(position, triangle.transform_at_time(uniform.cur_t()) * vertex.position);
						if (cur_sqr_distance < nearest_sqr_distance) {
							nearest_sqr_distance = cur_sqr_distance;
							selected_vertex = &vertex;
						}
					}
				}
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
		if (!is_pressed)
			return;

		switch (key) {
		case 'i':
		case 'o':
		case 'p':
		case 'c':
			if (selected_obj_id && key != 'o') {
				triangles[selected_obj_id].unselected();
				selected_obj_id = 0;
				viewer.redraw_next = true;
			}
			if (selected_vertex && key != 'c') {
				selected_vertex = NULL;
			}
			mode = key;
			break;
		case '\'': // insert new keyframe after the last keyframe
			uniform.cur_frame_id = (uniform.n_keyframes - 1) * uniform.n_inter_frames; // jump to the last frame
		case ';':  // insert new keyframe after the current keyframe
			{
				uniform.n_keyframes++;
				int cur_keyframe_id = uniform.cur_keyframe_id();
				for (auto& id_obj: triangles) {
					Triangle& triangle = id_obj.second;
					triangle.insert_keyframe(cur_keyframe_id);
				}
				cur_keyframe_id++;
				uniform.cur_frame_id = cur_keyframe_id * uniform.n_inter_frames; // jump to the new frame
			}
			viewer.redraw_next = true;
			break;
		case '/':  // remove the current keyframe
			if (uniform.n_keyframes > 1) {
				int cur_keyframe_id = uniform.cur_keyframe_id();
				for (auto& id_obj: triangles) {
					Triangle& triangle = id_obj.second;
					triangle.remove_keyframe(cur_keyframe_id);
				}
				uniform.cur_frame_id = cur_keyframe_id * uniform.n_inter_frames; // jump to the next frame
				uniform.n_keyframes--;
			}
			viewer.redraw_next = true;
			break;
		case '\\': // clear all (except the first) keyframes
			{
				for (auto& id_obj: triangles) {
					Triangle& triangle = id_obj.second;
					triangle.clear_keyframes();
				}
				uniform.cur_frame_id = 0;
				uniform.n_keyframes = 1;
			}
			viewer.redraw_next = true;
			break;
		case ' ':  // play/advance the animation using linear interpolation
		case '\t': // play/advance the animation using Bezier interpolation
			interpolation = key == ' ' ? linear : bezier;
			if (++uniform.cur_frame_id > (uniform.n_keyframes - 1) * uniform.n_inter_frames)
				uniform.cur_frame_id = 0; // jump to the first frame
			viewer.redraw_next = true;
			break;
		case 'h':
		case 'j':
		case 'k':
		case 'l':
			if (selected_obj_id) {
				Triangle& triangle = triangles[selected_obj_id];
				Transform4 transform;
				switch (key) {
					case 'h': transform = rotation_xy(-10/180.*pi, triangle.barycenter(uniform.cur_keyframe_id())); break;
					case 'j': transform = rotation_xy(+10/180.*pi, triangle.barycenter(uniform.cur_keyframe_id())); break;
					case 'k': transform = scaling(1+0.25, triangle.barycenter(uniform.cur_keyframe_id())); break;
					case 'l': transform = scaling(1-0.25, triangle.barycenter(uniform.cur_keyframe_id())); break;
				}
				triangle.append_transform(transform, uniform.cur_keyframe_id());
				viewer.redraw_next = true;
			}
			break;
		case '=': // temporary solution for not recognizing '+'
		case '+': uniform.view_transform = scaling(1+0.2) * uniform.view_transform; viewer.redraw_next = true; break;
		case '-': uniform.view_transform = scaling(1-0.2) * uniform.view_transform; viewer.redraw_next = true; break;
		case 'w': uniform.view_transform = translation(Position4(0,-2*0.2,0,0)) * uniform.view_transform; viewer.redraw_next = true; break;
		case 'a': uniform.view_transform = translation(Position4(+2*0.2,0,0,0)) * uniform.view_transform; viewer.redraw_next = true; break;
		case 's': uniform.view_transform = translation(Position4(0,+2*0.2,0,0)) * uniform.view_transform; viewer.redraw_next = true; break;
		case 'd': uniform.view_transform = translation(Position4(-2*0.2,0,0,0)) * uniform.view_transform; viewer.redraw_next = true; break;
		default:
			if (mode == 'c' && selected_vertex && key >= '0' && key <= '9') {
				selected_vertex->color = colors[key-'0'];
				viewer.redraw_next = true;
			}
			break;
		}
    };

    viewer.redraw = [&](SDLViewer &viewer) {
        // Clear the framebuffer
        for (size_t i = 0; i < frameBuffer.rows(); i++)
            for (size_t j = 0; j < frameBuffer.cols(); j++)
                frameBuffer(i,j) = FrameBufferAttributes();

		//std::cerr << "n_triangles = " << triangles.size() << std::endl;
		std::vector<VertexAttributes> triangle_vertices;
		for (const auto& id_obj: triangles) {
			const Triangle& triangle = id_obj.second;
			for (int i = 0; i < 3; i++) triangle_vertices.push_back(triangle.vas[i]);
		}
		//for (auto& va: triangle_vertices) std::cerr << va.position.transpose() << std::endl;
		rasterize_triangles(program, uniform, triangle_vertices, frameBuffer);

		std::vector<VertexAttributes> line_vertices;
		// draw a line from every vertex to every other vertex
		for (int i = 0; i <= n_new_vertices; i++)
			for (int j = i + 1; j <= n_new_vertices; j++) {
				line_vertices.push_back(new_vertices[i]);
				line_vertices.push_back(new_vertices[j]);
			}
		//std::cerr << "n_lines = " << line_vertices.size() / 2 << std::endl;
		//for (auto& va: line_vertices) std::cerr << va.position.transpose() << std::endl;
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
