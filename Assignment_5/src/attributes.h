#pragma once

#include <vector>
#include <cmath>
#include <Eigen/Core>

typedef Eigen::Vector4d            Color ;  // RGBA Color in [0, 1]
typedef Eigen::Matrix<uint8_t,4,1> Color8;  // RGBA Color in [0, 255]
typedef Eigen::Vector4d            Position;  // Homogeneous coordinates

const double eps = 1e-7, inf = 1/0., pi = acos(-1);

inline Position position3_to_position(const Eigen::Vector3d& position3) {
	return Position(position3[0], position3[1], position3[2], 1);
}

inline Color no_transparency(const Color& color) {
	return Color(color[0], color[1], color[2], 1);
}

inline Color blend(const Color& back_color, const Color& front_color) {
	const Color obscured_back_color = (1 - front_color[3]) * back_color[3] * no_transparency(back_color),
	            obscured_front_color = front_color[3] * no_transparency(front_color);
	Color new_color = obscured_front_color + obscured_back_color;
	new_color.topRows<3>() /= new_color(3);
	return new_color;
}

class VertexAttributes {
	public:
	Position position;
	Color color;
	Eigen::Vector3d normal;

	VertexAttributes(
		const Position& position=Position(0, 0, 0, 1),
		const Color& color=Color(0, 0, 0, 1),
		const Eigen::Vector3d& normal=Eigen::Vector3d(0, 0, 0)
	): position(position), color(color), normal(normal) {
	}

    // Interpolates the vertex attributes
    static VertexAttributes interpolate(
        const VertexAttributes& a,
        const VertexAttributes& b,
        const VertexAttributes& c,
        const double alpha, 
        const double beta, 
        const double gamma
    ) {
		return VertexAttributes(
			alpha * a.position + beta * b.position + gamma * c.position,
			alpha * a.color    + beta * b.color    + gamma * c.color,
			alpha * a.normal   + beta * b.normal   + gamma * c.normal
		);
    }
};

class FragmentAttributes {
	public:
	Color color;
	Position position;

	FragmentAttributes(
		const Color& color=Color(0, 0, 0, 1),
		const Position& position=Position(0, 0, 0, 1)
	): color(color), position(position) {
	}
};

class FrameBufferAttributes {
	public:

	class Layer {
	public:
		Color color;
		double depth;

		Layer(const Color& color, const double depth): color(color), depth(depth) {
		}
	};

	std::vector<Layer> layers;
	Color8 color;

	FrameBufferAttributes(
		const Color& bg_color=Color(0, 0, 0, 0)
	): color((bg_color * 255).cast<uint8_t>()) {
		layers.emplace_back(bg_color, inf);
	}

	void update_color() {
		auto it = layers.begin();
		Color cur_color = it->color;
		for (++it; it != layers.end(); ++it) {
			cur_color = blend(cur_color, it->color);
		}
		color = (cur_color * 255).cast<uint8_t>();
	}

	void blend_with(
		const Color& blended_color,
		const double blended_depth
	) {
		auto it = layers.begin();
		while (it != layers.end() && it->depth > blended_depth - 5e-3) ++it;
		layers.insert(it, Layer(blended_color, blended_depth));
		update_color();
	}
};

class UniformAttributes {
	public:
	Color color;
};
