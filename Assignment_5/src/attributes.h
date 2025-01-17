#pragma once

#include <vector>
#include <string>
#include <cmath>
#include <Eigen/Core>
#include <Eigen/Geometry>

typedef Eigen::Vector4d            Color ;  // RGBA Color in [0, 1]
typedef Eigen::Matrix<uint8_t,4,1> Color8;  // RGBA Color in [0, 255]
typedef Eigen::Vector4d            Position;  // Homogeneous coordinates

const double eps = 1e-7, inf = 1/0., pi = acos(-1);

////////////////////////////////////////////////////////////////////////////////
// Utility functions
////////////////////////////////////////////////////////////////////////////////

template<class T>
T sqr(const T x) {
	return x * x;
}

inline bool endsWith(const std::string& str, const std::string& suffix) {
	return str.size() >= suffix.size() && 0 == str.compare(str.size()-suffix.size(), suffix.size(), suffix);
}

////////////////////////////////////////////////////////////////////////////////
// Classes for the scene
////////////////////////////////////////////////////////////////////////////////

class Light {
public:
	Eigen::Vector3d position;
	Color color;
	Light(const Eigen::Vector3d &position, const Color &color): position(position), color(color) {
	}
};

////////////////////////////////////////////////////////////////////////////////

inline Position position3_to_position(const Eigen::Vector3d& position3) {
	return Position(position3[0], position3[1], position3[2], 1);
}

inline Eigen::Vector3d position_to_position3(const Position& position) {
	return position.head<3>() / position[3];
}

inline Position normalize_position(const Position& position) {
	return position / position[3];
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
	Eigen::Vector3d normal;
	Color color;

	VertexAttributes(
		const Position& position=Position(0, 0, 0, 1),
		const Eigen::Vector3d& normal=Eigen::Vector3d(0, 0, 0),
		const Color& color=Color(0, 0, 0, 1)
	): position(position), normal(normal), color(color) {
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
			alpha * normalize_position(a.position) + beta * normalize_position(b.position) + gamma * normalize_position(c.position),
			alpha * a.normal   + beta * b.normal   + gamma * c.normal,
			alpha * a.color    + beta * b.color    + gamma * c.color
		);
    }

	VertexAttributes transform(
		const Eigen::Matrix4d& affine,
		const Eigen::Matrix4d& rotation
	) const {
		return VertexAttributes(
			affine * (rotation * position),
			position_to_position3(rotation * position3_to_position(normal)),
			color
		);
	}

	VertexAttributes transform(
		const Eigen::Matrix4d& affine
	) const {
		return VertexAttributes(
			affine * position,
			normal,
			color
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
		double z;

		Layer(const Color& color, const double z): color(color), z(z) {
		}
	};

	std::vector<Layer> layers;
	Color8 color;

	FrameBufferAttributes(
		const Color& bg_color=Color(0, 0, 0, 0)
	): color((bg_color * 255).cast<uint8_t>()) {
		layers.emplace_back(bg_color, -inf);
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
		const double blended_z
	) {
		auto it = layers.begin();
		while (it != layers.end() && it->z < blended_z + 5e-3) ++it;
		layers.insert(it, Layer(blended_color, blended_z));
		update_color();
	}
};

class UniformAttributes {
	public:
	Color
		color,
		obj_ambient_color,
		obj_diffuse_color,
		obj_specular_color,
		ambient_light;
	double obj_specular_exponent, alpha;
	std::vector<Light> lights;
	Eigen::Matrix4d rotation, shift, view;
	UniformAttributes():
		rotation(Eigen::Matrix4d::Identity()),
		shift(Eigen::Matrix4d::Identity()),
		view(Eigen::Matrix4d::Identity())
	{
	}
};
