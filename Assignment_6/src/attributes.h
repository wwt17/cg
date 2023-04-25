#pragma once

#include <cmath>
#include <map>

#include <Eigen/Core>


const double inf = 1/0., pi = acos(-1);


typedef Eigen::Vector4d Position4;
typedef Eigen::Matrix4d Transform4;
typedef Eigen::Vector4d Color;
typedef Eigen::Matrix<uint8_t,4,1> Color8;


// position utils

inline Position4 normalize(const Position4& position) {
	return position / position[3];
}

inline double sqr_distance(const Position4& x, const Position4& y) {
	return (normalize(y) - normalize(x)).squaredNorm();
}

inline Transform4 translation(const Position4& v) {
	Transform4 transform; transform.setIdentity();
	transform.block<4,1>(0,3) += v;
	return transform;
}

inline Transform4 rotation_xy(const double theta) {
	Transform4 transform; transform.setZero();
	transform(0,0) = +cos(theta);
	transform(0,1) = -sin(theta);
	transform(1,0) = +sin(theta);
	transform(1,1) = +cos(theta);
	transform(2,2) = 1;
	transform(3,3) = 1;
	return transform;
}

inline Transform4 rotation_xy(const double theta, const Position4& center) {
	Position4 v = normalize(center); v[3] = 0;
	return translation(v) * rotation_xy(theta) * translation(-v);
}

inline Transform4 scaling(const double s) {
	Transform4 transform; transform.setZero();
	transform(0,0) = s;
	transform(1,1) = s;
	transform(2,2) = s;
	transform(3,3) = 1;
	return transform;
}

inline Transform4 scaling(const double s, const Position4& center) {
	Position4 v = normalize(center); v[3] = 0;
	return translation(v) * scaling(s) * translation(-v);
}


// color utils

inline Color8 color_to_color8(const Color& color) {
	return (255*color).cast<uint8_t>();
}

inline Color color8_to_color(const Color8& color) {
	return color.cast<double>() / 255;
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


const Color colors[11] = {
	color8_to_color(Color8(255, 255, 255, 255)),
	color8_to_color(Color8(158,   1,  66, 255)),
	color8_to_color(Color8(213,  62,  79, 255)),
	color8_to_color(Color8(244, 109,  67, 255)),
	color8_to_color(Color8(253, 174,  97, 255)),
	color8_to_color(Color8(254, 224, 139, 255)),
	color8_to_color(Color8(230, 245, 152, 255)),
	color8_to_color(Color8(171, 221, 164, 255)),
	color8_to_color(Color8(102, 194, 165, 255)),
	color8_to_color(Color8( 50, 136, 189, 255)),
	color8_to_color(Color8( 94,  79, 162, 255))
};


class VertexAttributes {
public:
	Position4 position;
	Color color;
	int obj_id;

	VertexAttributes(
		const Position4& position=Position4(0,0,0,1),
		const Color& color=Color(1,1,1,1),
		const int obj_id=0
	): position(position), color(color), obj_id(obj_id) {
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
			alpha*a.position + beta*b.position + gamma*c.position,
			alpha*a.color + beta*b.color + gamma*c.color,
			a.obj_id
		);
    }
};


class FragmentAttributes {
public:
	Color color;
	int obj_id;

	FragmentAttributes(
		const Color& color=Color(0,0,0,1),
		const int obj_id=0
	): color(color), obj_id(obj_id) {
	}
};

class FrameBufferAttributes {
public:
	Color8 color;
	int obj_id;

	FrameBufferAttributes(
		const Color8& color=Color8(0,0,0,255),
		const int obj_id=0
	): color(color), obj_id(obj_id) {
	}
};

class Triangle {
public:
	VertexAttributes vas[3];
	Transform4 transform;
	Triangle() {
	}
	Triangle(const VertexAttributes _vas[3]) {
		for (int i = 0; i < 3; i++) vas[i] = _vas[i];
		transform.setIdentity();
	}
	void selected() {
		for (int i = 0; i < 3; i++) vas[i].color[3] = 0.5;
	}
	void unselected() {
		for (int i = 0; i < 3; i++) vas[i].color[3] = 1;
	}
	void append_transform(const Transform4& new_transform) {
		transform = new_transform * transform;
	}
	Position4 barycenter() const {
		Position4 S(0, 0, 0, 0);
		for (int i = 0; i < 3; i++) S += normalize(vas[i].position);
		return transform * (S / 3);
	}
};

class UniformAttributes {
public:
	std::map<int, Triangle> triangles;
};
