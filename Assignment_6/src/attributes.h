#pragma once

#include <cmath>

#include <Eigen/Core>


const double pi = acos(-1);


typedef Eigen::Vector4d Position4;
typedef Eigen::Matrix4d Transform4;
typedef Eigen::Vector4d Color;
typedef Eigen::Matrix<uint8_t,4,1> Color8;


inline Position4 normalize(const Position4& position) {
	return position / position[3];
}


inline Transform4 movement(const Position4& v) {
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
	return movement(v) * rotation_xy(theta) * movement(-v);
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
	return movement(v) * scaling(s) * movement(-v);
}


inline Color8 color_to_color8(const Color& color) {
	return (255*color).cast<uint8_t>();
}


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

class UniformAttributes {
public:
};
