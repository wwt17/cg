#pragma once

#include <Eigen/Core>


typedef Eigen::Vector4d Position4;
typedef Eigen::Vector4d Color;
typedef Eigen::Matrix<uint8_t,4,1> Color8;


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
