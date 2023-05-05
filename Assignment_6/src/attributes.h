#pragma once

#include <cmath>
#include <map>
#include <vector>

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


// interpolation
enum InterpolationType {linear, bezier};


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
	std::vector<Transform4> transforms; // length of n_keyframes
	Triangle() {
	}
	Triangle(const VertexAttributes _vas[3], const int n_keyframes=1): transforms(std::vector<Transform4>(n_keyframes, Transform4::Identity())) {
		for (int i = 0; i < 3; i++) vas[i] = _vas[i];
	}
	void selected() {
		for (int i = 0; i < 3; i++) vas[i].color[3] = 0.5;
	}
	void unselected() {
		for (int i = 0; i < 3; i++) vas[i].color[3] = 1;
	}
	void insert_keyframe(const int keyframe_id) {  // insert a new keyframe after this keyframe
		transforms.insert(transforms.begin() + (keyframe_id + 1), transforms[keyframe_id]);
	}
	void remove_keyframe(const int keyframe_id) {
		transforms.erase(transforms.begin() + keyframe_id);
	}
	void clear_keyframes() {
		transforms.resize(1);
	}
	void append_transform(const Transform4& new_transform, const int keyframe_id) {
		Transform4 &transform = transforms[keyframe_id];
		transform = new_transform * transform;
	}
	Position4 barycenter(const int keyframe_id) const {
		Position4 S(0, 0, 0, 0);
		for (int i = 0; i < 3; i++) S += normalize(vas[i].position);
		return transforms[keyframe_id] * (S / 3);
	}
	Transform4 transform_at_time(const double t, const InterpolationType interpolation) const {
		if (interpolation == linear) {
			const double t_keyframes = t * (int(transforms.size()) - 1);
			const int keyframe_id = int(t_keyframes + 1e-7);
			const double t_inter_keyframe = t_keyframes - keyframe_id;
			return (1 - t_inter_keyframe) * transforms[keyframe_id] + t_inter_keyframe * transforms[keyframe_id + 1];
		}
		else { // Bezier
			const int n = transforms.size() - 1;
			// c[i] = (n choose i)
			double c[n+1];
			c[0] = 1;
			for (int i = 1; i * 2 <= n; i++) c[i] = c[i-1] * (n-i+1)/i;
			for (int i = 0; i * 2 < n; i++) c[n-i] = c[i];
			// p0[i] = t^i
			double p0[n+1];
			p0[0] = 1;
			for (int i = 1; i <= n; i++) p0[i] = p0[i-1] * t;
			// p1[i] = (1-t)^i
			double p1[n+1];
			p1[0] = 1;
			for (int i = 1; i <= n; i++) p1[i] = p1[i-1] * (1-t);
			// B(n,i) = (n choose i) * t^i * (1-t)^(n-i) = c[i] * p0[i] * p1[n-i]
			// ans = sum_i B(n,i) * a(i)
			Transform4 ans; ans.setZero();
			for (int i = 0; i <= n; i++) ans += (c[i] * p0[i] * p1[n-i]) * transforms[i];
			return ans;
		}
	}
};

class UniformAttributes {
public:
	int n_keyframes;
	int n_inter_frames;
	int cur_frame_id;  // in [0, (n_keyframes - 1) * n_inter_frames]
	InterpolationType interpolation;
	std::map<int, Triangle> triangles;
	Transform4 view_transform;
	UniformAttributes(const int n_keyframes=1, const int n_inter_frames=1, const int cur_frame_id=0, const InterpolationType interpolation=linear): n_keyframes(n_keyframes), n_inter_frames(n_inter_frames), cur_frame_id(cur_frame_id), interpolation(interpolation) {
		view_transform.setIdentity();
	}
	int cur_keyframe_id() const {
		return cur_frame_id / n_inter_frames;
	}
	double cur_t() const {
		const int n_frames = (n_keyframes - 1) * n_inter_frames;
		return n_frames ? cur_frame_id / double(n_frames) : 0.;
	}
};
