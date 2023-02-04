////////////////////////////////////////////////////////////////////////////////
#include <algorithm>
#include <complex>
#include <fstream>
#include <iostream>
#include <numeric>
#include <vector>
////////////////////////////////////////////////////////////////////////////////

typedef std::complex<double> Point;
typedef std::vector<Point> Polygon;

const double eps = 1e-7;

double inline det(const Point &u, const Point &v) {
	return u.real() * v.imag() - v.real() * u.imag();
}

struct Compare {
	Point p0; // Leftmost point of the poly
	bool operator ()(const Point &p1, const Point &p2) {
		Point v1 = p1 - p0, v2 = p2 - p0;
		double d = det(v1, v2);
		return (d > eps) || ((std::abs(d) <= eps) && (std::norm(v1) < std::norm(v2)));
	}
};

bool inline salientAngle(const Point &a, const Point &b, const Point &c) {
	return det(b - a, c - a) > eps;
}

////////////////////////////////////////////////////////////////////////////////

Polygon convex_hull(std::vector<Point> &points) {
	Compare order;
	// find the leftmost point (e.g., smallest x; if there are many, find the one with smallest y)
	order.p0 = points[0];
	for (const Point& p: points) {
		if ((p.real() < order.p0.real()) || ((p.real() == order.p0.real()) && p.imag() < order.p0.imag())) {
			order.p0 = p;
		}
	}
	std::sort(points.begin(), points.end(), order);
	Polygon hull;
	hull.push_back(order.p0);
	for (const Point& p: points) {
		// skip p0
		if (p == order.p0)
			continue;

		while ((int(hull.size()) >= 2) && ! salientAngle(hull[hull.size() - 2], hull[hull.size() - 1], p)) {
			hull.pop_back();
		}
		hull.push_back(p);
	}
	return hull;
}

////////////////////////////////////////////////////////////////////////////////

std::vector<Point> load_xyz(const std::string &filename) {
	std::vector<Point> points;
	std::ifstream in(filename);
	int n;
	in >> n;
	points = std::vector<Point>(n);
	for (Point& point: points) {
		double x, y, z;
		in >> x >> y >> z;
		point = Point(x, y);
	}
	return points;
}

void save_obj(const std::string &filename, Polygon &poly) {
	std::ofstream out(filename);
	if (!out.is_open()) {
		throw std::runtime_error("failed to open file " + filename);
	}
	out << std::fixed;
	for (const auto &v : poly) {
		out << "v " << v.real() << ' ' << v.imag() << " 0\n";
	}
	out << "f";
	for (size_t i = 0; i < poly.size(); ++i) {
		out << ' ' << i+1;
	}
	out << std::endl;
}

////////////////////////////////////////////////////////////////////////////////

int main(int argc, char * argv[]) {
	if (argc <= 2) {
		std::cerr << "Usage: " << argv[0] << " points.xyz output.obj" << std::endl;
	}
	std::vector<Point> points = load_xyz(argv[1]);
	Polygon hull = convex_hull(points);
	save_obj(argv[2], hull);
	return 0;
}
