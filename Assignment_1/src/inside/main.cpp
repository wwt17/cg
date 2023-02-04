////////////////////////////////////////////////////////////////////////////////
#include <algorithm>
#include <complex>
#include <fstream>
#include <iostream>
#include <numeric>
#include <vector>
#include <string>
#include <sstream>
#include <limits>
////////////////////////////////////////////////////////////////////////////////

typedef std::complex<double> Point;
typedef std::vector<Point> Polygon;
const double inf = std::numeric_limits<double>::infinity();
const double eps = 1e-7;

double inline det(const Point &u, const Point &v) {
	return u.real() * v.imag() - v.real() * u.imag();
}

// Return true iff [a,b] intersects [c,d], and store the intersection in ans
bool intersect_segment(const Point &a, const Point &b, const Point &c, const Point &d, Point &ans) {
	Point u = b - a, v = d - c;
	double det_ = det(u, -v);
	if (std::abs(det_) <= eps)
		return false;
	Point w = c - a;
	double lambda_u = det(w, -v) / det_, lambda_v = det(u, w) / det_;
	ans = a + lambda_u * u;
	return (lambda_u >= -eps) && (lambda_u <= 1 + eps) && (lambda_v >= -eps) && (lambda_v <= 1 + eps);
}

////////////////////////////////////////////////////////////////////////////////

bool is_inside(const Polygon &poly, const Point &query) {
	// 1. Compute bounding box and set coordinate of a point outside the polygon
	double x_min = inf, y_min = inf, x_max = -inf, y_max = -inf;
	for (const Point& p: poly) {
		x_min = std::min(x_min, p.real());
		x_max = std::max(x_max, p.real());
		y_min = std::min(y_min, p.imag());
		y_max = std::max(y_max, p.imag());
	}
	Point outside(x_min - 1, y_min - 1);
	// 2. Cast a ray from the query point to the 'outside' point, count number of intersections
	size_t n_is = 0;
	for (size_t i = 0; i < poly.size(); i++) {
		Point intersect;
		if (intersect_segment(query, outside, poly[i], poly[(i+(size_t)1) % poly.size()], intersect)) {
			n_is++;
		}
	}
	return n_is & 1;
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

Polygon load_obj(const std::string &filename) {
	std::ifstream in(filename);
	std::vector<Point> vertexes;
	Polygon polygon;
	std::string line;
	while (std::getline(in, line)) {
		std::istringstream sin(line);
		char type;
		sin >> type;
		if (type == 'v') {
			double x, y, z;
			sin >> x >> y >> z;
			vertexes.push_back(Point(x, y));
		}
		else if (type == 'f') {
			int id;
			while (sin >> id) {
				polygon.push_back(vertexes[id - 1]);
			}
		}
	}
	return polygon;
}

void save_xyz(const std::string &filename, const std::vector<Point> &points) {
	std::ofstream out(filename);
	out << points.size() << std::endl;
	for (const Point& p: points)
		out << p.real() << ' ' << p.imag() << ' ' << 0 << std::endl;
}

////////////////////////////////////////////////////////////////////////////////

int main(int argc, char * argv[]) {
	if (argc <= 3) {
		std::cerr << "Usage: " << argv[0] << " points.xyz poly.obj result.xyz" << std::endl;
	}
	std::vector<Point> points = load_xyz(argv[1]);
	Polygon poly = load_obj(argv[2]);
	std::vector<Point> result;
	for (size_t i = 0; i < points.size(); ++i) {
		if (is_inside(poly, points[i])) {
			result.push_back(points[i]);
		}
	}
	save_xyz(argv[3], result);
	return 0;
}
