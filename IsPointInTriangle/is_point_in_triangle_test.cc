
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>

constexpr double kEpsilon = 1e-9;

struct Point {
  double x;
  double y;

  Point() : x(0), y(0){};
  Point(const double xx, const double yy) : x(xx), y(yy) {}

  double length() const { return std::hypot(x, y); }
  double dot(const Point p) const { return x * p.x + y * p.y; }
  double cross(const Point p) const { return x * p.y - y * p.x; }
  Point operator-(const Point p) const { return Point(x - p.x, y - p.y); }
};

double RandomDouble() {
  const double r =
      static_cast<double>(std::rand()) / static_cast<double>(RAND_MAX) - 0.5;
  return r * 1e6;
}

void GetRandomPoints(Point& t1, Point& t2, Point& t3, Point& p, int k) {
  bool is_triangle = false;
  while (!is_triangle) {
    t1 = Point(RandomDouble(), RandomDouble());
    t2 = Point(RandomDouble(), RandomDouble());
    t3 = Point(RandomDouble(), RandomDouble());
    const Point vec1 = t2 - t1;
    const Point vec2 = t3 - t1;
    if (vec1.length() < 1e-3 || vec2.length() < 1e-3) {
      is_triangle = false;
    } else {
      is_triangle =
          std::abs(vec1.cross(vec2)) / (vec1.length() * vec2.length()) > 1e-3;
    }
    // std::cout << "cross: " << std::fixed << std::setprecision(3)
    //           << std::abs(vec1.cross(vec2)) << "\n";
  }
  // random points
  p = Point(RandomDouble(), RandomDouble());

  // p is inside triangle
  // double w1 = static_cast<double>(k - 95) / 10.0;
  // double w2 = 0.4 * (1 - w1);
  // double w3 = 1 - w1 - w2;
  // p = Point(w1 * t1.x + w2 * t2.x + w3 * t3.x,
  //           w1 * t1.y + w2 * t2.y + w3 * t3.y);
  // std::cout << "w1 " << w1 << " w2 " << w2 << " w3 " << w3 << "\n";

  // p is on a edge
  // p = Point((t1.x + t3.x) / 2, (t1.y + t3.y) / 2);

  // p is nearly on a edge
  // double w = 0.5;
  // p = Point(w * t1.x + (1 - w) * t2.x, w * t1.y + (1 - w) * t2.y);
  // w = 1 - 1e-5;
  // p = Point(w * p.x + (1 - w) * t3.x, w * p.y + (1 - w) * t3.y);

  // p is a vertex
  // p = t3;
  std::cout << "t1 = [" << t1.x << ", " << t1.y << "]\n";
  std::cout << "t2 = [" << t2.x << ", " << t2.y << "]\n";
  std::cout << "t3 = [" << t3.x << ", " << t3.y << "]\n";
  std::cout << "p = [" << p.x << ", " << p.y << "]\n";
}

bool IsPointIntriangle_ByAngle(const Point& t1, const Point& t2,
                               const Point& t3, const Point& p) {
  const Point v1 = t1 - p;
  const Point v2 = t2 - p;
  const Point v3 = t3 - p;
  std::cout << "v1: " << v1.x << " " << v1.y << " v2: " << v2.x << " " << v2.y
            << " v3: " << v3.x << " " << v3.y << "\n";
  if (v1.length() == 0 || v2.length() == 0 || v3.length() == 0) {
    return true;
  }

  const double cos_angle_1 = v1.dot(v2) / (v1.length() * v2.length());
  const double cos_angle_2 = v2.dot(v3) / (v2.length() * v3.length());
  const double cos_angle_3 = v3.dot(v1) / (v3.length() * v1.length());
  std::cout << "cos angles: " << std::fixed << std::setprecision(12)
            << v1.dot(v2) / (v1.length() * v2.length()) << " "
            << v2.dot(v3) / (v2.length() * v3.length()) << " "
            << v3.dot(v1) / (v3.length() * v1.length()) << "\n";
  if (cos_angle_1 <= -1.0 || cos_angle_2 <= -1.0 || cos_angle_3 <= -1.0)
    return true;

  const double angle_1 = std::acos(cos_angle_1);
  const double angle_2 = std::acos(cos_angle_2);
  const double angle_3 = std::acos(cos_angle_3);
  std::cout << "angles: " << angle_1 << " " << angle_2 << " " << angle_3
            << "\n";

  std::cout << "angle diff: " << std::fixed << std::setprecision(12)
            << std::abs(angle_1 + angle_2 + angle_3 - M_PI * 2) << "\n";

  return std::abs(angle_1 + angle_2 + angle_3 - M_PI * 2) < kEpsilon;
}

bool IsPointIntriangle_ByDirection(const Point& t1, const Point& t2,
                                   const Point& t3, const Point& p) {
  const Point tt1 = t1;

  const Point vv12 = t2 - t1;
  const Point vv13 = t3 - t1;
  const Point tt2 = vv12.cross(vv13) > 0 ? t2 : t3;
  const Point tt3 = vv12.cross(vv13) > 0 ? t3 : t2;

  const Point v1p = p - tt1;
  const Point v2p = p - tt2;
  const Point v3p = p - tt3;
  if (v1p.length() == 0 || v2p.length() == 0 || v3p.length() == 0) {
    return true;
  }

  const Point v12 = tt2 - tt1;
  const Point v23 = tt3 - tt2;
  const Point v31 = tt1 - tt3;

  // std::cout << "v1p: " << v1p.x << " " << v1p.y << " v2p: " << v2p.x << " "
  //           << v2p.y << " v3p: " << v3p.x << " " << v3p.y << "\n";
  // std::cout << "v12: " << v12.x << " " << v12.y << " v23: " << v23.x << " "
  //           << v23.y << " v31: " << v31.x << " " << v31.y << "\n";
  // std::cout << "cross: " << v12.cross(v1p) << " " << v23.cross(v2p) << " "
  //           << v31.cross(v3p) << "\n";

  if (v12.cross(v1p) >= -kEpsilon && v23.cross(v2p) >= -kEpsilon &&
      v31.cross(v3p) >= -kEpsilon) {
    return true;
  }
  return false;
}

bool IsPointIntriangle_ByProjection(const Point& t1, const Point& t2,
                                    const Point& t3, const Point& p) {
  const Point v1p = p - t1;
  const Point v12 = t2 - t1;
  const Point v13 = t3 - t1;

  const double tmp = v12.cross(v13);
  const double a = v1p.cross(v13) / tmp;
  const double b = v12.cross(v1p) / tmp;

  std::cout << "a " << std::fixed << std::setprecision(12) << a << " b " << b
            << " a + b - 1.0 " << a + b - 1.0 << std::endl;

  if (std::abs(a) < kEpsilon || std::abs(b) < kEpsilon) return true;
  if (a < -kEpsilon || b < -kEpsilon) return false;
  return (a + b) <= 1.0 + kEpsilon;
}

int main(int argc, char* argv[]) {
  std::cout << "Hello, world" << std::endl;

  // const Point t1(1.0, 1.0);
  // const Point t2(2.0, 1.0);
  // const Point t3(2.0, 2.0);
  // const Point p(1.5, 1.0);

  for (int i = 1; i <= 100; i++) {
    std::cout << "\ni = " << i << "\n";

    Point t1, t2, t3, p;
    // GetRandomPoints(t1, t2, t3, p, i);

    // test my data and answer
    std::string filename1 = "test_data/" + std::to_string(i) + ".in";
    std::ifstream in_file(filename1, std::ios::in);
    in_file >> t1.x >> t1.y;
    in_file >> t2.x >> t2.y;
    in_file >> t3.x >> t3.y;
    in_file >> p.x >> p.y;
    in_file.close();
    std::cout << "t1 = [" << t1.x << ", " << t1.y << "]\n";
    std::cout << "t2 = [" << t2.x << ", " << t2.y << "]\n";
    std::cout << "t3 = [" << t3.x << ", " << t3.y << "]\n";
    std::cout << "p = [" << p.x << ", " << p.y << "]\n";

    std::string filename2 = "test_data/" + std::to_string(i) + ".out";
    std::ifstream out_file(filename2, std::ios::in);
    int result;
    out_file >> result;
    in_file.close();
    std::cout << "result = " << result << "\n";

    const bool result1 = IsPointIntriangle_ByAngle(t1, t2, t3, p);
    const bool result2 = IsPointIntriangle_ByDirection(t1, t2, t3, p);
    const bool result3 = IsPointIntriangle_ByProjection(t1, t2, t3, p);

    std::cout << "IsPointIntriangle_ByAngle: " << result1 << "\n";
    std::cout << "IsPointIntriangle_ByDirection: " << result2 << "\n";
    std::cout << "IsPointIntriangle_ByProjection: " << result3 << "\n";

    if ((result1 == result2) && (result1 == result3)) {
      std::cout << "My result: " << result1 << " answer: " << result
                << std::endl;
    } else {
      std::cout << "Something is wrong!" << std::endl;
      break;
    }

    // generate test data
    // if ((result1 == result2) && (result1 == result3)) {
    //   std::string filename1 = "test_data/" + std::to_string(i) + ".in";
    //   std::ofstream in_file(filename1, std::ios::out | std::ios::trunc);
    //   in_file << std::fixed << std::setprecision(6) << t1.x << " " << t1.y
    //           << "\n"
    //           << t2.x << " " << t2.y << "\n"
    //           << t3.x << " " << t3.y << "\n"
    //           << p.x << " " << p.y << "\n";
    //   in_file.close();
    //   std::string filename2 = "test_data/" + std::to_string(i) + ".out";
    //   std::ofstream out_file(filename2, std::ios::out | std::ios::trunc);
    //   out_file << result1 << "\n";
    //   out_file.close();
    // } else {
    //   std::cout << "Something is wrong!" << std::endl;
    //   break;
    // }
  }

  return 0;
}
