
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>

struct Point {
  double x;
  double y;

  Point() : x(0), y(0){};
  Point(const double xx, const double yy) : x(xx), y(yy) {}

  double cross(const Point p) const { return x * p.y - y * p.x; }
  Point operator-(const Point p) const { return Point(x - p.x, y - p.y); }
};

bool IsPointIntriangle(const Point& t1, const Point& t2, const Point& t3,
                       const Point& p) {
  constexpr double kEpsilon = 1e-9;
  const Point v1p = p - t1;
  const Point v12 = t2 - t1;
  const Point v13 = t3 - t1;

  const double tmp = v12.cross(v13);
  const double a = v1p.cross(v13) / tmp;
  const double b = v12.cross(v1p) / tmp;

  if (std::abs(a) < kEpsilon || std::abs(b) < kEpsilon) return true;
  if (a < -kEpsilon || b < -kEpsilon) return false;
  return (a + b) <= 1.0 + kEpsilon;
}

int main(int argc, char* argv[]) {
  std::cout << "Hello, world" << std::endl;

  for (int i = 1; i <= 100; i++) {
    std::cout << "\ni = " << i << "\n";

    Point t1, t2, t3, p;

    std::string filename1 = "test_data/" + std::to_string(i) + ".in";
    std::ifstream in_file(filename1, std::ios::in);
    in_file >> t1.x >> t1.y;
    in_file >> t2.x >> t2.y;
    in_file >> t3.x >> t3.y;
    in_file >> p.x >> p.y;
    in_file.close();

    std::string filename2 = "test_data/" + std::to_string(i) + ".out";
    std::ifstream out_file(filename2, std::ios::in);
    bool answer;
    out_file >> answer;
    in_file.close();

    const bool result = IsPointIntriangle(t1, t2, t3, p);
    std::cout << "My result: " << result << "  The answer: " << answer
              << std::endl;

    if (answer != result) {
      std::cout << "The result is Wrong!" << std::endl;
      break;
    }
  }

  return 0;
}
