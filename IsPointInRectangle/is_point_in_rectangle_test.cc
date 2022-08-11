
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <vector>

constexpr double kEpsilon = 1e-9;
const double ww[4][2] = {{0.5, 0.5}, {0.5, -0.5}, {-0.5, -0.5}, {-0.5, 0.5}};

struct Vec2d {
  double x;
  double y;

  Vec2d() = default;
  Vec2d(const double xx, const double yy) : x(xx), y(yy) {}

  double length() const { return std::hypot(x, y); }
  double dot(const Vec2d p) const { return x * p.x + y * p.y; }
  double cross(const Vec2d p) const { return x * p.y - y * p.x; }
  Vec2d operator-(const Vec2d p) const { return Vec2d(x - p.x, y - p.y); }
  Vec2d operator+(const Vec2d p) const { return Vec2d(x + p.x, y + p.y); }
};

class Rect {
 public:
  Rect() = default;
  Rect(const double x, const double y, const double l, const double w,
       const double yaw)
      : x_(x), y_(y), l_(l), w_(w), yaw_(yaw) {
    std::cout << "center: " << x_ << " " << y_ << " l " << l_ << " w_ " << w_
              << " yaw " << yaw_ << std::endl;
  }

  std::vector<Vec2d> GetCorners() const {
    std::vector<Vec2d> corners;
    Vec2d corner1, corner2;
    for (int i = 0; i < 4; ++i) {
      corner1 = Vec2d(l_ * ww[i][0], w_ * ww[i][1]);
      corner2.x = corner1.x * std::cos(yaw_) - corner1.y * std::sin(yaw_);
      corner2.y = corner1.x * std::sin(yaw_) + corner1.y * std::cos(yaw_);
      corner2 = corner2 + Vec2d(x_, y_);
      corners.push_back(corner2);
    }
    return corners;
  }

  Vec2d RotateToRec(const Vec2d p) const {
    const Vec2d pp = p - Vec2d(x_, y_);
    return Vec2d(pp.x * std::cos(yaw_) + pp.y * std::sin(yaw_),
                 -pp.x * std::sin(yaw_) + pp.y * std::cos(yaw_));
  }

  bool IsPointInRect(const Vec2d p) const {
    const Vec2d pp = RotateToRec(p);
    // std::cout << "Point: pp " << pp.x << " " << pp.y << std::endl;
    return std::abs(pp.x) <= l_ / 2 && std::abs(pp.y) <= w_ / 2;
  }

  bool IsCircleOverlapRect(const Vec2d center, const double r) const {
    if (r < 0) return false;
    const Vec2d center_to_circle = RotateToRec(center);
    const Vec2d center_to_circle_abs(std::abs(center_to_circle.x),
                                     std::abs(center_to_circle.y));
    const Vec2d center_to_corner(l_ / 2, w_ / 2);
    const Vec2d circle_to_rect(
        std::max(center_to_circle_abs.x - center_to_corner.x, 0.0),
        std::max(center_to_circle_abs.y - center_to_corner.y, 0.0));
    return circle_to_rect.length() <= r;
  }

  bool IsCircleOverlapRect2(const Vec2d center, const double r) const {
    if (r < 0) return false;
    const Vec2d center_to_corner(l_ / 2, w_ / 2);
    const double half_diagonal_length = center_to_corner.length();
    const Vec2d center_to_circle(center.x - x_, center.y - y_);
    const double c2c_length = center_to_circle.length();

    if (c2c_length > half_diagonal_length + r) return false;
    if (c2c_length <= std::min(l_, w_) / 2) return true;

    const Vec2d center_in_rect = RotateToRec(center);
    if (std::abs(center_in_rect.x) <= l_ / 2) {
      return std::abs(center_in_rect.y) <= w_ / 2 + r;
    } else if (std::abs(center_in_rect.y) <= w_ / 2) {
      return std::abs(center_in_rect.x) <= l_ / 2 + r;
    } else {
      const Vec2d circle_to_corner(std::abs(center_in_rect.x) - l_ / 2,
                                   std::abs(center_in_rect.y) - w_ / 2);
      return circle_to_corner.length() <= r;
    }

    return false;
  }

 private:
  double x_;
  double y_;
  double l_;
  double w_;
  double yaw_;
};

double RandomDouble(const int scale) {
  const double r =
      (static_cast<double>(std::rand()) / static_cast<double>(RAND_MAX) - 0.5) *
      2;
  return r * scale;
}

int main(int argc, char* argv[]) {
  std::cout << "Hello, world" << std::endl;

  // const Rect rect(10.0, 10.0, 2.0, 1.0, M_PI / 4);
  // const Vec2d p(11.0, 10);
  // std::cout << "Point: " << p.x << " " << p.y << std::endl;
  // std::cout << "IsPointInRect: " << rect.IsPointInRect(p) << std::endl;

  // const Rect rect(0.0, 0.0, 2.0, 1.0, 0.0);
  // const Vec2d center(0.0, 1.5);
  // const double radius = 1.0;
  // const bool result = rect.IsCircleOverlapRect(center, radius);
  // std::cout << "Circle center: " << center.x << " " << center.y
  //           << " radius: " << radius << std::endl;
  // std::cout << "IsCircleOverlapRect: " << result << std::endl;

  for (int i = 0; i < 1e7; ++i) {
    std::cout << "\ni = " << i << std::endl;
    const Rect rect(RandomDouble(1000) + 500, RandomDouble(1000) + 500,
                    std::abs(RandomDouble(200) * 2),
                    std::abs(RandomDouble(200) * 2), RandomDouble(M_PI / 2));

    const Vec2d center(RandomDouble(1000) + 500, RandomDouble(1000) + 500);
    const double radius = RandomDouble(1000) + 1000;
    const bool result1 = rect.IsCircleOverlapRect(center, radius);
    const bool result2 = rect.IsCircleOverlapRect2(center, radius);
    std::cout << "Circle center: " << center.x << " " << center.y
              << " radius: " << radius << std::endl;
    std::cout << "result1: " << result1 << " result2: " << result2 << std::endl;

    if (result1 != result2) {
      std::cout << "\nSomething is wrong!" << std::endl;
      std::vector<Vec2d> corners = rect.GetCorners();
      std::cout << "Rect corners: " << corners[0].x << ", " << corners[0].y
                << " " << corners[1].x << ", " << corners[1].y << " "
                << corners[2].x << ", " << corners[2].y << " " << corners[3].x
                << ", " << corners[3].y << std::endl;
      std::cout
          << "<circle cx=\"" << center.x << "\" cy=\"" << center.y << "\" r=\""
          << radius
          << "\" stroke=\" black \" fill=\"transparent\" stroke-width=\"2\"/>"
          << std::endl;
      break;
    }
  }

  return 0;
}
