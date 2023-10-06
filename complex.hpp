#pragma once
#include <cmath>
#include <istream>  // << overload
#include <ostream>  // >> overload
#include <cstdlib>  // is digit

namespace nl {
  struct polar {
    double r;
    double theta;

    friend std::ostream& operator<< (std::ostream& out, polar const& p) {
      out << "(r=" << p.r << ", theta=" << p.theta << ')';
      return out;
    }
  };
  struct complex {
    double r;
    double i;

    complex() { r = i = 0.0; }
    complex(double real) { r = real; i = 0.0; }
    complex(double real, double imag) { r = real; i = imag; }
    complex(std::string str) { *this = complex::fromString(str); }
    complex(const char str[]) { *this = complex(std::string{str}); }

    friend complex operator+ (complex x, complex y) {
      return complex{x.r + y.r, x.i + y.i};
    }

    friend complex operator- (complex x, complex y) {
      return complex {x.r - y.r, x.i - y.i};
    }

    friend complex operator* (complex x, complex y) {
      return complex {
        x.r * y.r - x.i * y.i,
        x.i * y.r + x.r * y.i
      };
    }

    friend complex operator/ (complex x, complex y) {
      const double div = y.norm();
      return complex {
        (x.r * y.r + x.i * y.i) / div,
        (x.i * y.r - x.r * y.i) / div
      };
    }

    complex operator+ () {
      return complex {r, i};
    };

    complex operator- () {
      return complex {-r, -i};
    }

    complex& operator+= (complex x) {
      *this = *this + x;
      return *this;
    }

    complex& operator-= (complex x) {
      *this = *this - x;
      return *this;
    }

    complex& operator*= (complex x) {
      *this = *this * x;
      return *this;
    }

    complex& operator/= (complex x) {
      *this = *this / x;
      return *this;
    }

    friend bool operator== (complex x, complex y) {
      return x.r == y.r && x.i == y.i;
    }

    friend bool operator!= (complex x, complex y) {
      return x.r != y.r || x.i != y.i;
    }

    ////// ! IO: start

    void _formatImaginary(std::ostream& out, char prefix = '\0') const {
      if (i ==  0) return;
      if (i ==  1) { out << prefix << 'i'; return; }
      if (i == -1) { out << "-i"; return; }

      out << (i > 0 ? prefix : '\0') << i << 'i';
    }

    static double _parseDouble(std::string& str, int& index) {
      if (index >= (int) str.length()) return 0.0;
      // rule:
      // (+|-) [digit] (. [digit]) (e (+|-) [digit])
      std::string number;

      auto _parseUnary = [&str, &index, &number]() {
        char c = str[index];
        if (c == '+' || c == '-') {
          ++index;
          number += c;
        }
      };

      auto _parseDigit = [&str, &index, &number]() {
        char c = str[index];
        while (std::isdigit(c)) {
          number += c;
          c = str[++index];
        }
      };

      auto _scanChar = [&str, &index, &number](char x) -> bool {
        if (str[index] == x) {
          ++index;
          number += x;
          return true;
        }
        return false;
      };

      _parseUnary();
      _parseDigit();
      if (_scanChar('.')) {
        _parseDigit();
      }
      if (_scanChar('e')) {
        _parseUnary();
        _parseDigit();
      }
      if (number.length() == 0) return 0.0;
      if (number == std::string{'-'}) return -1.0;
      if (number == std::string{'+'}) return 1.0;
      return std::stod(number);
    }

    static complex fromString(std::string& str) {
      try {
        if (str == std::string{'i'}) return complex {0,1};
        if (str == std::string{"-i"}) return complex {0,-1};
        int index = 0;

        auto _scanChar = [&str, &index](char x) -> bool {
          if (str[index] == x) {
            ++index;
            return true;
          }
          return false;
        };

        double one = _parseDouble(str, index);
        if (index == 0) throw std::invalid_argument(NULL); // nothing can be matched!
        if (_scanChar('i')) return complex{0,one};

        _scanChar('+'); // ignore '+' because numbers are positive by default.
        double two = _parseDouble(str, index);
        return complex{one, two};
      } catch (...) {
        throw std::invalid_argument(std::string{"invalid complex number: "} +  str);
      }
    }


    friend std::ostream& operator<< (std::ostream& out, complex const& x) {
      if (x.i == 0) out << x.r;
      else if (x.r == 0) x._formatImaginary(out);
      else { out << x.r; x._formatImaginary(out, '+'); }
      return out;
    }
    friend std::istream& operator>> (std::istream& in, complex& x) {
      std::string str; in >> str;
      x = complex::fromString(str);
      return in;
    }

    //////// ! IO: end

    // check if two number is close together.
    static bool isClose(complex x, complex y, double epsilon = 1e-10) {
      return std::abs(x.r - y.r) < epsilon && std::abs(x.i - y.i) < epsilon;
    }

    // Absloute value of a complex number.
    // sqrt(a^2 + b^2)
    double abs() {
      return std::sqrt(norm());
    }

    // Argument of complex number.
    // atan2(b,a)
    double arg() {
      return std::atan2(i, r);
    }

    // a^2 + b^2
    double norm() {
      return r * r + i * i;
    }

    // Conjugate of a complex number.
    // a-bi
    complex conj() {
      return complex {r, -i};
    }

    // Absloute value of two parts, R and I.
    // abs(a) + abs(b)i
    complex pabs() {
      return complex {
        std::abs(r),
        std::abs(i)
      };
    }

    // projection of a complex number to Riemann sphere
    // z / (1 + norm(z))
    complex proj() {
      return *this / (norm() + 1);
    }

    // Polar to rectangular plane.
    // (r * cos(theta),r * sin(theta))
    static complex fromPolar(double r, double theta) {
      return complex {
        r * std::cos(theta),
        r * std::sin(theta)
      };
    }

    // Polar to rectangular plane.
    // (r * cos(theta),r * sin(theta))
    static complex fromPolar(struct polar p) {
      return fromPolar(p.r, p.theta);
    }

    // Rectangular plane to polar.
    // r = abs(z), θ = arg(z)
    static struct polar fromRectangular(double R, double I) {
      return (struct polar) {
        std::sqrt(R * R + I * I),
        std::atan2(I, R)
      };
    }

    // Rectangular plane to polar.
    // r = abs(z), θ = arg(z)
    static struct polar fromRectangular(complex r) {
      return fromRectangular(r.r, r.i);
    }

    // Inverse of complex number.
    // 1 / z
    complex inv() {
      const double n = norm();
      return complex {
        r  / n,
        -i / n
      };
    }

    // 0 subtract z
    // -z
    complex neg() {
      return complex {
        -r,
        -i
      };
    }


    // z multiplied by i.
    // iz
    complex muli() {
      return complex {-i, r};
    }

    // e^z, where z is a complex number.
    // e^a * (cos(b) + sin(b)i)
    complex exp() {
      const double e = std::exp(r);
      if (i == 0) return complex {e};
      return complex {
        e * std::cos(i),
        e * std::sin(i)
      };
    }

    // Natural logarithm of complex number. Base `e`.
    // 1/2 * ln(a^2 + b^2) + arg(z)i
    complex log() {
      return complex {
        0.5 * std::log(norm()),
        arg()
      };
    }

    // Natural logarithm of complex number. Same as `log()`.
    // 1/2 * ln(a^2 + b^2) + arg(z)i
    inline complex ln() {
      return log();
    }

    // Logarithm of base `10`.
    // ln(z) / ln(10)
    complex log10() {
      const double INV_LN10 = 0.43429448190325182;
      return log() * INV_LN10;
    }

    // Logarithm of base `2`.
    // ln(z) / ln(2)
    complex log2() {
      const double INV_LN2 = 1.442695040888963407;
      return log() * INV_LN2;
    }

    // Logarithm of base `z2`.
    // ln(z) / ln(z2)
    complex logZ(complex z2) {
      return log() / z2.log();
    }

    // `z` pow `z2`.
    // exp(ln(z) * z2)
    complex pow(complex z2) {
      return (ln() * z2).exp();
    }

    // z^2
    // z * z
    complex square() {
      return (*this) * (*this);
    }

    // Square root of complex number.
    // z ^ 0.5
    complex sqrt() {
      return (ln() * 0.5).exp();
    }

    // z2 th root of complex number.
    // z ^ (1 / z2)
    complex rootZ(complex z2) {
      return pow(z2.inv());
    }

    // sine.
    // sin(a)cosh(b)+cos(a)sinh(b)i
    complex sin() {
      return complex {
        std::sin(r) * std::cosh(i),
        std::cos(r) * std::sinh(i)
      };
    }

    // cosine.
    // cos(a)cosh(b)-sin(a)sinh(b)i
    complex cos() {
      return complex {
          std::cos(r) * std::cosh(i),
          -std::sin(r) * std::sinh(i)
      };
    }

    // tangent.
    // (tan(a)+tanh(b)i) / (1-tan(a)tanh(b)i)
    complex tan() {
      const double tan = std::tan(r);
      const double tanh = std::tanh(i);
      return complex {tan, tanh} / complex {1, -tan*tanh};
    }

    // arcsin. Inverse function of `sin`.
    // -i*ln(i*x + sqrt(1 - x^2))
    complex asin() {
      return ((1.0 - square()).sqrt() + muli()).log().muli().neg();
    }

    // arccos. Inverse function of `cos`.
    // 1/2 (π - 2*asin(z))
    complex acos() {
      return 0.5 * (3.1415926535897932 - 2.0 * asin());
    }

    // arctan. Inverse function of `tan`.
    // 1/2*i*ln(1 - iz) - 1/2*i*ln(1 + iz)
    complex atan() {
      complex iz = muli();
      complex hi {0,0.5};
      return hi * (1.0 - iz).log() - hi * (1.0 + iz).log();
    }

    // sinh. Hyperbolic function.
    // -e^(-z)/2 + e^z/2
    complex sinh() {
      return ((neg().exp() - exp()) / 2.0).neg();
    }

    // cosh. Hyperbolic function.
    // e^(-z)/2 + e^z/2
    complex cosh() {
      return (neg().exp() + exp()) / 2.0;
    }

    // tanh. Hyperbolic function.
    // -e^(-z)/(e^(-z) + e^z) + e^z/(e^(-z) + e^z)
    // (e^z - e^(-z)) / (e^z + e^(-z))
    complex tanh() {
      complex ez  = exp();
      complex enz = neg().exp();
      return (ez - enz) / (ez + enz);
    }

    // arsinh. Inverse hyperbolic function.
    // ln(z + sqrt(1 + z^2))
    complex asinh() {
      return ((1 + square()).sqrt() + *this).log();
    }

    // arcosh. Inverse hyperbolic function.
    // ln(z + sqrt(z - 1) * sqrt(z + 1))
    complex acosh() {
      return ((*this - 1).sqrt() * (*this + 1).sqrt() + *this).log();
    }

    // artanh. Inverse hyperbolic function.
    // 1/2 ln(1 + z) - 1/2 ln(1 - z)
    complex atanh() {
      return ((1 + *this).log() - (1 - *this).log()) / 2.0;
    }

    // round function
    complex round(double to = 1.0) {
      return complex {
        std::round(r * to) / to,
        std::round(i * to) / to
      };
    }

    // floor function
    complex floor(double to = 1.0) {
      return complex {
        std::floor(r * to) / to,
        std::floor(i * to) / to
      };
    }

    // ceil function
    complex ceil(double to = 1.0) {
      return complex {
        std::ceil(r * to) / to,
        std::ceil(i * to) / to
      };
    }

    // trunc function
    complex trunc() {
      return complex {
        std::trunc(r),
        std::trunc(i)
      };
    }

    // round off floating points errors.
    complex rndErr(double to = 1E14) {
      return round(to);
    }

    // gamma function
    // uses Lanczos approximation for gamma where N = 9 and g = 7.
    static complex gamma(complex z) {
      if (z == 0.0) return INFINITY;
      if (z == 1.0) return 1;
      if (z == 2.0) return 2;
      const double pi = 3.1415926535897932;
      const double s2pi = 2.50662827463100050241; // sqrt(2 * pi)
      const double g = 7.0 + 0.5;  // g = 7
      const double coeff[] = {
        0.99999999999980993,
        676.520368121885100,
        -1259.13921672240280,
        771.323428777653130,
        -176.615029162140590,
        12.5073432786869050,
        -0.13857109526572012,
        9.9843695780195716e-6,
        1.5056327351493116e-7
      };
      
      if (z.r < 0.5) {
        // reflection formula
        // Γ(z)Γ(1-z) = pi / sin(pi * z)
        // Γ(z) = pi / sin(pi * z) * Γ(1-z)
        return pi / (pi * z).sin() * complex::gamma(1 - z);
      } else {
        // Lanczos approximation
        // Γ(z) = sqrt(2pi) * (z+g+0.5)^(z+0.5) * exp(-(z+g+0.5)) * Ag(z)
        // where sqrt(2pi) = s2pi,
        // (z + g + 0.5) = t,
        // Ag(z) = Agz
        // final formula: Γ(z) = s2pi * t^(z+0.5) * exp(-t) * Agz
        z -= 1;

        complex Agz = coeff[0];
        for (int i = 1; i < 9; i++) Agz += coeff[i] / (z + i);

        complex t = g + z; // z + g + 0.5
        return s2pi * t.pow(z + 0.5) * (-t).exp() * Agz;
      }
    }

    // gamma function
    // uses Lanczos approximation for gamma where N = 9 and g = 7.
    complex gamma() {
      return complex::gamma(*this);
    }

    // distance of two complex numbers in rectangular plane.
    static double dist(complex x, complex y) {
      return std::sqrt((x.r-y.r)*(x.r-y.r) + (x.i-y.i)*(x.i-y.i));
    }

    // distance of two complex numbers in rectangular plane.
    double dist(complex y) {
      return dist(*this, y);
    }

  };

  // complex 1+0i
  complex CPLX_R{1,0};
  // complex 0+1i
  complex CPLX_I{0,1};
  // complex pi
  complex CPLX_PI{3.1415926535897932,0};
  // complex e
  complex CPLX_E{2.71828182845904523,0};
}