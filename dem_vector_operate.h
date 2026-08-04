#pragma once

struct Vec3 { double x, y, z; };

//ÖØÔØËã·û
inline Vec3 operator+(const Vec3& a, const Vec3& b) { return { a.x + b.x, a.y + b.y, a.z + b.z }; }
inline Vec3 operator-(const Vec3& a, const Vec3& b) { return { a.x - b.x, a.y - b.y, a.z - b.z }; }
inline Vec3 operator*(double s, const Vec3& v) { return { s * v.x, s * v.y, s * v.z }; }
inline Vec3 operator*(const Vec3& v, double s) { return s * v; }
inline Vec3 operator/(const Vec3& v, double s) { return { v.x / s, v.y / s, v.z / s }; }

inline double dot(const Vec3& a, const Vec3& b) { return a.x * b.x + a.y * b.y + a.z * b.z; }

inline double len(const Vec3& a){ return sqrt(dot(a, a)); }

inline Vec3 cross(const Vec3& a, const Vec3& b)
{
	return { a.y * b.z - a.z * b.y,
			  a.z * b.x - a.x * b.z,
			  a.x * b.y - a.y * b.x };
}

inline Vec3 norm(const Vec3& a){ return a / len(a); }

