// Minimal 2D vector utilities
#pragma once
#include <cmath>

struct Vec2 {
    float x{0}, y{0};
    Vec2() = default;
    Vec2(float _x, float _y) : x(_x), y(_y) {}
    Vec2 operator+(const Vec2& o) const { return Vec2(x + o.x, y + o.y); }
    Vec2 operator-(const Vec2& o) const { return Vec2(x - o.x, y - o.y); }
    Vec2 operator*(float s) const { return Vec2(x * s, y * s); }
    Vec2 operator/(float s) const { return Vec2(x / s, y / s); }
    Vec2& operator+=(const Vec2& o) { x += o.x; y += o.y; return *this; }
    Vec2& operator-=(const Vec2& o) { x -= o.x; y -= o.y; return *this; }
    Vec2& operator*=(float s) { x *= s; y *= s; return *this; }
    Vec2& operator/=(float s) { x /= s; y /= s; return *this; }
};

inline float dot(const Vec2& a, const Vec2& b) { return a.x*b.x + a.y*b.y; }
inline float length(const Vec2& v) { return std::sqrt(dot(v,v)); }
inline Vec2 normalize(const Vec2& v) {
    float len = length(v);
    if (len > 1e-8f) return v / len; else return Vec2(0,0);
}
inline Vec2 clampVec(const Vec2& v, const Vec2& lo, const Vec2& hi) {
    float x = v.x < lo.x ? lo.x : (v.x > hi.x ? hi.x : v.x);
    float y = v.y < lo.y ? lo.y : (v.y > hi.y ? hi.y : v.y);
    return Vec2(x,y);
}
inline float clampf(float x, float lo, float hi) {
    return x < lo ? lo : (x > hi ? hi : x);
}

