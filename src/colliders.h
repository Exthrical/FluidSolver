#pragma once
#include <vector>
#include <algorithm>
#include "vec2.h"

enum class ColliderType { Rectangle, Circle };

struct Collider {
    ColliderType type{ColliderType::Rectangle};
    // For Rect: a=min, b=max; For Circle: a=center, r=radius
    Vec2 a{0,0};
    Vec2 b{1,1};
    float r{0.1f};

    static Collider Rect(const Vec2& minPt, const Vec2& maxPt) {
        Collider c; c.type = ColliderType::Rectangle; c.a = minPt; c.b = maxPt; return c;
    }
    static Collider Circle(const Vec2& center, float rad) {
        Collider c; c.type = ColliderType::Circle; c.a = center; c.r = rad; return c;
    }

    bool contains(const Vec2& p) const {
        if (type == ColliderType::Rectangle) {
            return (p.x >= a.x && p.x <= b.x && p.y >= a.y && p.y <= b.y);
        } else {
            Vec2 d = p - a;
            return dot(d,d) <= r*r;
        }
    }

    float sdf(const Vec2& p) const {
        if (type == ColliderType::Rectangle) {
            // Signed distance to axis-aligned rectangle [a,b]
            Vec2 c = (a + b) * 0.5f;
            Vec2 e = (b - a) * 0.5f; // extents
            Vec2 d = Vec2(std::abs(p.x - c.x), std::abs(p.y - c.y)) - e;
            float outside = length(Vec2(std::max(d.x, 0.0f), std::max(d.y, 0.0f)));
            float inside = std::min(std::max(d.x, d.y), 0.0f);
            return outside + inside;
        } else {
            return length(p - a) - r;
        }
    }

    Vec2 normal(const Vec2& p) const {
        // Gradient of SDF, approximated for rect
        if (type == ColliderType::Rectangle) {
            // Approximate via finite difference
            const float eps = 1e-3f;
            float dx = sdf(Vec2(p.x + eps, p.y)) - sdf(Vec2(p.x - eps, p.y));
            float dy = sdf(Vec2(p.x, p.y + eps)) - sdf(Vec2(p.x, p.y - eps));
            return normalize(Vec2(dx, dy));
        } else {
            return normalize(p - a);
        }
    }
};

struct ColliderSet {
    std::vector<Collider> items;

    void clear() { items.clear(); }
    void addRect(const Vec2& minPt, const Vec2& maxPt) { items.push_back(Collider::Rect(minPt, maxPt)); }
    void addCircle(const Vec2& c, float r) { items.push_back(Collider::Circle(c, r)); }
    bool contains(const Vec2& p) const {
        for (const auto& c : items) if (c.contains(p)) return true; return false;
    }
    float sdf(const Vec2& p) const {
        float d = 1e9f;
        for (const auto& c : items) d = std::min(d, c.sdf(p));
        return d;
    }
    Vec2 normal(const Vec2& p) const {
        float dmin = 1e9f; Vec2 n{0,1};
        for (const auto& c : items) {
            float d = c.sdf(p);
            if (d < dmin) { dmin = d; n = c.normal(p); }
        }
        return n;
    }
};

