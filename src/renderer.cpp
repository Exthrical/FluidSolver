#include "renderer.h"
#include <cmath>
#include <algorithm>

#ifdef __APPLE__
#include <OpenGL/gl.h>
#else
#include <GL/gl.h>
#endif

static float clampf_local(float x, float lo, float hi) { return x < lo ? lo : (x > hi ? hi : x); }

void Renderer2D::computeViewport(int fbWidth, int fbHeight, Viewport& vp) {
    vp.width = fbWidth; vp.height = fbHeight;
    // Reserve a fixed sidebar on the left for UI and minimal margins
    const float kSidebar = 360.0f; // pixels
    const float kMargin = 16.0f;   // pixels
    float availW = std::max(0.0f, (float)fbWidth - kSidebar - 2.0f * kMargin);
    float availH = std::max(0.0f, (float)fbHeight - 2.0f * kMargin);
    float s = std::min(availW, availH);
    vp.scale = s;
    vp.offset = Vec2(kSidebar + kMargin + 0.5f * (availW - s), kMargin + 0.5f * (availH - s));
}

Vec2 Renderer2D::worldToScreen(const Vec2& wp, const Viewport& vp) {
    // world [0,1]^2 to screen pixel coords (origin top-left OpenGL screen y-up to be handled by matrix)
    return Vec2(vp.offset.x + wp.x * vp.scale, vp.offset.y + (1.0f - wp.y) * vp.scale);
}

Vec2 Renderer2D::screenToWorld(const Vec2& sp, const Viewport& vp) {
    float x = (sp.x - vp.offset.x) / vp.scale;
    float y = 1.0f - (sp.y - vp.offset.y) / vp.scale;
    return Vec2(x, y);
}

void Renderer2D::drawBackground() {
    glClearColor(0.08f, 0.09f, 0.11f, 1.0f);
    glClear(GL_COLOR_BUFFER_BIT);
}

void Renderer2D::drawDomain(const Viewport& vp) {
    // Draw domain border rectangle
    glMatrixMode(GL_PROJECTION); glLoadIdentity();
    glOrtho(0, vp.width, vp.height, 0, -1, 1);
    glMatrixMode(GL_MODELVIEW); glLoadIdentity();

    glLineWidth(2.0f);
    glColor3f(0.18f, 0.22f, 0.28f);
    Vec2 o = vp.offset; float s = vp.scale;
    glBegin(GL_LINE_LOOP);
    glVertex2f(o.x, o.y);
    glVertex2f(o.x + s, o.y);
    glVertex2f(o.x + s, o.y + s);
    glVertex2f(o.x, o.y + s);
    glEnd();
}

void Renderer2D::drawParticles(const FlipSolver2D& solver, const Viewport& vp, const RenderSettings& rs) {
    glEnable(GL_POINT_SMOOTH);
    glPointSize(clampf_local(rs.particleSize, 1.0f, 12.0f));
    glBegin(GL_POINTS);
    int n = solver.particleCount();
    const auto& parts = solver.particles();
    const auto& dens = solver.pbfDensities();
    for (int i = 0; i < n; ++i) {
        const Vec2 sp = worldToScreen(parts[i].p, vp);
        if (rs.showDensityHeatmap && (int)dens.size() == n) {
            float rho = dens[i];
            float err = std::fabs(rho - 1.0f); // rho0=1
            // Map error to red-yellow-white
            float r = clampf_local(0.5f + err, 0.2f, 1.0f);
            float g = clampf_local(0.2f + 0.8f * (1.0f - std::fabs(err - 0.5f)), 0.2f, 1.0f);
            float b = 0.2f;
            glColor3f(r, g, b);
        } else {
            float t = clampf_local(parts[i].p.y, 0.0f, 1.0f);
            // subtle blue gradient
            glColor3f(0.2f, 0.5f + 0.4f * (1.0f - t), 0.95f);
        }
        glVertex2f(sp.x, sp.y);
    }
    glEnd();
    glDisable(GL_POINT_SMOOTH);
}

void Renderer2D::drawColliders(const FlipSolver2D& solver, const Viewport& vp, const RenderSettings& rs) {
    if (!rs.showColliders) return;
    glColor4f(0.9f, 0.3f, 0.25f, 0.25f);
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

    for (const auto& col : solver.colliders().items) {
        if (col.type == ColliderType::Rectangle) {
            Vec2 s0 = worldToScreen(col.a, vp);
            Vec2 s1 = worldToScreen(col.b, vp);
            glBegin(GL_QUADS);
            glVertex2f(s0.x, s0.y);
            glVertex2f(s1.x, s0.y);
            glVertex2f(s1.x, s1.y);
            glVertex2f(s0.x, s1.y);
            glEnd();
            glColor3f(0.9f, 0.35f, 0.3f);
            glBegin(GL_LINE_LOOP);
            glVertex2f(s0.x, s0.y);
            glVertex2f(s1.x, s0.y);
            glVertex2f(s1.x, s1.y);
            glVertex2f(s0.x, s1.y);
            glEnd();
            glColor4f(0.9f, 0.3f, 0.25f, 0.25f);
        } else {
            // Approximate circle with poly
            const int seg = 48;
            glBegin(GL_TRIANGLE_FAN);
            Vec2 c = worldToScreen(col.a, vp);
            glVertex2f(c.x, c.y);
            for (int k = 0; k <= seg; ++k) {
                float ang = (float)k / seg * 6.2831853f;
                Vec2 wp = Vec2(col.a.x + std::cos(ang)*col.r, col.a.y + std::sin(ang)*col.r);
                Vec2 sp = worldToScreen(wp, vp);
                glVertex2f(sp.x, sp.y);
            }
            glEnd();
            glColor3f(0.9f, 0.35f, 0.3f);
            glBegin(GL_LINE_LOOP);
            for (int k = 0; k < seg; ++k) {
                float ang = (float)k / seg * 6.2831853f;
                Vec2 wp = Vec2(col.a.x + std::cos(ang)*col.r, col.a.y + std::sin(ang)*col.r);
                Vec2 sp = worldToScreen(wp, vp);
                glVertex2f(sp.x, sp.y);
            }
            glEnd();
        }
    }
    glDisable(GL_BLEND);
}

void Renderer2D::drawGrid(const FlipSolver2D& solver, const Viewport& vp, const RenderSettings& rs) {
    if (!rs.showGrid) return;
    glColor3f(0.15f, 0.17f, 0.2f);
    glLineWidth(1.0f);
    int nx = solver.nx(); int ny = solver.ny();
    for (int i = 0; i <= nx; ++i) {
        Vec2 s0 = worldToScreen(Vec2(i/(float)nx, 0), vp);
        Vec2 s1 = worldToScreen(Vec2(i/(float)nx, 1), vp);
        glBegin(GL_LINES); glVertex2f(s0.x, s0.y); glVertex2f(s1.x, s1.y); glEnd();
    }
    for (int j = 0; j <= ny; ++j) {
        Vec2 s0 = worldToScreen(Vec2(0, j/(float)ny), vp);
        Vec2 s1 = worldToScreen(Vec2(1, j/(float)ny), vp);
        glBegin(GL_LINES); glVertex2f(s0.x, s0.y); glVertex2f(s1.x, s1.y); glEnd();
    }
}

void Renderer2D::drawVelocity(const FlipSolver2D& solver, const Viewport& vp, const RenderSettings& rs) {
    if (!rs.showVel) return;
    glColor3f(0.3f, 0.9f, 0.4f);
    glLineWidth(1.0f);
    int nx = solver.nx(); int ny = solver.ny();
    int stride = std::max(1, rs.velStride);
    // Sample at cell centers for visualization
    for (int j = stride/2; j < ny; j += stride) {
        for (int i = stride/2; i < nx; i += stride) {
            Vec2 wp = Vec2((i + 0.5f) / (float)nx, (j + 0.5f) / (float)ny);
            Vec2 vel = solver.sampleMAC(wp) * 0.02f; // scale for display
            Vec2 s0 = worldToScreen(wp, vp);
            Vec2 s1 = worldToScreen(Vec2(wp.x + vel.x, wp.y + vel.y), vp);
            glBegin(GL_LINES); glVertex2f(s0.x, s0.y); glVertex2f(s1.x, s1.y); glEnd();
        }
    }
}

void Renderer2D::drawParticleVel(const FlipSolver2D& solver, const Viewport& vp, const RenderSettings& rs) {
    if (!rs.showParticleVel) return;
    glColor3f(0.95f, 0.8f, 0.3f);
    glLineWidth(1.0f);
    const auto& parts = solver.particles();
    const float scale = 0.02f; // visual scale
    for (const auto& prt : parts) {
        Vec2 v = solver.sampleMAC(prt.p) * scale;
        Vec2 s0 = worldToScreen(prt.p, vp);
        Vec2 s1 = worldToScreen(Vec2(prt.p.x + v.x, prt.p.y + v.y), vp);
        glBegin(GL_LINES); glVertex2f(s0.x, s0.y); glVertex2f(s1.x, s1.y); glEnd();
    }
}

void Renderer2D::drawPbfCorrections(const FlipSolver2D& solver, const Viewport& vp, const RenderSettings& rs) {
    if (!rs.showPbfCorrections) return;
    glColor3f(0.9f, 0.9f, 0.2f);
    glLineWidth(1.0f);
    const auto& parts = solver.particles();
    const auto& corr = solver.pbfCorrections();
    const float scale = 1.0f; // already in world units per-iteration; small
    int n = solver.particleCount();
    if ((int)corr.size() != n) return;
    for (int i = 0; i < n; ++i) {
        Vec2 d = corr[i] * scale;
        if (d.x == 0.0f && d.y == 0.0f) continue;
        Vec2 s0 = worldToScreen(parts[i].p, vp);
        Vec2 s1 = worldToScreen(Vec2(parts[i].p.x + d.x, parts[i].p.y + d.y), vp);
        glBegin(GL_LINES); glVertex2f(s0.x, s0.y); glVertex2f(s1.x, s1.y); glEnd();
    }
}
