#pragma once
#include <vector>
#include "vec2.h"
#include "flip_solver.h"

struct RenderSettings {
    bool showGrid = false;
    bool showVel = false;
    int velStride = 8;
    bool showColliders = true;
    float particleSize = 3.0f; // pixels
    bool showParticleVel = false; // draw MAC-sampled velocity at particle positions
    bool showDensityHeatmap = false; // color particles by PBF density error
    bool showPbfCorrections = false; // draw PBF correction vectors
};

struct Viewport {
    int width{1280};
    int height{720};
    float scale{1.0f};
    Vec2 offset{0,0}; // screen space offset in pixels
};

class Renderer2D {
public:
    static void computeViewport(int fbWidth, int fbHeight, Viewport& vp);
    static Vec2 worldToScreen(const Vec2& wp, const Viewport& vp);
    static Vec2 screenToWorld(const Vec2& sp, const Viewport& vp);

    static void drawBackground();
    static void drawDomain(const Viewport& vp);
    static void drawParticles(const FlipSolver2D& solver, const Viewport& vp, const RenderSettings& rs);
    static void drawColliders(const FlipSolver2D& solver, const Viewport& vp, const RenderSettings& rs);
    static void drawGrid(const FlipSolver2D& solver, const Viewport& vp, const RenderSettings& rs);
    static void drawVelocity(const FlipSolver2D& solver, const Viewport& vp, const RenderSettings& rs);
    static void drawParticleVel(const FlipSolver2D& solver, const Viewport& vp, const RenderSettings& rs);
    static void drawPbfCorrections(const FlipSolver2D& solver, const Viewport& vp, const RenderSettings& rs);
};
