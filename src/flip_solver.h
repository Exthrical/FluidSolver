// 2D FLIP solver with PBF/XPBD volume preservation
//
// This refactor replaces ad-hoc min-separation and density relax passes
// with a robust Position-Based Fluids (PBF) / XPBD-style positional solver.
// Each substep performs:
//  - Neighbor search via grid bins using anisotropic spacing (hx, hy)
//  - Density constraint C_i = (rho_i/rho0 - 1) using 2D poly6 kernel
//  - Solve Lagrange multipliers lambda with XPBD compliance (alpha = compliance/dt^2)
//  - Position corrections: dx_i = (1/rho0) sum_j (lambda_i + lambda_j + s_corr) * gradW_ij
//  - Update particle velocities from position differences v += dx/dt to preserve momentum
// Kernels use world-space support radius r = kernelRadiusCells * max(hx, hy).
// Small per-iteration displacement cap improves stability.
// Diagnostics (density heatmap, correction vectors, simple metrics) are stored for UI.
#pragma once
#include <vector>
#include <cstdint>
#include <limits>
#include <random>
#include <algorithm>
#include "vec2.h"
#include "colliders.h"

struct Particle {
    Vec2 p;
    Vec2 v;
};

struct FlipParams {
    int nx = 128;
    int ny = 96;
    float dt = 1.0f / 60.0f;
    int substeps = 1;
    float gravity = -9.8f; // y-direction downwards
    float flipRatio = 0.7f; // 1 = pure FLIP, 0 = pure PIC
    int pressureIters = 60;
    int maxParticles = 200000;
    int particlesPerCell = 8; // for initial fill; helps stability at fine grids
    // Volume preservation (PBF/XPBD) controls
    // enableMinSeparation: interpreted as enabling PBF/XPBD volume preservation
    bool enableMinSeparation = true;
    // minSepIterations: number of solver iterations per substep
    int minSepIterations = 2;
    // minSepRelax: interpreted as XPBD compliance in [0..1]; higher = softer
    float minSepRelax = 0.0f;
    // densityStrength: interpreted as s_corr strength (artificial pressure) [0..0.5]
    float densityStrength = 0.1f;
    // densityBlur: interpreted as kernel support radius in cells (>=1)
    int densityBlur = 2;
};

class FlipSolver2D {
public:
    explicit FlipSolver2D(const FlipParams& params = {});

    void resize(int nx, int ny);

    void clearFluid();
    void clearColliders();

    void addFluidRect(const Vec2& minPt, const Vec2& maxPt, float jitter = 0.3f);
    void addFluidCircle(const Vec2& c, float r, float jitter = 0.3f);
    void addParticle(const Vec2& p, const Vec2& v = Vec2(0,0));

    void addColliderRect(const Vec2& minPt, const Vec2& maxPt);
    void addColliderCircle(const Vec2& c, float r);
    void eraseFluidCircle(const Vec2& c, float r);

    void step(float dt);

    // Accessors
    const std::vector<Particle>& particles() const { return particles_; }
    const ColliderSet& colliders() const { return colliders_; }
    int nx() const { return nx_; }
    int ny() const { return ny_; }
    float h() const { return h_; }
    float hx() const { return hx_; }
    float hy() const { return hy_; }
    float flipRatio() const { return params_.flipRatio; }
    void setFlipRatio(float r) { params_.flipRatio = r; }
    void setGravity(float g) { params_.gravity = g; }
    void setPressureIters(int n) { params_.pressureIters = n; }
    void setSubsteps(int n) { params_.substeps = std::max(1, n); }
    // Volume preservation setters (reinterpreted)
    void setMinSeparationEnabled(bool b) { params_.enableMinSeparation = b; }
    void setMinSeparationIters(int it) { params_.minSepIterations = std::max(0, it); }
    void setMinSeparationRelax(float r) { params_.minSepRelax = std::max(0.0f, std::min(1.0f, r)); }
    void setDensityRelaxEnabled(bool /*unused*/) { /* deprecated */ }
    void setDensityStrength(float k) { params_.densityStrength = std::max(0.0f, k); }
    void setDensityBlur(int r) { params_.densityBlur = std::max(1, r); }
    int particleCount() const { return (int)particles_.size(); }

    // For debugging / visualization
    const std::vector<float>& gridU() const { return u_; }
    const std::vector<float>& gridV() const { return v_; }
    // Canonical MAC sampler (used by G2P and debug rendering)
    Vec2 sampleMAC(const Vec2& pos) const;

    // Diagnostics
    const std::vector<float>& pbfDensities() const { return pbfDensity_; }
    const std::vector<Vec2>& pbfCorrections() const { return pbfDelta_; }
    float metricTopBandVelVar() const { return metricTopBandVelVar_; }
    float metricMaxOverlap() const { return metricMaxOverlap_; }
    float metricMeanDensityErr() const { return metricMeanDensityErr_; }
    float metricMomentumDrift() const { return metricMomentumDrift_; }

private:
    FlipParams params_;
    int nx_{0}, ny_{0};
    float h_{1.0f}; // representative cell size (min of hx,hy) for UI/brush spacing
    float hx_{1.0f}, hy_{1.0f}; // anisotropic spacings for x and y

    // Particles
    std::vector<Particle> particles_;
    std::mt19937 rng_{1234u};

    // Colliders
    ColliderSet colliders_;

    // Grid face velocities and helper arrays
    // u lives on (i in [0..nx], j in [0..ny-1]) at positions (i*h, (j+0.5)*h)
    // v lives on (i in [0..nx-1], j in [0..ny]) at positions ((i+0.5)*h, j*h)
    std::vector<float> u_, v_, uOld_, vOld_, uW_, vW_;
    // Cell-centered pressure and markers
    std::vector<float> p_, div_;
    std::vector<uint8_t> marker_; // 0=air, 1=fluid, 2=solid

    inline int idxU(int i, int j) const { return i + (nx_ + 1) * j; }
    inline int idxV(int i, int j) const { return i + nx_ * j; }
    inline int idxC(int i, int j) const { return i + nx_ * j; }

    void allocate();
    void clearGrid();
    void buildMarkers();
    void particlesToGrid();
    void applySolidBoundariesOnGrid();
    void computeDivergence();
    void pressureSolve();
    void subtractPressureGradient();
    void gridToParticles();
    void advectParticles(float dt);
    void pushOutOfColliders();
    void enforceParticleCollisions(Particle& p) const;
    void volumePreservationStep(float dt);

    // Interpolation helpers
    Vec2 sampleGridVelocity(const Vec2& pos, const std::vector<float>& u, const std::vector<float>& v) const;
    Vec2 sampleGridDelta(const Vec2& pos, const std::vector<float>& uNew, const std::vector<float>& vNew, const std::vector<float>& uOld, const std::vector<float>& vOld) const;

    // PBF storage for diagnostics
    std::vector<float> pbfDensity_;
    std::vector<float> pbfLambda_;
    std::vector<Vec2>  pbfDelta_;
    // Metrics
    float metricTopBandVelVar_{0.0f};
    float metricMaxOverlap_{0.0f};
    float metricMeanDensityErr_{0.0f};
    float metricMomentumDrift_{0.0f};
};
