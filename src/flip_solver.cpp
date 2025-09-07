#include "flip_solver.h"
#include <cassert>
#include <cmath>

FlipSolver2D::FlipSolver2D(const FlipParams& params) : params_(params) {
    resize(params_.nx, params_.ny);
}

void FlipSolver2D::resize(int nx, int ny) {
    nx_ = std::max(4, nx);
    ny_ = std::max(4, ny);
    hx_ = 1.0f / (float)nx_;
    hy_ = 1.0f / (float)ny_;
    h_ = std::min(hx_, hy_);
    allocate();
}

void FlipSolver2D::allocate() {
    u_.assign((nx_ + 1) * ny_, 0.0f);
    v_.assign(nx_ * (ny_ + 1), 0.0f);
    uOld_ = u_;
    vOld_ = v_;
    uW_.assign((nx_ + 1) * ny_, 0.0f);
    vW_.assign(nx_ * (ny_ + 1), 0.0f);
    p_.assign(nx_ * ny_, 0.0f);
    div_.assign(nx_ * ny_, 0.0f);
    marker_.assign(nx_ * ny_, 0);
}

void FlipSolver2D::clearGrid() {
    std::fill(u_.begin(), u_.end(), 0.0f);
    std::fill(v_.begin(), v_.end(), 0.0f);
    std::fill(uW_.begin(), uW_.end(), 0.0f);
    std::fill(vW_.begin(), vW_.end(), 0.0f);
}

void FlipSolver2D::clearFluid() {
    particles_.clear();
}

void FlipSolver2D::clearColliders() {
    colliders_.clear();
}

void FlipSolver2D::addFluidRect(const Vec2& minPt, const Vec2& maxPt, float jitter) {
    std::uniform_real_distribution<float> uni(-jitter, jitter);
    int count = 0;
    float spacing = h_ / std::sqrt(std::max(1, params_.particlesPerCell));
    for (float y = minPt.y; y <= maxPt.y; y += spacing) {
        for (float x = minPt.x; x <= maxPt.x; x += spacing) {
            if ((int)particles_.size() >= params_.maxParticles) break;
            Vec2 p = Vec2(x, y) + Vec2(uni(rng_), uni(rng_)) * spacing;
            if (p.x < 0.0f || p.x > 1.0f || p.y < 0.0f || p.y > 1.0f) continue;
            if (colliders_.contains(p)) continue;
            particles_.push_back({p, Vec2(0,0)});
            count++;
        }
    }
}

void FlipSolver2D::addFluidCircle(const Vec2& c, float r, float jitter) {
    std::uniform_real_distribution<float> uni(-jitter, jitter);
    float spacing = h_ / std::sqrt(std::max(1, params_.particlesPerCell));
    int count = 0;
    int steps = std::max(1, (int)(2*r/spacing));
    for (int j = -steps; j <= steps; ++j) {
        for (int i = -steps; i <= steps; ++i) {
            if ((int)particles_.size() >= params_.maxParticles) break;
            Vec2 p = c + Vec2(i * spacing, j * spacing) + Vec2(uni(rng_), uni(rng_)) * spacing;
            if (dot(p - c, p - c) > r*r) continue;
            if (p.x < 0.0f || p.x > 1.0f || p.y < 0.0f || p.y > 1.0f) continue;
            if (colliders_.contains(p)) continue;
            particles_.push_back({p, Vec2(0,0)});
            count++;
        }
    }
}

void FlipSolver2D::addParticle(const Vec2& p, const Vec2& v) {
    if ((int)particles_.size() >= params_.maxParticles) return;
    if (p.x < 0.0f || p.x > 1.0f || p.y < 0.0f || p.y > 1.0f) return;
    if (colliders_.contains(p)) return;
    particles_.push_back({p, v});
}

void FlipSolver2D::addColliderRect(const Vec2& minPt, const Vec2& maxPt) {
    colliders_.addRect(minPt, maxPt);
}

void FlipSolver2D::addColliderCircle(const Vec2& c, float r) {
    colliders_.addCircle(c, r);
}

void FlipSolver2D::eraseFluidCircle(const Vec2& c, float r) {
    const float rr = r * r;
    particles_.erase(std::remove_if(particles_.begin(), particles_.end(), [&](const Particle& p){
        return dot(p.p - c, p.p - c) <= rr;
    }), particles_.end());
}

void FlipSolver2D::step(float dt) {
    if (particles_.empty()) return;
    // Adaptive substepping based on CFL to avoid particles crossing many cells per step
    int sub = std::max(1, params_.substeps);
    float maxSpeed = 0.0f;
    for (const auto& p : particles_) {
        float s = std::sqrt(p.v.x * p.v.x + p.v.y * p.v.y);
        if (s > maxSpeed) maxSpeed = s;
    }
    const float cfl = 1.0f; // target <= 1 cell per substep
    if (maxSpeed > 0.0f) {
        int need = (int)std::ceil((maxSpeed * dt) / (cfl * h_));
        sub = std::max(sub, std::max(1, need));
        sub = std::min(sub, 128); // hard cap
    }
    float sdt = dt / (float)sub;
    for (int k = 0; k < sub; ++k) {
        // Gravity to particles (more stable than grid add for FLIP)
        for (auto& p : particles_) {
            p.v.y += params_.gravity * sdt;
        }
        // Advect positions only; do not clamp before separation
        advectParticles(sdt);
        // Push particles out of colliders (but do not clamp to domain yet)
        pushOutOfColliders();
        // Particle-based volume preservation (PBF/XPBD)
        if (params_.enableMinSeparation) volumePreservationStep(sdt);
        // After separation, clamp to domain and colliders so P2G sees valid positions
        for (auto& p : particles_) enforceParticleCollisions(p);
        clearGrid();
        buildMarkers();
        particlesToGrid();
        applySolidBoundariesOnGrid();
        uOld_ = u_; vOld_ = v_;
        computeDivergence();
        pressureSolve();
        subtractPressureGradient();
        gridToParticles();
        // Old density relax removed; velocity already updated via PBF position solve
        // Final clamp for this substep to keep particles within domain even if density relax disabled
        for (auto& p : particles_) enforceParticleCollisions(p);
    }
}

void FlipSolver2D::buildMarkers() {
    std::fill(marker_.begin(), marker_.end(), 0);
    // Mark solid cells from colliders (domain boundaries handled via velocity BCs)
    for (int j = 0; j < ny_; ++j) {
        for (int i = 0; i < nx_; ++i) {
            Vec2 c = Vec2((i + 0.5f) * hx_, (j + 0.5f) * hy_);
            if (colliders_.sdf(c) < 0.0f) marker_[idxC(i,j)] = 2;
        }
    }
    // Mark fluid cells where particles reside
    for (const auto& p : particles_) {
        int i = (int)std::floor(p.p.x / hx_);
        int j = (int)std::floor(p.p.y / hy_);
        if (i >= 0 && j >= 0 && i < nx_ && j < ny_) {
            if (marker_[idxC(i,j)] != 2) marker_[idxC(i,j)] = 1; // fluid unless solid
        }
    }
    // Dilate fluid region by 2 cells to avoid under-sampling holes on fine grids
    for (int pass = 0; pass < 2; ++pass) {
        std::vector<uint8_t> mark2 = marker_;
        for (int j = 0; j < ny_; ++j) {
            for (int i = 0; i < nx_; ++i) {
                int c = idxC(i,j);
                if (marker_[c] != 0) continue; // keep air only
                bool nearFluid = false;
                if (i > 0) nearFluid |= (marker_[idxC(i-1,j)] == 1);
                if (i < nx_-1) nearFluid |= (marker_[idxC(i+1,j)] == 1);
                if (j > 0) nearFluid |= (marker_[idxC(i,j-1)] == 1);
                if (j < ny_-1) nearFluid |= (marker_[idxC(i,j+1)] == 1);
                if (nearFluid) mark2[c] = 1;
            }
        }
        // Do not overwrite solids
        for (int j = 0; j < ny_; ++j) for (int i = 0; i < nx_; ++i) if (marker_[idxC(i,j)] == 2) mark2[idxC(i,j)] = 2;
        marker_ = std::move(mark2);
    }
}

void FlipSolver2D::particlesToGrid() {
    // Accumulate particle velocities to face centers with bilinear weights
    for (const auto& prt : particles_) {
        const Vec2& P = prt.p; const Vec2& V = prt.v;
        // U faces at (i*hx, (j+0.5)*hy)
        float x = P.x / hx_; float y = P.y / hy_ - 0.5f;
        int i0 = (int)std::floor(x); int j0 = (int)std::floor(y);
        float fx = x - i0; float fy = y - j0;
        for (int dj = 0; dj <= 1; ++dj) for (int di = 0; di <= 1; ++di) {
            int i = i0 + di; int j = j0 + dj;
            if (i >= 0 && i <= nx_ && j >= 0 && j < ny_) {
                float w = (di ? fx : 1.0f - fx) * (dj ? fy : 1.0f - fy);
                int idx = idxU(i,j);
                u_[idx] += V.x * w; uW_[idx] += w;
            }
        }
        // V faces at ((i+0.5)*hx, j*hy)
        x = P.x / hx_ - 0.5f; y = P.y / hy_;
        i0 = (int)std::floor(x); j0 = (int)std::floor(y);
        fx = x - i0; fy = y - j0;
        for (int dj = 0; dj <= 1; ++dj) for (int di = 0; di <= 1; ++di) {
            int i = i0 + di; int j = j0 + dj;
            if (i >= 0 && i < nx_ && j >= 0 && j <= ny_) {
                float w = (di ? fx : 1.0f - fx) * (dj ? fy : 1.0f - fy);
                int idx = idxV(i,j);
                v_[idx] += V.y * w; vW_[idx] += w;
            }
        }
    }
    // Normalize by weights
    for (size_t i = 0; i < u_.size(); ++i) if (uW_[i] > 0.0f) u_[i] /= uW_[i];
    for (size_t i = 0; i < v_.size(); ++i) if (vW_[i] > 0.0f) v_[i] /= vW_[i];
}

void FlipSolver2D::applySolidBoundariesOnGrid() {
    // Zero out normal component at solid faces and domain borders
    for (int j = 0; j < ny_; ++j) {
        for (int i = 0; i <= nx_; ++i) {
            int idx = idxU(i,j);
            bool solid = (i == 0 || i == nx_);
            Vec2 fc = Vec2(i * hx_, (j + 0.5f) * hy_);
            if (colliders_.sdf(fc) < 0.0f) solid = true;
            if (solid) u_[idx] = 0.0f;
        }
    }
    for (int j = 0; j <= ny_; ++j) {
        for (int i = 0; i < nx_; ++i) {
            int idx = idxV(i,j);
            bool solid = (j == 0 || j == ny_);
            Vec2 fc = Vec2((i + 0.5f) * hx_, j * hy_);
            if (colliders_.sdf(fc) < 0.0f) solid = true;
            if (solid) v_[idx] = 0.0f;
        }
    }
}

void FlipSolver2D::computeDivergence() {
    for (int j = 0; j < ny_; ++j) {
        for (int i = 0; i < nx_; ++i) {
            int c = idxC(i,j);
            if (marker_[c] == 2) { div_[c] = 0.0f; continue; }
            float uR = u_[idxU(i+1, j)];
            float uL = u_[idxU(i,   j)];
            float vT = v_[idxV(i,   j+1)];
            float vB = v_[idxV(i,   j  )];
            div_[c] = (uR - uL) / hx_ + (vT - vB) / hy_;
        }
    }
}

void FlipSolver2D::pressureSolve() {
    // Gauss-Seidel iterations for anisotropic grid spacing
    // Discrete Poisson: (pL - pC)/hx^2 + (pR - pC)/hx^2 + (pB - pC)/hy^2 + (pT - pC)/hy^2 = div
    // Handle boundary/air/solid by weighting per-direction.
    const float wx = 1.0f / (hx_ * hx_);
    const float wy = 1.0f / (hy_ * hy_);
    for (int it = 0; it < params_.pressureIters; ++it) {
        for (int j = 0; j < ny_; ++j) {
            for (int i = 0; i < nx_; ++i) {
                int c = idxC(i,j);
                if (marker_[c] == 2 || marker_[c] == 0) { p_[c] = 0.0f; continue; }
                float sum = 0.0f; float diag = 0.0f;
                // Left
                if (i > 0) {
                    uint8_t m = marker_[idxC(i-1,j)];
                    if (m == 1) { sum += wx * p_[idxC(i-1,j)]; diag += wx; }
                    else if (m == 0) { diag += wx; } // air: Dirichlet p=0 contributes to diag only
                    else { /* solid: Neumann -> skip (no contribution) */ }
                } else { /* domain boundary: Neumann -> skip */ }
                // Right
                if (i < nx_-1) {
                    uint8_t m = marker_[idxC(i+1,j)];
                    if (m == 1) { sum += wx * p_[idxC(i+1,j)]; diag += wx; }
                    else if (m == 0) { diag += wx; }
                    else { }
                } else { }
                // Bottom
                if (j > 0) {
                    uint8_t m = marker_[idxC(i,j-1)];
                    if (m == 1) { sum += wy * p_[idxC(i,j-1)]; diag += wy; }
                    else if (m == 0) { diag += wy; }
                    else { }
                } else { }
                // Top
                if (j < ny_-1) {
                    uint8_t m = marker_[idxC(i,j+1)];
                    if (m == 1) { sum += wy * p_[idxC(i,j+1)]; diag += wy; }
                    else if (m == 0) { diag += wy; }
                    else { }
                } else { }

                if (diag > 0.0f) p_[c] = (sum - div_[c]) / diag; else p_[c] = 0.0f;
            }
        }
    }
}

void FlipSolver2D::subtractPressureGradient() {
    // U faces
    for (int j = 0; j < ny_; ++j) {
        for (int i = 1; i < nx_; ++i) { // interior faces
            int cR = idxC(i, j);
            int cL = idxC(i-1, j);
            if (marker_[cR] == 2 && marker_[cL] == 2) continue; // solid-solid
            float pR = (marker_[cR] == 1) ? p_[cR] : 0.0f;
            float pL = (marker_[cL] == 1) ? p_[cL] : 0.0f;
            // Neumann at solids: copy neighbor pressure to avoid artificial gradient
            if (marker_[cR] == 2 && marker_[cL] == 1) pR = pL;
            if (marker_[cL] == 2 && marker_[cR] == 1) pL = pR;
            u_[idxU(i,j)] -= (pR - pL) / hx_;
        }
    }
    // V faces
    for (int j = 1; j < ny_; ++j) {
        for (int i = 0; i < nx_; ++i) {
            int cT = idxC(i, j);
            int cB = idxC(i, j-1);
            if (marker_[cT] == 2 && marker_[cB] == 2) continue;
            float pT = (marker_[cT] == 1) ? p_[cT] : 0.0f;
            float pB = (marker_[cB] == 1) ? p_[cB] : 0.0f;
            if (marker_[cT] == 2 && marker_[cB] == 1) pT = pB;
            if (marker_[cB] == 2 && marker_[cT] == 1) pB = pT;
            v_[idxV(i,j)] -= (pT - pB) / hy_;
        }
    }
    applySolidBoundariesOnGrid();
}

Vec2 FlipSolver2D::sampleMAC(const Vec2& P) const {
    return sampleGridVelocity(P, u_, v_);
}

Vec2 FlipSolver2D::sampleGridVelocity(const Vec2& P, const std::vector<float>& u, const std::vector<float>& v) const {
    // Sample u
    float x = P.x / hx_; float y = P.y / hy_ - 0.5f;
    int i0 = (int)std::floor(x); int j0 = (int)std::floor(y);
    float fx = x - i0; float fy = y - j0;
    float ux = 0.0f; float sumWu = 0.0f;
    for (int dj = 0; dj <= 1; ++dj) for (int di = 0; di <= 1; ++di) {
        int i = i0 + di; int j = j0 + dj;
        if (i >= 0 && i <= nx_ && j >= 0 && j < ny_) {
            float w = (di ? fx : 1.0f - fx) * (dj ? fy : 1.0f - fy);
            ux += u[idxU(i,j)] * w; sumWu += w;
        }
    }
    // Sample v
    x = P.x / hx_ - 0.5f; y = P.y / hy_;
    i0 = (int)std::floor(x); j0 = (int)std::floor(y);
    fx = x - i0; fy = y - j0;
    float vy = 0.0f; float sumWv = 0.0f;
    for (int dj = 0; dj <= 1; ++dj) for (int di = 0; di <= 1; ++di) {
        int i = i0 + di; int j = j0 + dj;
        if (i >= 0 && i < nx_ && j >= 0 && j <= ny_) {
            float w = (di ? fx : 1.0f - fx) * (dj ? fy : 1.0f - fy);
            vy += v[idxV(i,j)] * w; sumWv += w;
        }
    }
    if (sumWu > 0.0f) ux /= sumWu; if (sumWv > 0.0f) vy /= sumWv;
    return Vec2(ux, vy);
}

Vec2 FlipSolver2D::sampleGridDelta(const Vec2& pos, const std::vector<float>& uNew, const std::vector<float>& vNew, const std::vector<float>& uOld, const std::vector<float>& vOld) const {
    Vec2 vN = sampleGridVelocity(pos, uNew, vNew);
    Vec2 vO = sampleGridVelocity(pos, uOld, vOld);
    return vN - vO;
}

void FlipSolver2D::gridToParticles() {
    for (auto& p : particles_) {
        Vec2 vpic = sampleGridVelocity(p.p, u_, v_);
        Vec2 vdelta = sampleGridDelta(p.p, u_, v_, uOld_, vOld_);
        Vec2 vflip = p.v + vdelta;
        // Adaptive FLIP-PIC blending to tame excessive FLIP energy injection
        float magPic = std::sqrt(vpic.x*vpic.x + vpic.y*vpic.y) + 1e-6f;
        float magDel = std::sqrt(vdelta.x*vdelta.x + vdelta.y*vdelta.y);
        float ratio = magDel / magPic; // how large the update is relative to PIC
        // If the grid correction is very large, lean towards PIC.
        float t = (ratio - 1.5f) / (3.0f - 1.5f); // 0 at 1.5x, 1 at 3x
        if (t < 0.0f) t = 0.0f; if (t > 1.0f) t = 1.0f;
        float w = params_.flipRatio * (1.0f - 0.8f * t); // reduce up to 80%
        if (w < 0.0f) w = 0.0f; if (w > 1.0f) w = 1.0f;
        p.v = vpic * (1.0f - w) + vflip * w;
    }
}

void FlipSolver2D::advectParticles(float dt) {
    for (auto& p : particles_) {
        // Simple explicit Euler
        p.p += p.v * dt;
    }
}

void FlipSolver2D::pushOutOfColliders() {
    for (auto& prt : particles_) {
        float d = colliders_.sdf(prt.p);
        if (d < 0.0f) {
            Vec2 n = colliders_.normal(prt.p);
            prt.p -= n * d; // push out
            float vn = dot(prt.v, n);
            if (vn < 0.0f) prt.v -= n * vn; // kill inward normal component
        }
    }
}

void FlipSolver2D::enforceParticleCollisions(Particle& prt) const {
    // Domain walls
    const float pad = 1e-4f;
    if (prt.p.x < pad) { prt.p.x = pad; if (prt.v.x < 0) prt.v.x = 0; }
    if (prt.p.x > 1.0f - pad) { prt.p.x = 1.0f - pad; if (prt.v.x > 0) prt.v.x = 0; }
    if (prt.p.y < pad) { prt.p.y = pad; if (prt.v.y < 0) prt.v.y = 0; }
    if (prt.p.y > 1.0f - pad) { prt.p.y = 1.0f - pad; if (prt.v.y > 0) prt.v.y = 0; }

    // Colliders: project out if inside
    float d = colliders_.sdf(prt.p);
    if (d < 0.0f) {
        Vec2 n = colliders_.normal(prt.p);
        prt.p -= n * d; // push out
        float vn = dot(prt.v, n);
        if (vn < 0.0f) prt.v -= n * vn; // kill inward normal component
    }
}

// Position-Based Fluids (PBF/XPBD) density constraint solve
void FlipSolver2D::volumePreservationStep(float dt) {
    const int n = (int)particles_.size();
    if (n == 0) return;

    // Ensure diagnostic storage is sized
    if ((int)pbfDensity_.size() != n) pbfDensity_.assign(n, 0.0f);
    if ((int)pbfLambda_.size()  != n) pbfLambda_.assign(n, 0.0f);
    if ((int)pbfDelta_.size()   != n) pbfDelta_.assign(n, Vec2(0,0));

    // Kernel support radius (world units)
    const int supportCells = std::max(1, params_.densityBlur); // reinterpretation
    const float hKer = std::max(hx_, hy_) * (float)supportCells;
    const float hKer2 = hKer * hKer;

    // Poly6 (2D) for density, Spiky (2D) for gradient
    auto W_poly6 = [&](float r2) -> float {
        if (r2 >= hKer2) return 0.0f;
        float x = hKer2 - r2;
        const float k = 4.0f / (3.14159265358979323846f * hKer2 * hKer2 * hKer2 * hKer2); // 4/(pi h^8)
        return k * x * x * x;
    };
    auto gradW_spiky = [&](const Vec2& rij) -> Vec2 {
        float r2 = rij.x*rij.x + rij.y*rij.y;
        if (r2 <= 1e-12f || r2 >= hKer2) return Vec2(0,0);
        float r = std::sqrt(r2);
        Vec2 dir = rij / r;
        // 2D spiky grad: -30/(pi h^5) * (h - r)^2 * dir
        const float k = -30.0f / (3.14159265358979323846f * hKer*hKer*hKer*hKer*hKer);
        float s = (hKer - r);
        float mag = k * s * s;
        return dir * mag;
    };

    // Build spatial bins for neighbors
    std::vector<std::vector<int>> bins(nx_ * ny_);
    auto binIndex = [&](const Vec2& p) -> int {
        int ix = (int)std::floor(p.x / hx_);
        int iy = (int)std::floor(p.y / hy_);
        if (ix < 0 || iy < 0 || ix >= nx_ || iy >= ny_) return -1;
        return idxC(ix, iy);
    };
    for (int i = 0; i < n; ++i) {
        int b = binIndex(particles_[i].p);
        if (b >= 0) bins[b].push_back(i);
    }

    // Neighbor cell range based on support radius
    const int rx = std::max(1, (int)std::ceil(hKer / hx_));
    const int ry = std::max(1, (int)std::ceil(hKer / hy_));

    // Rest density (dimensionless due to normalized kernel, mass=1)
    const float rho0 = 1.0f;
    const float eps = 1e-6f; // to avoid division by zero in lambda
    const float compliance = params_.minSepRelax; // [0..1]
    const float alpha = (compliance > 0.0f && dt > 0.0f) ? (compliance / (dt*dt)) : 0.0f;
    const float scorrK = params_.densityStrength; // small artificial pressure to reduce clumping
    const float scorrN = 4.0f; // exponent
    const float wq = W_poly6(0.25f * hKer2); // reference weight at q = 0.5 h

    // Save old positions for velocity update and compute momentum drift
    std::vector<Vec2> posPrev(n);
    for (int i = 0; i < n; ++i) posPrev[i] = particles_[i].p;

    // Solver iterations
    const int iters = std::max(0, params_.minSepIterations);
    for (int it = 0; it < iters; ++it) {
        // 1) Compute densities and lambdas
        std::fill(pbfDensity_.begin(), pbfDensity_.end(), 0.0f);
        std::fill(pbfLambda_.begin(),  pbfLambda_.end(),  0.0f);
        std::fill(pbfDelta_.begin(),   pbfDelta_.end(),   Vec2(0,0));

        // Densities
        for (int jy = 0; jy < ny_; ++jy) {
            for (int ix = 0; ix < nx_; ++ix) {
                int c = idxC(ix, jy);
                const auto& cell = bins[c];
                if (cell.empty()) continue;
                // For each particle in cell, search neighbor cells within rx,ry inclusively
                for (int aIdx = 0; aIdx < (int)cell.size(); ++aIdx) {
                    int i = cell[aIdx];
                    const Vec2& xi = particles_[i].p;
                    float rho = 0.0f;
                    for (int dj = -ry; dj <= ry; ++dj) {
                        int nyi = jy + dj; if (nyi < 0 || nyi >= ny_) continue;
                        for (int di = -rx; di <= rx; ++di) {
                            int nxi = ix + di; if (nxi < 0 || nxi >= nx_) continue;
                            const auto& nb = bins[idxC(nxi, nyi)];
                            for (int bIdx = 0; bIdx < (int)nb.size(); ++bIdx) {
                                int j = nb[bIdx];
                                Vec2 rij = xi - particles_[j].p;
                                float r2 = rij.x*rij.x + rij.y*rij.y;
                                rho += W_poly6(r2);
                            }
                        }
                    }
                    pbfDensity_[i] = rho; // m=1
                }
            }
        }

        // Lambdas
        for (int jy = 0; jy < ny_; ++jy) {
            for (int ix = 0; ix < nx_; ++ix) {
                int c = idxC(ix, jy);
                const auto& cell = bins[c]; if (cell.empty()) continue;
                for (int aIdx = 0; aIdx < (int)cell.size(); ++aIdx) {
                    int i = cell[aIdx];
                    const Vec2& xi = particles_[i].p;
                    float Ci = pbfDensity_[i] / rho0 - 1.0f;
                    // Compute sum of gradient norms squared
                    Vec2 gradCi_i(0,0);
                    float sumGrad2 = 0.0f;
                    for (int dj = -ry; dj <= ry; ++dj) {
                        int nyi = jy + dj; if (nyi < 0 || nyi >= ny_) continue;
                        for (int di = -rx; di <= rx; ++di) {
                            int nxi = ix + di; if (nxi < 0 || nxi >= nx_) continue;
                            const auto& nb = bins[idxC(nxi, nyi)];
                            for (int bIdx = 0; bIdx < (int)nb.size(); ++bIdx) {
                                int j = nb[bIdx];
                                Vec2 rij = xi - particles_[j].p;
                                Vec2 gradW = gradW_spiky(rij) / rho0; // m=1
                                if (j == i) continue;
                                sumGrad2 += gradW.x*gradW.x + gradW.y*gradW.y;
                                gradCi_i -= gradW; // negative sum of neighbor grads
                            }
                        }
                    }
                    sumGrad2 += gradCi_i.x*gradCi_i.x + gradCi_i.y*gradCi_i.y;
                    float denom = sumGrad2 + eps + alpha;
                    pbfLambda_[i] = (denom > 0.0f) ? (-Ci / denom) : 0.0f;
                }
            }
        }

        // 2) Compute position corrections
        for (int jy = 0; jy < ny_; ++jy) {
            for (int ix = 0; ix < nx_; ++ix) {
                int c = idxC(ix, jy);
                const auto& cell = bins[c]; if (cell.empty()) continue;
                for (int aIdx = 0; aIdx < (int)cell.size(); ++aIdx) {
                    int i = cell[aIdx];
                    const Vec2& xi = particles_[i].p;
                    Vec2 delta(0,0);
                    for (int dj = -ry; dj <= ry; ++dj) {
                        int nyi = jy + dj; if (nyi < 0 || nyi >= ny_) continue;
                        for (int di = -rx; di <= rx; ++di) {
                            int nxi = ix + di; if (nxi < 0 || nxi >= nx_) continue;
                            const auto& nb = bins[idxC(nxi, nyi)];
                            for (int bIdx = 0; bIdx < (int)nb.size(); ++bIdx) {
                                int j = nb[bIdx];
                                if (j == i) continue;
                                Vec2 rij = xi - particles_[j].p;
                                float r2 = rij.x*rij.x + rij.y*rij.y;
                                if (r2 >= hKer2) continue;
                                Vec2 gradW = gradW_spiky(rij) / rho0;
                                // Artificial pressure term (s_corr)
                                float sc = 0.0f;
                                if (scorrK > 0.0f && wq > 0.0f) {
                                    float w = W_poly6(r2);
                                    sc = -scorrK * std::pow(w / wq, scorrN);
                                }
                                delta += (pbfLambda_[i] + pbfLambda_[j] + sc) * gradW;
                            }
                        }
                    }
                    pbfDelta_[i] = delta / rho0;
                }
            }
        }

        // 3) Apply corrections with cap for stability
        const float maxDisp = 0.1f * hKer; // per-iteration cap
        for (int i = 0; i < n; ++i) {
            Vec2 d = pbfDelta_[i];
            float dl = std::sqrt(d.x*d.x + d.y*d.y);
            if (dl > maxDisp) d *= (maxDisp / dl);
            if (d.x == 0.0f && d.y == 0.0f) continue;
            particles_[i].p += d;
            // Keep within domain/colliders softly during iterations
            enforceParticleCollisions(particles_[i]);
        }

        // Rebuild bins after moving particles for next iteration
        for (auto& b : bins) b.clear();
        for (int i = 0; i < n; ++i) { int b = binIndex(particles_[i].p); if (b >= 0) bins[b].push_back(i); }
    }

    // Update velocities from position changes (momentum-consistent)
    Vec2 totalMomentumBefore(0,0), totalMomentumAfter(0,0);
    const float m = 1.0f; // unit mass per particle
    for (int i = 0; i < n; ++i) totalMomentumBefore += particles_[i].v * m;
    for (int i = 0; i < n; ++i) {
        Vec2 dp = particles_[i].p - posPrev[i];
        if (dt > 0.0f) particles_[i].v += dp / dt;
    }
    for (int i = 0; i < n; ++i) totalMomentumAfter += particles_[i].v * m;
    Vec2 dP = totalMomentumAfter - totalMomentumBefore;
    metricMomentumDrift_ = std::sqrt(dP.x*dP.x + dP.y*dP.y);

    // Metrics: density error mean and max overlap proxy (based on density)
    float sumErr = 0.0f; float maxErr = 0.0f;
    for (int i = 0; i < n; ++i) { float e = std::abs(pbfDensity_[i]/rho0 - 1.0f); sumErr += e; if (e > maxErr) maxErr = e; }
    metricMeanDensityErr_ = (n > 0) ? (sumErr / n) : 0.0f;
    metricMaxOverlap_ = maxErr * h_; // scale by cell to give length-ish proxy

    // Metric: horizontal-velocity variance in top band (y>0.9) using particle velocities
    float mean = 0.0f; float var = 0.0f; int cnt = 0;
    for (int i = 0; i < n; ++i) if (particles_[i].p.y > 0.9f) { mean += particles_[i].v.x; cnt++; }
    if (cnt > 0) mean /= cnt;
    for (int i = 0; i < cnt; ++i) {}
    if (cnt > 0) {
        for (int i = 0; i < n; ++i) if (particles_[i].p.y > 0.9f) { float d = particles_[i].v.x - mean; var += d*d; }
        var /= cnt;
    }
    metricTopBandVelVar_ = var;
}
