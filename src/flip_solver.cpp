#include "flip_solver.h"
#include <cassert>

FlipSolver2D::FlipSolver2D(const FlipParams& params) : params_(params) {
    resize(params_.nx, params_.ny);
}

void FlipSolver2D::resize(int nx, int ny) {
    nx_ = std::max(4, nx);
    ny_ = std::max(4, ny);
    h_ = 1.0f / std::max(nx_, ny_);
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
        advectParticles(sdt);
        clearGrid();
        buildMarkers();
        particlesToGrid();
        applySolidBoundariesOnGrid();
        uOld_ = u_; vOld_ = v_;
        computeDivergence();
        pressureSolve();
        subtractPressureGradient();
        gridToParticles();
    }
}

void FlipSolver2D::buildMarkers() {
    std::fill(marker_.begin(), marker_.end(), 0);
    // Mark solid cells from colliders (domain boundaries handled via velocity BCs)
    for (int j = 0; j < ny_; ++j) {
        for (int i = 0; i < nx_; ++i) {
            Vec2 c = Vec2((i + 0.5f) * h_, (j + 0.5f) * h_);
            if (colliders_.sdf(c) < 0.0f) marker_[idxC(i,j)] = 2;
        }
    }
    // Mark fluid cells where particles reside
    for (const auto& p : particles_) {
        int i = (int)std::floor(p.p.x / h_);
        int j = (int)std::floor(p.p.y / h_);
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
        // U faces at (i*h, (j+0.5)*h)
        float x = P.x / h_; float y = P.y / h_ - 0.5f;
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
        // V faces at ((i+0.5)*h, j*h)
        x = P.x / h_ - 0.5f; y = P.y / h_;
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
            Vec2 fc = Vec2(i * h_, (j + 0.5f) * h_);
            if (colliders_.sdf(fc) < 0.0f) solid = true;
            if (solid) u_[idx] = 0.0f;
        }
    }
    for (int j = 0; j <= ny_; ++j) {
        for (int i = 0; i < nx_; ++i) {
            int idx = idxV(i,j);
            bool solid = (j == 0 || j == ny_);
            Vec2 fc = Vec2((i + 0.5f) * h_, j * h_);
            if (colliders_.sdf(fc) < 0.0f) solid = true;
            if (solid) v_[idx] = 0.0f;
        }
    }
}

void FlipSolver2D::computeDivergence() {
    const float invh = 1.0f / h_;
    for (int j = 0; j < ny_; ++j) {
        for (int i = 0; i < nx_; ++i) {
            int c = idxC(i,j);
            if (marker_[c] == 2) { div_[c] = 0.0f; continue; }
            float uR = u_[idxU(i+1, j)];
            float uL = u_[idxU(i,   j)];
            float vT = v_[idxV(i,   j+1)];
            float vB = v_[idxV(i,   j  )];
            div_[c] = invh * (uR - uL + vT - vB);
        }
    }
}

void FlipSolver2D::pressureSolve() {
    // Gauss-Seidel iterations
    // Discrete Poisson: (sum_neighbors - N * p_c) / h^2 = div
    // => p_c = (sum_neighbors - div * h^2) / N
    for (int it = 0; it < params_.pressureIters; ++it) {
        for (int j = 0; j < ny_; ++j) {
            for (int i = 0; i < nx_; ++i) {
                int c = idxC(i,j);
                if (marker_[c] == 2 || marker_[c] == 0) { p_[c] = 0.0f; continue; }
                float sum = 0.0f; int count = 0;
                // Left
                if (i > 0) {
                    uint8_t m = marker_[idxC(i-1,j)];
                    if (m == 1) { sum += p_[idxC(i-1,j)]; count++; }
                    else if (m == 0) { count++; } // air: Dirichlet p=0
                    else { count++; } // solid: Neumann
                } else { count++; } // domain boundary: Neumann
                // Right
                if (i < nx_-1) {
                    uint8_t m = marker_[idxC(i+1,j)];
                    if (m == 1) { sum += p_[idxC(i+1,j)]; count++; }
                    else if (m == 0) { count++; }
                    else { count++; }
                } else { count++; }
                // Bottom
                if (j > 0) {
                    uint8_t m = marker_[idxC(i,j-1)];
                    if (m == 1) { sum += p_[idxC(i,j-1)]; count++; }
                    else if (m == 0) { count++; }
                    else { count++; }
                } else { count++; }
                // Top
                if (j < ny_-1) {
                    uint8_t m = marker_[idxC(i,j+1)];
                    if (m == 1) { sum += p_[idxC(i,j+1)]; count++; }
                    else if (m == 0) { count++; }
                    else { count++; }
                } else { count++; }
                if (count > 0) {
                    p_[c] = (sum - div_[c] * (h_ * h_)) / (float)count;
                } else {
                    p_[c] = 0.0f;
                }
            }
        }
    }
}

void FlipSolver2D::subtractPressureGradient() {
    const float invh = 1.0f / h_;
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
            u_[idxU(i,j)] -= (pR - pL) * invh;
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
            v_[idxV(i,j)] -= (pT - pB) * invh;
        }
    }
    applySolidBoundariesOnGrid();
}

Vec2 FlipSolver2D::sampleGridVelocity(const Vec2& P, const std::vector<float>& u, const std::vector<float>& v) const {
    // Sample u
    float x = P.x / h_; float y = P.y / h_ - 0.5f;
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
    x = P.x / h_ - 0.5f; y = P.y / h_;
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
        enforceParticleCollisions(p);
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
