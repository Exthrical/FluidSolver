#include "flip_solver.h"
#include <cassert>

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
        // Resolve excessive clustering to preserve volume thickness
        if (params_.enableMinSeparation) enforceMinimumSeparation(sdt);
        clearGrid();
        buildMarkers();
        particlesToGrid();
        applySolidBoundariesOnGrid();
        uOld_ = u_; vOld_ = v_;
        computeDivergence();
        pressureSolve();
        subtractPressureGradient();
        gridToParticles();
        if (params_.enableDensityRelax) enforceCellDensity(sdt);
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

void FlipSolver2D::enforceMinimumSeparation(float dt) {
    if (particles_.empty()) return;
    // Target spacing based on desired particles-per-cell
    float spacing = h_ / std::sqrt(std::max(1, params_.particlesPerCell));
    float rSep = 0.9f * spacing; // minimum allowed separation
    float r2 = rSep * rSep;
    const int iterations = std::max(0, params_.minSepIterations);
    const float relax = params_.minSepRelax; // move fraction of computed displacement per iteration

    // Build cell bins for neighbor search
    auto buildBins = [&](std::vector<std::vector<int>>& bins){
        bins.clear();
        bins.resize(nx_ * ny_);
        for (int i = 0; i < (int)particles_.size(); ++i) {
            const Vec2& p = particles_[i].p;
            int ix = (int)std::floor(p.x / hx_);
            int iy = (int)std::floor(p.y / hy_);
            if (ix < 0 || iy < 0 || ix >= nx_ || iy >= ny_) continue;
            bins[idxC(ix, iy)].push_back(i);
        }
    };

    std::vector<Vec2> disp(particles_.size(), Vec2(0,0));
    std::vector<std::vector<int>> bins;

    for (int it = 0; it < iterations; ++it) {
        std::fill(disp.begin(), disp.end(), Vec2(0,0));
        buildBins(bins);
        for (int jy = 0; jy < ny_; ++jy) {
            for (int ix = 0; ix < nx_; ++ix) {
                int c = idxC(ix, jy);
                const auto& cell = bins[c];
                if (cell.empty()) continue;
                for (int dj = -1; dj <= 1; ++dj) {
                    int nyi = jy + dj; if (nyi < 0 || nyi >= ny_) continue;
                    for (int di = -1; di <= 1; ++di) {
                        int nxi = ix + di; if (nxi < 0 || nxi >= nx_) continue;
                        const auto& nb = bins[idxC(nxi, nyi)];
                        for (int aIdx = 0; aIdx < (int)cell.size(); ++aIdx) {
                            int a = cell[aIdx];
                            // To avoid double counting, only pair with b where (nb cell > cell) or (same cell and b index > a index)
                            for (int bIdx = 0; bIdx < (int)nb.size(); ++bIdx) {
                                int b = nb[bIdx];
                                if (&nb == &cell) { if (bIdx <= aIdx) continue; }
                                else {
                                    // If neighbor cell index is lower in linear order, skip to avoid double counting
                                    int thisLin = idxC(ix, jy);
                                    int nbLin = idxC(nxi, nyi);
                                    if (nbLin < thisLin) continue;
                                }
                                Vec2 d = particles_[a].p - particles_[b].p;
                                float d2 = d.x*d.x + d.y*d.y;
                                if (d2 >= r2 || d2 <= 1e-12f) continue;
                                float dist = std::sqrt(d2);
                                float overlap = rSep - dist;
                                Vec2 n = d / dist;
                                Vec2 corr = n * (0.5f * overlap);
                                disp[a] += corr;
                                disp[b] -= corr;
                            }
                        }
                    }
                }
            }
        }
        // Apply displacements with relaxation and update velocities accordingly
        for (size_t i = 0; i < particles_.size(); ++i) {
            if (disp[i].x == 0.0f && disp[i].y == 0.0f) continue;
            Vec2 delta = disp[i] * relax;
            particles_[i].p += delta;
            // Update velocity to be consistent with position change
            if (dt > 0.0f) particles_[i].v += delta / dt;
            enforceParticleCollisions(particles_[i]);
        }
    }
}

void FlipSolver2D::enforceCellDensity(float dt) {
    if (particles_.empty()) return;
    // Target per-cell occupancy
    const float n0 = std::max(1, params_.particlesPerCell);
    std::vector<float> occ(nx_ * ny_, 0.0f);
    for (const auto& p : particles_) {
        int ix = (int)std::floor(p.p.x / hx_);
        int iy = (int)std::floor(p.p.y / hy_);
        if (ix < 0 || iy < 0 || ix >= nx_ || iy >= ny_) continue;
        occ[idxC(ix, iy)] += 1.0f;
    }
    // Optional blur to control size of density forces
    int br = std::max(0, params_.densityBlur);
    if (br > 0) {
        std::vector<float> tmp = occ;
        // box blur in x
        for (int j = 0; j < ny_; ++j) {
            for (int i = 0; i < nx_; ++i) {
                float sum = 0.0f; int cnt = 0;
                for (int di = -br; di <= br; ++di) {
                    int ii = i + di; if (ii < 0 || ii >= nx_) continue; sum += tmp[idxC(ii, j)]; cnt++;
                }
                occ[idxC(i,j)] = (cnt > 0) ? sum / cnt : tmp[idxC(i,j)];
            }
        }
        tmp = occ;
        // box blur in y
        for (int j = 0; j < ny_; ++j) {
            for (int i = 0; i < nx_; ++i) {
                float sum = 0.0f; int cnt = 0;
                for (int dj = -br; dj <= br; ++dj) {
                    int jj = j + dj; if (jj < 0 || jj >= ny_) continue; sum += tmp[idxC(i, jj)]; cnt++;
                }
                occ[idxC(i,j)] = (cnt > 0) ? sum / cnt : tmp[idxC(i,j)];
            }
        }
    }
    // Compute simple gradients of occupancy to steer particles toward lower density
    std::vector<Vec2> grad(nx_ * ny_, Vec2(0,0));
    for (int j = 0; j < ny_; ++j) {
        for (int i = 0; i < nx_; ++i) {
            float l = (i > 0) ? occ[idxC(i-1, j)] : occ[idxC(i, j)];
            float r = (i < nx_-1) ? occ[idxC(i+1, j)] : occ[idxC(i, j)];
            float b = (j > 0) ? occ[idxC(i, j-1)] : occ[idxC(i, j)];
            float t = (j < ny_-1) ? occ[idxC(i, j+1)] : occ[idxC(i, j)];
            grad[idxC(i,j)] = Vec2(0.5f * (r - l), 0.5f * (t - b));
        }
    }
    const float k = params_.densityStrength; // strength per application
    for (auto& p : particles_) {
        int ix = (int)std::floor(p.p.x / hx_);
        int iy = (int)std::floor(p.p.y / hy_);
        if (ix < 0 || iy < 0 || ix >= nx_ || iy >= ny_) continue;
        int c = idxC(ix, iy);
        float oc = occ[c] - n0;
        if (oc <= 0.0f) continue;
        Vec2 g = grad[c];
        float gl = std::sqrt(g.x*g.x + g.y*g.y);
        if (gl < 1e-6f) {
            // No gradient info: bias a tiny push up to avoid infinite flattening
            g = Vec2(0.0f, 1.0f); gl = 1.0f;
        }
        Vec2 dir = g / gl; // toward higher density, so move opposite
        Vec2 delta = dir * (-k * (oc / n0) * h_);
        p.p += delta;
        if (dt > 0.0f) p.v += delta / dt;
        enforceParticleCollisions(p);
    }
}
