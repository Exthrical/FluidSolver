// Simple deterministic test harness for the 2D FLIP + PBF solver
#include <cstdio>
#include <vector>
#include <cmath>
#include "../src/flip_solver.h"

static void run_steps(FlipSolver2D& solver, float dt, int steps) {
    for (int i = 0; i < steps; ++i) solver.step(dt);
}

static float collider_pile_test_overlap(const FlipSolver2D& solver) {
    // Approximate max overlap using neighbor search with a target spacing
    const auto& parts = solver.particles();
    int nx = solver.nx(), ny = solver.ny();
    float hx = solver.hx(), hy = solver.hy();
    float h = solver.h();
    float spacing = h / std::sqrt( (float)std::max(1, 8) ); // assume ~8 ppc typical
    float rSep = spacing; float r2 = rSep * rSep;
    std::vector<std::vector<int>> bins(nx*ny);
    auto idx = [&](int i,int j){return i + nx*j;};
    auto binOf = [&](const Vec2& p){int ix=(int)std::floor(p.x/hx);int iy=(int)std::floor(p.y/hy);if(ix<0||iy<0||ix>=nx||iy>=ny)return -1;return idx(ix,iy);};
    for (int i = 0; i < (int)parts.size(); ++i) { int b=binOf(parts[i].p); if (b>=0) bins[b].push_back(i); }
    int rx = std::max(1, (int)std::ceil(rSep / hx));
    int ry = std::max(1, (int)std::ceil(rSep / hy));
    float maxOverlap = 0.0f;
    for (int j = 0; j < ny; ++j) for (int i = 0; i < nx; ++i) {
        const auto& cell = bins[idx(i,j)]; if (cell.empty()) continue;
        for (int aIdx = 0; aIdx < (int)cell.size(); ++aIdx) {
            int a = cell[aIdx];
            for (int dj = -ry; dj <= ry; ++dj) {int jj=j+dj; if(jj<0||jj>=ny) continue; for (int di = -rx; di <= rx; ++di){int ii=i+di; if(ii<0||ii>=nx) continue; const auto& nb=bins[idx(ii,jj)]; for (int bIdx = 0; bIdx < (int)nb.size(); ++bIdx){int b=nb[bIdx]; if (b<=a && &nb==&cell) continue; Vec2 d = parts[a].p - parts[b].p; float d2 = d.x*d.x + d.y*d.y; if (d2 < r2) { float dist = std::sqrt(d2); maxOverlap = std::max(maxOverlap, rSep - dist); }}}}
        }
    }
    return maxOverlap / h; // in units of h
}

int main() {
    // Common params
    FlipParams params;
    params.nx = 128; params.ny = 96; // will override per test
    params.substeps = 2;
    params.gravity = -9.8f;
    params.flipRatio = 0.7f;
    params.minSepIterations = 4; // more robust for tests
    params.minSepRelax = 0.0f;   // stiff by default
    params.densityStrength = 0.1f;
    params.densityBlur = 2;      // kernel radius in cells

    // Test 1: NX=400,NY=96, water at top and bottom
    {
        params.nx = 400; params.ny = 96;
        FlipSolver2D solver(params);
        solver.addFluidRect(Vec2(0.05f, 0.80f), Vec2(0.95f, 0.90f));
        solver.addFluidRect(Vec2(0.05f, 0.05f), Vec2(0.95f, 0.15f));
        const float dt = 1.0f/60.0f; run_steps(solver, dt, (int)(2.0f/dt));
        float var = solver.metricTopBandVelVar();
        std::printf("Test1 NX=400 NY=96 top-band u variance: %.6e -> %s\n", var, (var > 1e-4f ? "PASS" : "FAIL"));
    }

    // Test 2: Collider pile test
    {
        params.nx = 128; params.ny = 96;
        FlipSolver2D solver(params);
        // Flat floor collider
        solver.addColliderRect(Vec2(0.0f, 0.0f), Vec2(1.0f, 0.05f));
        // Drop a blob
        solver.addFluidRect(Vec2(0.25f, 0.7f), Vec2(0.75f, 0.9f));
        const float dt = 1.0f/60.0f; run_steps(solver, dt, (int)(3.0f/dt));
        float overlapH = collider_pile_test_overlap(solver);
        std::printf("Test2 Collider pile max overlap: %.3f h -> %s\n", overlapH, (overlapH < 0.2f ? "PASS" : "FAIL"));
    }

    // Test 3 & 4: Mass and momentum conservation (proxy)
    {
        params.nx = 200; params.ny = 120;
        FlipSolver2D solver(params);
        solver.addFluidRect(Vec2(0.1f, 0.1f), Vec2(0.9f, 0.3f));
        int n0 = solver.particleCount();
        const float dt = 1.0f/60.0f; run_steps(solver, dt, (int)(2.0f/dt));
        int n1 = solver.particleCount();
        float densErr = solver.metricMeanDensityErr();
        float momDrift = solver.metricMomentumDrift();
        std::printf("Test3 Mass particles %d -> %d %s | mean density err=%.3e | momentum drift=%.3e\n",
            n0, n1, (n0==n1?"PASS":"FAIL"), densErr, momDrift);
    }

    return 0;
}

