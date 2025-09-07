#include <cstdio>
#include <vector>
#include <string>
#include <chrono>

#ifdef __APPLE__
#include <OpenGL/gl.h>
#else
#include <GL/gl.h>
#endif

#include <GLFW/glfw3.h>

#include <imgui.h>
#include <imgui/backends/imgui_impl_glfw.h>
#include <imgui/backends/imgui_impl_opengl2.h>

#include "flip_solver.h"
#include "renderer.h"

enum class ToolMode {
    None,
    DrawWater,
    AddWaterRect,
    AddColliderRect,
    AddColliderCircle,
    EraseFluid
};

struct UIState {
    bool running = false;
    float dt = 1.0f/60.0f;
    float gravity = -9.8f;
    float flipRatio = 0.7f;
    int pressureIters = 60;
    int substeps = 4;
    float brushRadius = 0.04f;
    ToolMode tool = ToolMode::DrawWater;
    RenderSettings render;
    bool showDemo = false;
    // Volume preservation
    bool enableMinSep = true;
    int minSepIters = 2;
    float minSepRelax = 0.6f;
    bool enableDensity = true;
    float densityStrength = 0.15f;
    int densityBlur = 0;
};

static void glfwErrorCallback(int error, const char* description) {
    std::fprintf(stderr, "GLFW Error %d: %s\n", error, description);
}

int main() {
    glfwSetErrorCallback(glfwErrorCallback);
    if (!glfwInit()) return 1;
    GLFWwindow* window = glfwCreateWindow(1280, 720, "2D FLIP Water - FluidSolver", nullptr, nullptr);
    if (!window) { glfwTerminate(); return 1; }
    glfwMakeContextCurrent(window);
    glfwSwapInterval(1);

    // Setup ImGui
    IMGUI_CHECKVERSION();
    ImGui::CreateContext();
    ImGuiIO& io = ImGui::GetIO();
    io.IniFilename = nullptr; // avoid writing ini to disk in sandbox
    ImGui::StyleColorsDark();
    ImGui_ImplGlfw_InitForOpenGL(window, true);
    ImGui_ImplOpenGL2_Init();

    // Nice styling tweaks
    ImGuiStyle& style = ImGui::GetStyle();
    style.FrameRounding = 5.0f;
    style.GrabRounding = 4.0f;
    style.WindowRounding = 7.0f;
    style.WindowTitleAlign = ImVec2(0.5f, 0.5f);

    FlipParams params;
    FlipSolver2D solver(params);
    UIState ui;

    // Initial water block for a nice start
    solver.addFluidRect(Vec2(0.1f, 0.1f), Vec2(0.9f, 0.35f));

    bool dragging = false;
    Vec2 dragStart{0,0};

    auto lastTime = std::chrono::high_resolution_clock::now();
    while (!glfwWindowShouldClose(window)) {
        glfwPollEvents();
        int fbW, fbH; glfwGetFramebufferSize(window, &fbW, &fbH);

        // Update solver parameters live
        solver.setGravity(ui.gravity);
        solver.setFlipRatio(ui.flipRatio);
        solver.setPressureIters(ui.pressureIters);
        solver.setSubsteps(ui.substeps);
        // Volume preservation params
        solver.setMinSeparationEnabled(ui.enableMinSep);
        solver.setMinSeparationIters(ui.minSepIters);
        solver.setMinSeparationRelax(ui.minSepRelax);
        solver.setDensityRelaxEnabled(ui.enableDensity);
        solver.setDensityStrength(ui.densityStrength);
        solver.setDensityBlur(ui.densityBlur);

        // Step simulation
        auto now = std::chrono::high_resolution_clock::now();
        float realDt = std::chrono::duration<float>(now - lastTime).count();
        lastTime = now;
        if (ui.running) {
            // simulate with fixed dt for stability
            solver.step(ui.dt);
        }

        // Render scene background
        Renderer2D::drawBackground();
        Viewport vp; Renderer2D::computeViewport(fbW, fbH, vp);
        Renderer2D::drawDomain(vp);
        Renderer2D::drawColliders(solver, vp, ui.render);
        if (ui.render.showGrid) Renderer2D::drawGrid(solver, vp, ui.render);
        if (ui.render.showVel) Renderer2D::drawVelocity(solver, vp, ui.render);
        Renderer2D::drawParticles(solver, vp, ui.render);

        // Start ImGui frame
        ImGui_ImplOpenGL2_NewFrame();
        ImGui_ImplGlfw_NewFrame();
        ImGui::NewFrame();

        // Left sidebar (fixed, non-overlapping) combining all controls
        ImGui::SetNextWindowPos(ImVec2(0, 0), ImGuiCond_Always);
        int fbW_ui = fbW, fbH_ui = fbH;
        ImGui::SetNextWindowSize(ImVec2(360.0f, (float)fbH_ui), ImGuiCond_Always);
        ImGuiWindowFlags sidebarFlags = ImGuiWindowFlags_NoResize | ImGuiWindowFlags_NoMove |
                                        ImGuiWindowFlags_NoSavedSettings | ImGuiWindowFlags_NoCollapse;
        if (ImGui::Begin("Controls", nullptr, sidebarFlags)) {
            // Simulation Controls
            if (ImGui::Button(ui.running ? "Pause" : "Play")) ui.running = !ui.running;
            ImGui::SameLine();
            if (ImGui::Button("Step")) { solver.step(ui.dt); }
            ImGui::SameLine();
            if (ImGui::Button("Reset")) { solver.clearFluid(); solver.clearColliders(); }
            ImGui::Separator();
            ImGui::SliderFloat("Time Step", &ui.dt, 1.0f/240.0f, 1.0f/30.0f, "%.4f s");
            ImGui::SliderInt("Substeps (min)", &ui.substeps, 1, 8);
            ImGui::SliderFloat("Gravity Y", &ui.gravity, -25.0f, 0.0f);
            ImGui::SliderFloat("FLIP Ratio", &ui.flipRatio, 0.0f, 1.0f);
            ImGui::SliderInt("Pressure Iters", &ui.pressureIters, 5, 200);
            ImGui::Separator();
            ImGui::Text("Particles: %d", solver.particleCount());
            ImGui::Checkbox("Show Grid", &ui.render.showGrid); ImGui::SameLine();
            ImGui::Checkbox("Vel", &ui.render.showVel);
            if (ui.render.showVel) ImGui::SliderInt("Vel Stride", &ui.render.velStride, 4, 32);
            ImGui::SliderFloat("Particle Size", &ui.render.particleSize, 1.0f, 8.0f);
            ImGui::Checkbox("Show Colliders", &ui.render.showColliders);
            ImGui::Checkbox("ImGui Demo", &ui.showDemo);

            ImGui::Separator();
            if (ImGui::CollapsingHeader("Initialization & Tools", ImGuiTreeNodeFlags_DefaultOpen)) {
                int mode = (int)ui.tool;
                const char* modes[] = { "None", "Draw Water", "Water Rect", "Collider Rect", "Collider Circle", "Erase" };
                ImGui::Combo("Tool", &mode, modes, IM_ARRAYSIZE(modes));
                ui.tool = (ToolMode)mode;
                ImGui::SliderFloat("Brush Radius", &ui.brushRadius, 0.005f, 0.15f);
                if (ImGui::Button("Clear Fluid")) solver.clearFluid();
                ImGui::SameLine(); if (ImGui::Button("Clear Colliders")) solver.clearColliders();
                ImGui::Separator();
                static int newNx = solver.nx(); static int newNy = solver.ny();
                ImGui::InputInt("Grid NX", &newNx); ImGui::SameLine(); ImGui::InputInt("NY", &newNy);
                if (ImGui::Button("Resize Grid")) { solver.resize(std::max(4, newNx), std::max(4, newNy)); }
            }
            ImGui::Separator();
            if (ImGui::CollapsingHeader("Volume Preservation", ImGuiTreeNodeFlags_DefaultOpen)) {
                ImGui::Checkbox("Min Separation", &ui.enableMinSep);
                ImGui::SameLine(); ImGui::SliderInt("Iters", &ui.minSepIters, 0, 5);
                ImGui::SliderFloat("Relax", &ui.minSepRelax, 0.0f, 1.0f);
                ImGui::Separator();
                ImGui::Checkbox("Density Relax", &ui.enableDensity);
                ImGui::SliderFloat("Strength", &ui.densityStrength, 0.0f, 0.5f);
                ImGui::SliderInt("Size (blur)", &ui.densityBlur, 0, 4);
            }
        }
        ImGui::End();

        // Handle mouse interactions
        if (!ImGui::GetIO().WantCaptureMouse) {
            double mx, my; glfwGetCursorPos(window, &mx, &my);
            Vec2 mp = Vec2((float)mx, (float)my);
            Vec2 wp = Renderer2D::screenToWorld(mp, vp);
            bool inDomain = (wp.x >= 0 && wp.x <= 1 && wp.y >= 0 && wp.y <= 1);
            int left = glfwGetMouseButton(window, GLFW_MOUSE_BUTTON_LEFT);
            if (left == GLFW_PRESS && inDomain) {
                if (!dragging) { dragging = true; dragStart = wp; }
                // Continuous draw while pressed
                if (ui.tool == ToolMode::DrawWater) {
                    // Spawn water in brush circle as scattered single particles
                    float r = ui.brushRadius;
                    float spacing = solver.h() * 0.6f; // slightly denser than grid faces
                    int steps = std::max(1, (int)std::ceil(2*r/spacing));
                    for (int j = -steps; j <= steps; ++j) {
                        for (int i = -steps; i <= steps; ++i) {
                            Vec2 p = wp + Vec2(i*spacing, j*spacing);
                            if (dot(p - wp, p - wp) <= r*r) solver.addParticle(p);
                        }
                    }
                } else if (ui.tool == ToolMode::EraseFluid) {
                    solver.eraseFluidCircle(wp, ui.brushRadius);
                }
            } else if (left == GLFW_RELEASE && dragging) {
                // Commit shapes on release
                if (inDomain) {
                    if (ui.tool == ToolMode::AddWaterRect) {
                        Vec2 a = Vec2(std::min(dragStart.x, wp.x), std::min(dragStart.y, wp.y));
                        Vec2 b = Vec2(std::max(dragStart.x, wp.x), std::max(dragStart.y, wp.y));
                        solver.addFluidRect(a, b);
                    } else if (ui.tool == ToolMode::AddColliderRect) {
                        Vec2 a = Vec2(std::min(dragStart.x, wp.x), std::min(dragStart.y, wp.y));
                        Vec2 b = Vec2(std::max(dragStart.x, wp.x), std::max(dragStart.y, wp.y));
                        solver.addColliderRect(a, b);
                    } else if (ui.tool == ToolMode::AddColliderCircle) {
                        float r = length(wp - dragStart);
                        solver.addColliderCircle(dragStart, r);
                    }
                }
                dragging = false;
            }
        }

        // Shape preview overlay
        if (dragging) {
            double mx, my; glfwGetCursorPos(window, &mx, &my);
            Vec2 mp = Vec2((float)mx, (float)my); Vec2 wp = Renderer2D::screenToWorld(mp, vp);
            if (ui.tool == ToolMode::AddWaterRect || ui.tool == ToolMode::AddColliderRect) {
                Vec2 a = Vec2(std::min(dragStart.x, wp.x), std::min(dragStart.y, wp.y));
                Vec2 b = Vec2(std::max(dragStart.x, wp.x), std::max(dragStart.y, wp.y));
                Vec2 s0 = Renderer2D::worldToScreen(a, vp);
                Vec2 s1 = Renderer2D::worldToScreen(b, vp);
                glEnable(GL_BLEND); glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
                glColor4f(0.25f, 0.7f, 1.0f, 0.2f);
                glBegin(GL_QUADS);
                glVertex2f(s0.x, s0.y); glVertex2f(s1.x, s0.y); glVertex2f(s1.x, s1.y); glVertex2f(s0.x, s1.y);
                glEnd();
                glColor3f(0.25f, 0.7f, 1.0f);
                glBegin(GL_LINE_LOOP);
                glVertex2f(s0.x, s0.y); glVertex2f(s1.x, s0.y); glVertex2f(s1.x, s1.y); glVertex2f(s0.x, s1.y);
                glEnd(); glDisable(GL_BLEND);
            } else if (ui.tool == ToolMode::AddColliderCircle) {
                float r = length(wp - dragStart);
                const int seg = 48;
                glEnable(GL_BLEND); glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
                glColor4f(0.25f, 0.7f, 1.0f, 0.2f);
                glBegin(GL_TRIANGLE_FAN);
                Vec2 c = Renderer2D::worldToScreen(dragStart, vp);
                glVertex2f(c.x, c.y);
                for (int k = 0; k <= seg; ++k) {
                    float ang = (float)k / seg * 6.2831853f;
                    Vec2 wp2 = Vec2(dragStart.x + std::cos(ang)*r, dragStart.y + std::sin(ang)*r);
                    Vec2 sp = Renderer2D::worldToScreen(wp2, vp);
                    glVertex2f(sp.x, sp.y);
                }
                glEnd();
                glColor3f(0.25f, 0.7f, 1.0f);
                glBegin(GL_LINE_LOOP);
                for (int k = 0; k < seg; ++k) {
                    float ang = (float)k / seg * 6.2831853f;
                    Vec2 wp2 = Vec2(dragStart.x + std::cos(ang)*r, dragStart.y + std::sin(ang)*r);
                    Vec2 sp = Renderer2D::worldToScreen(wp2, vp);
                    glVertex2f(sp.x, sp.y);
                }
                glEnd(); glDisable(GL_BLEND);
            }
        }

        if (ui.showDemo) ImGui::ShowDemoWindow(&ui.showDemo);

        // Render ImGui
        ImGui::Render();
        ImGui_ImplOpenGL2_RenderDrawData(ImGui::GetDrawData());

        glfwSwapBuffers(window);
    }

    // Cleanup
    ImGui_ImplOpenGL2_Shutdown();
    ImGui_ImplGlfw_Shutdown();
    ImGui::DestroyContext();
    glfwDestroyWindow(window);
    glfwTerminate();
    return 0;
}
