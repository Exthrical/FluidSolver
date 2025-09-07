# FluidSolver 2D FLIP Water Simulation

A compact C++ OpenGL/GLFW app that simulates and displays a 2D FLIP water simulation with an initialization phase for drawing water and adding collision shapes. UI is built with Dear ImGui (GLFW + OpenGL2 backends).

## Features

- 2D FLIP solver (particles + MAC grid) with PIC/FLIP blending
- Incompressibility via iterative pressure projection
- Solid boundaries at the domain edges (window rectangle)
- User-defined colliders (rectangles and circles)
- Initialization tools: draw water brush, add water rectangles, add colliders
- Live controls: play/pause, step, dt, substeps, gravity, FLIP ratio, iterations
- Visual toggles: grid, velocity field, particle size, colliders

## Build

Dependencies you provide on your end:
- OpenGL (system)
- GLFW 3
- Dear ImGui headers and backends for GLFW/OpenGL2

The project expects ImGui headers to be available via an include path. You can point IMGUI_DIR to your local Dear ImGui directory (containing imgui.h and backends/).

cmake -S . -B build -DIMGUI_DIR=/path/to/imgui
cmake --build build -j

If glfw3 is not discovered automatically, make sure your environment exposes a CMake glfw/glfw3 target or pkg-config entry. On Windows, you can link your installed GLFW by adding it to your toolchain or CMake prefix path.

Run the app:

./build/FluidSolver

## Usage

- Initialization phase:
  - Use the Tools panel to choose a mode: Draw Water, Water Rect, Collider Rect, Collider Circle.
  - Click-drag in the domain to add water rectangles or collider shapes.
  - Draw Water continuously adds particles within a circular brush.
  - Clear Fluid / Clear Colliders available.
- Simulation:
  - Play/Pause to run or stop.
  - Step advances a single fixed time step.
  - Adjust dt, substeps, gravity, FLIP ratio, pressure iterations live.
- Display:
  - Toggle grid, velocity field; adjust velocity stride and particle size.

Notes:
- The domain is the large square in the window; its borders act as collision walls.
- For stability and performance, start with moderate grid sizes (e.g., 128x96) and adjust particle counts via brush usage.

## Implementation Notes

- MAC grid layout: u on vertical faces, v on horizontal faces; pressure cell-centered.
- Pressure solve: simple Gauss-Seidel iterations with solid boundary handling.
- PIC/FLIP blending: particle velocity updated by FLIP delta and PIC velocity with a tunable ratio.
- Colliders: axis-aligned rectangles and circles with signed distance projection for particle collision; grid velocities zeroed at solid faces.

Enjoy experimenting!

