cmake_minimum_required(VERSION 3.20)
project(RBFs LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# ── Glob your sources (re-run cmake when files are added) ──────────────────────
file(GLOB_RECURSE PROJECT_SOURCES CONFIGURE_DEPENDS
        "src/*.cpp"
)
file(GLOB_RECURSE PROJECT_HEADERS CONFIGURE_DEPENDS
        "include/*.hpp"
)

# ── ImGui + ImGui-SFML (vendored) ──────────────────────────────────────────────
set(IMGUI_DIR ${CMAKE_SOURCE_DIR}/external/imgui)
set(IMGUI_SFML_DIR ${CMAKE_SOURCE_DIR}/external/imgui-sfml)
set(IMPLOT_DIR ${CMAKE_SOURCE_DIR}/external/implot)


file(GLOB IMGUI_SOURCES ${IMGUI_DIR}/*.cpp)
file(GLOB IMPLOT_SOURCES ${IMPLOT_DIR}/*.cpp)

# ── SFML (system or set SFML_DIR to local install) ────────────────────────────
find_package(SFML 3 COMPONENTS Graphics REQUIRED)
find_package(OpenGL REQUIRED)
find_package(OpenMP REQUIRED)

#target_link_libraries(your_target PRIVATE )
# ── Executable ────────────────────────────────────────────────────────────────
add_executable(${PROJECT_NAME}
    ${PROJECT_SOURCES}
    ${PROJECT_HEADERS}
        include/function.hpp
        include/constants.hpp   # optional: makes headers visible in IDEs
        include/mesh.hpp
        include/fluid-solver.hpp
        include/RBF.hpp
        include/cells.hpp
        include/constants.hpp
        ${IMGUI_SOURCES}
        ${IMGUI_SFML_DIR}/imgui-SFML.cpp
        ${IMPLOT_SOURCES}
#        include/renderer.hpp
#        include/constants.hpp
#        include/vec.hpp
#        include/flag.hpp
#        src/renderer.cpp
#        src/main.cpp
#        include/utils.hpp
#        include/menu.hpp
#        src/menu.cpp
#        include/fluid.hpp
#        src/fluid.cpp
#        include/profiler.hpp
)

#target_compile_definitions(${PROJECT_NAME} PRIVATE
#        PARTICLE_TEXTURE_PATH="${CMAKE_SOURCE_DIR}/resources/textures/circle.png"
#)


target_include_directories(${PROJECT_NAME} PRIVATE
        include
        ${IMGUI_DIR}
        ${IMGUI_SFML_DIR}
        ${IMPLOT_DIR}
)

target_link_libraries(${PROJECT_NAME} PRIVATE
        SFML::Graphics
        SFML::Window
        SFML::System
        OpenGL::GL
)

# ── Optional: OpenGL (imgui-sfml may need it) ─────────────────────────────────
#find_package(OpenGL REQUIRED)
target_link_libraries(${PROJECT_NAME} PRIVATE
        #        OpenGL::GL
        OpenMP::OpenMP_CXX
)
