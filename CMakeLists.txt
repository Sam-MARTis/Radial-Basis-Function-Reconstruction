cmake_minimum_required(VERSION 3.20)
project(RBFs LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
enable_testing()

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
find_package(Eigen3 REQUIRED)


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
        include/matrix.hpp
#        src/RBF.cpp
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
        Eigen3::Eigen
)

# ── Optional: OpenGL (imgui-sfml may need it) ─────────────────────────────────
#find_package(OpenGL REQUIRED)
target_link_libraries(${PROJECT_NAME} PRIVATE
        #        OpenGL::GL
        OpenMP::OpenMP_CXX
)

# ── Headless kernel checks and comparison study ───────────────────────────────
option(RBF_BUILD_TOOLS "Build the headless RBF kernel checks and study" ON)
if (RBF_BUILD_TOOLS)
    add_executable(kernel_check tools/kernel_check.cpp)
    target_include_directories(kernel_check PRIVATE include)
    target_link_libraries(kernel_check PRIVATE Eigen3::Eigen)
    add_test(NAME kernel_check COMMAND kernel_check)

    add_executable(rbf_kernel_study tools/rbf_kernel_study.cpp src/mesh.cpp)
    target_include_directories(rbf_kernel_study PRIVATE include)
    target_link_libraries(rbf_kernel_study PRIVATE Eigen3::Eigen OpenMP::OpenMP_CXX)
    # The study runs many cells/steps; do not let a Debug configuration make it unusable.
    target_compile_options(rbf_kernel_study PRIVATE -O3)
    # The study must survive an inadmissible step so it can report divergence instead
    # of aborting inside the solver; runSod() checks rho/p positivity after every step.
    target_compile_definitions(rbf_kernel_study PRIVATE NDEBUG)

    # Same study with the limiter bracket reduced to the immediate stencil (self + one
    # hop each way). This is the A/B used in rbfKernelStudy.md.
    add_executable(rbf_kernel_study_lim1 tools/rbf_kernel_study.cpp src/mesh.cpp)
    target_include_directories(rbf_kernel_study_lim1 PRIVATE include)
    target_link_libraries(rbf_kernel_study_lim1 PRIVATE Eigen3::Eigen OpenMP::OpenMP_CXX)
    target_compile_options(rbf_kernel_study_lim1 PRIVATE -O3)
    target_compile_definitions(rbf_kernel_study_lim1 PRIVATE
            NDEBUG
            RBF_LIMITER_NEIGHBOURHOOD_DEPTH_OVERRIDE=1)
endif()
