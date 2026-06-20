cmake_minimum_required(VERSION 3.16)
project(FluidSim LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# ---------------------------------------------------------------------------
# Simulation core library (portable — no OS or SFML headers)
# ---------------------------------------------------------------------------
add_library(sim_core INTERFACE)

target_include_directories(sim_core INTERFACE
    ${CMAKE_CURRENT_SOURCE_DIR}/src/sim
)

# ---------------------------------------------------------------------------
# Desktop target (SFML visualiser)
# ---------------------------------------------------------------------------
find_package(SFML 2 COMPONENTS graphics window system REQUIRED)

add_executable(desktop
    src/desktop/main.cpp
)

target_link_libraries(desktop PRIVATE
    sim_core
    sfml-graphics
    sfml-window
    sfml-system
)

target_compile_options(desktop PRIVATE
    $<$<CONFIG:Release>:-march=native -ffast-math>
)

# ---------------------------------------------------------------------------
# ESP32 target (future — not wired up yet)
# ---------------------------------------------------------------------------
# When cross-compiling for ESP32, define a target here that links sim_core
# but replaces the desktop entry point with an ESP-IDF app_main.
#
# Example sketch (uncomment and adapt when the toolchain file is ready):
#
# if(ESP32_TARGET)
#     idf_component_register(
#         SRCS "src/esp32/main.cpp"
#         INCLUDE_DIRS "src/sim"
#     )
#     target_link_libraries(${COMPONENT_LIB} PRIVATE sim_core)
# endif()
