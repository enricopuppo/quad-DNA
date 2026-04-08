cmake_minimum_required(VERSION 3.16)

# avoid redowloading libs in every project
# set(FETCHCONTENT_BASE_DIR "C:/cmakeLibs" CACHE PATH "")
#set(FETCHCONTENT_UPDATES_DISCONNECTED ON)

message(STATUS "FetchContent base: ${FETCHCONTENT_BASE_DIR}")

project(dna_quad LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

add_executable(dna_quad
    main.cpp
    mesh.h mesh.cpp
    obj_import_export.cpp
    mesh_motorcycle_graph.cpp
    union_find.h
    union_find.cpp
    mesh_cutgraph_opt.cpp
    primary_structure.h primary_structure.cpp
    primary_to_secondary.cpp
    RNA.h RNA.cpp
    DNA.h DNA.cpp
    gl_viewer.h gl_viewer.cpp)

include(FetchContent)
FetchContent_Declare(
    glfw
    GIT_REPOSITORY https://github.com/glfw/glfw.git
    GIT_TAG        3.3.9
)
FetchContent_Declare(
    glad
    GIT_REPOSITORY https://github.com/Dav1dde/glad.git
    GIT_TAG v0.1.36
)
FetchContent_Declare(
    glm
    GIT_REPOSITORY https://github.com/g-truc/glm.git
    GIT_TAG        1.0.1
)
FetchContent_MakeAvailable(glm)
FetchContent_MakeAvailable(glad)
FetchContent_MakeAvailable(glfw)
find_package(OpenGL REQUIRED)
target_link_libraries(dna_quad PRIVATE glad glfw glm OpenGL::GL)
