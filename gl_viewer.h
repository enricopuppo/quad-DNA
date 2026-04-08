#pragma once
#include <functional>
#include "mesh.h"


typedef unsigned int uint;

class GLFWwindow;

class GlViewer
{

public:
    GlViewer();
    static void update();
    bool openWindow();
    void waitForClose();
    void renderMesh( const Mesh& m);

    bool isRunning() const;
    void setKeyCallback(std::function<bool(int,int)> cb) { keyCallback = cb; }
    void close();
    void resetBuffers();
    void applySpinRotation(float degrees);

private:
    GLFWwindow* window = nullptr;
    void createBuffers( const Mesh& m);
    bool buffersReady = false;
    std::function<bool(int,int)> keyCallback;

    static void glfwKeyCallback(GLFWwindow* w, int key, int scancode, int action, int mode);
    static void glfwMouseButtonCallback(GLFWwindow* w, int button, int action, int mods);
    static void glfwCursorPosCallback(GLFWwindow* w, double x, double y);

    void renderBuffers();



};
