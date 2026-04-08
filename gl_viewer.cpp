#include <glad/glad.h>
#include <GLFW/glfw3.h>
#include <iostream>
#include <glm/mat4x4.hpp>
#define GLM_ENABLE_EXPERIMENTAL
#include <glm/gtx/quaternion.hpp>
#include <glm/ext/matrix_transform.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp>

#include "gl_viewer.h"
#include "mesh.h"

class Trackball
{
public:
    Trackball() = default;

    void mouseDown(float x, float y, int fbWidth, int fbHeight);
    bool mouseMove(float x, float y, int fbWidth, int fbHeight);
    void mouseUp();

    glm::mat4 rotationMatrix() const;
    void applyRotation(float angleDegrees, const glm::vec3 &axis);

private:
    glm::quat currentRotation = glm::quat(1.0f, 0.0f, 0.0f, 0.0f);
    glm::quat dragRotation = glm::quat(1.0f, 0.0f, 0.0f, 0.0f);

    glm::vec3 lastPoint{};
    bool dragging = false;

    static glm::vec3 projectOnSphere(float x, float y, int fbWidth, int fbHeight);
};

// too lazy to make them fields:
static GLuint vao, vbo, eboFill, eboWire, eboFillN, eboWireN, shaderProgramFill, shaderProgramWire;
static int fbWidth, fbHeight;
static Trackball trackball;
static bool upToDate = false;

static BoundingCube boundingCube;

void GlViewer::glfwKeyCallback(GLFWwindow *w, int key, int scancode, int action, int mode)
{
    if (action != GLFW_PRESS)
        return;
    if (key == GLFW_KEY_ESCAPE)
        glfwTerminate();

    GlViewer *self = static_cast<GlViewer *>(glfwGetWindowUserPointer(w));

    if (self && self->keyCallback)
    {
        if (self->keyCallback(key, action))
            update();
    }
}

void GlViewer::resetBuffers()
{
    glDeleteBuffers(1, &eboFill);
    glDeleteBuffers(1, &eboWire);
    glDeleteBuffers(1, &vbo);
    glDeleteVertexArrays(1, &vao);
    buffersReady = false;
}

void GlViewer::glfwMouseButtonCallback(GLFWwindow *w, int button, int action, int mods)
{

    auto *self = static_cast<GlViewer *>(glfwGetWindowUserPointer(w));
    if (button != GLFW_MOUSE_BUTTON_LEFT)
        return;

    double x, y;
    glfwGetCursorPos(w, &x, &y);

    if (action == GLFW_PRESS)
        trackball.mouseDown((float)x, (float)y, fbWidth, fbHeight);
    else if (action == GLFW_RELEASE)
        trackball.mouseUp();
}

void GlViewer::glfwCursorPosCallback(GLFWwindow *w, double x, double y)
{
    // auto* self = static_cast<GlViewer*>(glfwGetWindowUserPointer(w));
    if (trackball.mouseMove((float)x, (float)y, fbWidth, fbHeight))
        update();
}

void GlViewer::applySpinRotation(float degrees)
{
    trackball.applyRotation(degrees, glm::vec3(0.0f, 1.0f, 0.0f));
}

static const char *vertexShaderSrc = R"glsl(
#version 330 core
layout(location = 0) in vec3 inPos;
out vec3 worldPos;
uniform vec3 center;
uniform float scale;
uniform mat4 modelMatrix;

void main() {
    //gl_Position = MVP * vec4(inPos, 1.0);

    worldPos = inPos;
    vec3 scale = vec3(scale,scale,-scale);
    gl_Position = (modelMatrix * vec4( (inPos-center) * scale ,1 ) );

    //gl_Position.w = 1.0;
}
)glsl";

static const char *fragmentShaderSrc = R"glsl(
#version 330 core
in vec3 worldPos;
out vec4 FragColor;
uniform mat4 modelMatrix;
void main() {
    vec3 dx = dFdx(worldPos);
    vec3 dy = dFdy(worldPos);
    vec3 normal = normalize(cross(dx, dy));
    vec3 lightDir = vec3( modelMatrix[0].z, modelMatrix[1].z, modelMatrix[2].z );
    float lambert = abs( dot( lightDir, normal ) );
    FragColor.rgb = vec3(lambert * 0.7 + 0.3 );
    FragColor.a = 0.0;
}
)glsl";

static const char *fragmentShaderWireSrc = R"glsl(
#version 330 core
out vec4 FragColor;
void main() {
    FragColor.rgb = vec3(0.0);
    FragColor.a = 0.8;
}
)glsl";

void centerView(GLuint shaderId)
{
    GLint locC = glGetUniformLocation(shaderId, "center");
    GLint locR = glGetUniformLocation(shaderId, "scale");
    glUniform3fv(locC, 1, &boundingCube.center.x);
    glUniform1f(locR, 0.98f / boundingCube.radius);
}

GlViewer::GlViewer() {}

glm::vec3 Trackball::projectOnSphere(const float x, const float y, const int fbWidth, const int fbHeight)
{
    // from  pixel to NDC [-1, 1]
    glm::vec3 p(
        (2.0f * x / static_cast<float>(fbWidth)) - 1.0f,
        -(2.0f * y / static_cast<float>(fbHeight)) + 1.0f,
        0.0f);

    float distanceFromOrigin = glm::length(glm::vec2(p));
    distanceFromOrigin = std::min(distanceFromOrigin, 1.0f);          // clamp
    p.z = -std::sqrt(1.0f - distanceFromOrigin * distanceFromOrigin); // project on front
    return glm::normalize(p);
}

void Trackball::mouseDown(const float x, const float y, const int fbWidth, const int fbHeight)
{
    lastPoint = projectOnSphere(x, y, fbWidth, fbHeight); // drag starting point
    dragRotation = glm::quat(1.0f, 0.0f, 0.0f, 0.0f);
    dragging = true;
}

bool Trackball::mouseMove(const float x, const float y, const int fbWidth, const int fbHeight)
{
    if (!dragging)
    {
        return false;
    }

    const glm::vec3 currentPoint = projectOnSphere(x, y, fbWidth, fbHeight);

    const glm::vec3 rotationAxis = glm::cross(lastPoint, currentPoint);

    const float angle = 1.0 * std::asin(std::min(glm::length(rotationAxis), 1.0f));

    if (glm::length(rotationAxis) > 1e-6f)
    {
        dragRotation = glm::angleAxis(angle, glm::normalize(rotationAxis));
        currentRotation = dragRotation * currentRotation;
    }
    lastPoint = currentPoint;
    return true;
}

void Trackball::mouseUp()
{
    dragging = false;
}

glm::mat4 Trackball::rotationMatrix() const
{
    return glm::toMat4(currentRotation);
}

void Trackball::applyRotation(const float angleDegrees, const glm::vec3 &axis)
{
    const glm::quat r = glm::angleAxis(glm::radians(angleDegrees), axis);
    currentRotation = r * currentRotation;
}

void GlViewer::createBuffers(const Mesh &m)
{
    glGenVertexArrays(1, &vao);
    glGenBuffers(1, &vbo);
    glGenBuffers(1, &eboFill);
    glGenBuffers(1, &eboWire);

    glBindVertexArray(vao);

    // Vertex buffer (dynamic positions)
    glBindBuffer(GL_ARRAY_BUFFER, vbo);
    glBufferData(GL_ARRAY_BUFFER,
                 m.verts.size() * sizeof(Vert),
                 nullptr,
                 GL_DYNAMIC_DRAW);

    // Index buffer (static)
    {
        std::vector<uint32_t> indices;
        indices.reserve(m.quads.size() * 6);

        for (const Quad &q : m.quads)
        {
            indices.push_back(q.vi[0]);
            indices.push_back(q.vi[1]);
            indices.push_back(q.vi[2]);
            indices.push_back(q.vi[0]);
            indices.push_back(q.vi[2]);
            indices.push_back(q.vi[3]);
        }
        glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, eboFill);
        glBufferData(GL_ELEMENT_ARRAY_BUFFER,
                     indices.size() * sizeof(uint32_t),
                     indices.data(),
                     GL_STATIC_DRAW);

        eboFillN = m.quads.size() * 6;
    }

    {
        std::vector<uint32_t> indices;
        indices.reserve(m.edges.size() * 2);

        for (const Edge &e : m.edges)
        {
            indices.push_back(m.viOf(e));
            indices.push_back(m.vjOf(e));
        }
        glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, eboWire);
        glBufferData(GL_ELEMENT_ARRAY_BUFFER,
                     indices.size() * sizeof(uint32_t),
                     indices.data(),
                     GL_STATIC_DRAW);

        eboWireN = m.edges.size() * 2;
    }

    // Attribute: position only
    glEnableVertexAttribArray(0);
    glVertexAttribPointer(
        0,
        3, GL_FLOAT,
        GL_FALSE,
        sizeof(Vert),
        (void *)0);

    glBindVertexArray(0);

    boundingCube = m.computeBoundingCube();

    glBindBuffer(GL_ARRAY_BUFFER, vbo);
    glBufferSubData(GL_ARRAY_BUFFER, 0, m.verts.size() * sizeof(Vert), &m.verts[0].pos);

    buffersReady = true;
}

void GlViewer::renderMesh(const Mesh &m)
{
    if (!window)
        return;
    if (glfwWindowShouldClose(window))
        return;
    if (!buffersReady)
        createBuffers(m);

    update();
    renderBuffers();
}

void GlViewer::renderBuffers()
{
    if (!buffersReady)
        return;
    if (upToDate)
        return;

    // std::cout << "RENDER" << std::endl;

    glClearColor(0.85, 0.85, 0.85, 1);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    glBindBuffer(GL_ARRAY_BUFFER, vbo);

    glEnable(GL_DEPTH_TEST);

    glBindVertexArray(vao);

    glEnable(GL_POLYGON_OFFSET_FILL);
    glPolygonOffset(+1, +1);

    glm::mat4 modelMatrix = trackball.rotationMatrix();

    glUseProgram(shaderProgramFill);
    centerView(shaderProgramFill);
    glUniformMatrix4fv(glGetUniformLocation(shaderProgramFill, "modelMatrix"), 1, GL_FALSE, glm::value_ptr(modelMatrix));

    glDisable(GL_BLEND);
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, eboFill);
    glDrawElements(GL_TRIANGLES, eboFillN, GL_UNSIGNED_INT, nullptr);

    // Wireframe
    glDisable(GL_POLYGON_OFFSET_FILL);
    glEnable(GL_BLEND);
    glBlendFunc(GL_ONE_MINUS_SRC_ALPHA, GL_SRC_ALPHA);
    glUseProgram(shaderProgramWire);
    centerView(shaderProgramWire);
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, eboWire);
    glUniformMatrix4fv(glGetUniformLocation(shaderProgramWire, "modelMatrix"), 1, GL_FALSE, glm::value_ptr(modelMatrix));
    glDrawElements(GL_LINES, eboWireN, GL_UNSIGNED_INT, nullptr);

    glBindVertexArray(0);

    glfwSwapBuffers(window);

    upToDate = true;
}

static GLuint compileShader(GLenum type, const char *source)
{
    GLuint shader = glCreateShader(type);
    glShaderSource(shader, 1, &source, nullptr);
    glCompileShader(shader);

    // Check compile status
    GLint success;
    glGetShaderiv(shader, GL_COMPILE_STATUS, &success);
    if (!success)
    {
        char infoLog[512];
        glGetShaderInfoLog(shader, 512, nullptr, infoLog);
        throw std::runtime_error(
            std::string("Shader compilation failed: ") + infoLog);
    }
    return shader;
}

static GLuint createShaderProgram(const char *vertexSrc, const char *fragmentSrc)
{
    GLuint vertexShader = compileShader(GL_VERTEX_SHADER, vertexSrc);
    GLuint fragmentShader = compileShader(GL_FRAGMENT_SHADER, fragmentSrc);

    GLuint program = glCreateProgram();
    glAttachShader(program, vertexShader);
    glAttachShader(program, fragmentShader);
    glLinkProgram(program);

    // Check link status
    GLint success;
    glGetProgramiv(program, GL_LINK_STATUS, &success);
    if (!success)
    {
        char infoLog[512];
        glGetProgramInfoLog(program, 512, nullptr, infoLog);
        throw std::runtime_error(
            std::string("Shader linking failed: ") + infoLog);
    }

    // Shaders can be deleted after linking
    glDeleteShader(vertexShader);
    glDeleteShader(fragmentShader);

    return program;
}

bool GlViewer::openWindow()
{

    if (!glfwInit())
        return false;

    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 4);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 1);
    glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
    glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE);

    window = glfwCreateWindow(640, 640, "Folding", NULL, NULL);
    if (!window)
    {
        glfwTerminate();
        return false;
    }
    glfwMakeContextCurrent(window);

    if (!gladLoadGLLoader((GLADloadproc)glfwGetProcAddress))
    {
        throw std::runtime_error("Failed to initialize GLAD");
    }

    glViewport(0, 0, 640, 640);

    shaderProgramFill = createShaderProgram(vertexShaderSrc, fragmentShaderSrc);
    shaderProgramWire = createShaderProgram(vertexShaderSrc, fragmentShaderWireSrc);

    glfwSetWindowUserPointer(window, this);
    glfwSetKeyCallback(window, glfwKeyCallback);
    glfwSetMouseButtonCallback(window, glfwMouseButtonCallback);
    glfwSetCursorPosCallback(window, glfwCursorPosCallback);
    glfwSetFramebufferSizeCallback(window, [](GLFWwindow *w, int width, int height)
                                   {
        glViewport(0, 0, width, height);
        GlViewer* self = static_cast<GlViewer*>(glfwGetWindowUserPointer(w));
        fbWidth  = width;
        fbHeight = height; });

    glfwGetFramebufferSize(window, &fbWidth, &fbHeight);
    glViewport(0, 0, fbWidth, fbHeight);

    return true;
}

void GlViewer::update()
{
    upToDate = false;
}

void GlViewer::waitForClose()
{
    if (!window)
        return;

    while (!glfwWindowShouldClose(window))
    {
        glfwWaitEvents();
        renderBuffers();
    }
    glfwTerminate();
}
