#include <iostream>
#include <GL/glew.h>
#include <GLFW/glfw3.h>
#include <GL/glut.h>
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp>
#include "glm/mat4x4.hpp"
#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"

//#include "Matrix.h"
#include "soil.h"
#include <vector>
#include <cmath>


using namespace glm;
using namespace std;






#define PI 3.14159265 
#define EPSILON 1.0e-5

#define IS_ZERO(v) (abs(v) < EPSILON)

#define SIGN(v) (int)(((v) > EPSILON) - ((v) < -EPSILON))

#define RESOLUTION 8

#define C 2.0

#pragma comment(lib, "glfw3.lib")
#pragma comment(lib, "glew32.lib")
#pragma comment(lib, "opengl32.lib")
using namespace glm;
using namespace std;






bool turn = false;
bool makeLine = false;
bool makeBody = false;

GLFWwindow* g_window;
GLuint g_shaderProgram1;
GLint g_uMVP;
GLint g_uMV;
GLint g_uMVP1;
GLint g_uMV1;


GLfloat yaw_var = -90.0f;
GLfloat pitch_var = 0.0f;

GLfloat lastX = 400, lastY = 300;

GLuint texture1_ID;
GLuint texture2_ID;
GLint maplocation1;
GLint maplocation2;



int texture1_width;
int texture1_height;
int texture1_cent;
int texture2_width;
int texture2_height;
int texture2_cent;
int cur_g;
int cnt_g;

int width;
int height;

class Model
{
public:
    GLuint vbo;
    GLuint ibo;
    GLuint vao;
    GLsizei indexCount;
};
class PointCoords
{
public:
    float x;
    float y;
    float z;

    PointCoords(float x_coord, float y_coord)
    {
        x = x_coord;
        y = y_coord;
        z = 0;
    }
    PointCoords()
    {
        x = -1;
        y = -1;
        z = -1;
    }
};

class Point2D
{
public:

    double x, y;


    Point2D() { x = y = 0.0; };
   
    Point2D(double _x, double _y) { x = _x; y = _y; };

    
    Point2D operator +(const Point2D& p) const { return Point2D(x + p.x, y + p.y); };
    
    Point2D operator -(const Point2D& p) const { return Point2D(x - p.x, y - p.y); };
    
    Point2D operator *(double v) const { return Point2D(x * v, y * v); };

    
    void normalize()
    {
        double l = sqrt(x * x + y * y);
        if (IS_ZERO(l))
            x = y = 0.0;
        else
        {
            x /= l;
            y /= l;
        }
    };


    static Point2D absMin(const Point2D& p1, const Point2D& p2)
    {
        return Point2D(abs(p1.x) < abs(p2.x) ? p1.x : p2.x, abs(p1.y) < abs(p2.y) ? p1.y : p2.y);
    };
};
class Segment
{
public:
    
    Point2D points[4];

    
    Point2D calc(double t) const
    {
        double t2 = t * t;
        double t3 = t2 * t;
        double nt = 1.0 - t;
        double nt2 = nt * nt;
        double nt3 = nt2 * nt;
        return Point2D(nt3 * points[0].x + 3.0 * t * nt2 * points[1].x + 3.0 * t2 * nt * points[2].x + t3 * points[3].x,
            nt3 * points[0].y + 3.0 * t * nt2 * points[1].y + 3.0 * t2 * nt * points[2].y + t3 * points[3].y);
    };
};

PointCoords* pointsArray = new PointCoords[100000];
int cur_i = 0;

Model g_model1;
Model g_model2;



GLuint createShader(const GLchar* code, GLenum type)
{
    GLuint result = glCreateShader(type);

    glShaderSource(result, 1, &code, NULL);
    glCompileShader(result);

    GLint compiled;
    glGetShaderiv(result, GL_COMPILE_STATUS, &compiled);

    if (!compiled)
    {
        GLint infoLen = 0;
        glGetShaderiv(result, GL_INFO_LOG_LENGTH, &infoLen);
        if (infoLen > 0)
        {
            //char infoLog[infoLen];
            //glGetShaderInfoLog(result, infoLen, NULL, infoLog);
            cout << "Shader compilation error" << endl << " " << endl;
        }
        glDeleteShader(result);
        return 0;
    }

    return result;
}

GLuint createProgram(GLuint vsh, GLuint fsh)
{
    GLuint result = glCreateProgram();

    glAttachShader(result, vsh);
    glAttachShader(result, fsh);

    glLinkProgram(result);

    GLint linked;
    glGetProgramiv(result, GL_LINK_STATUS, &linked);

    if (!linked)
    {
        GLint infoLen = 0;
        glGetProgramiv(result, GL_INFO_LOG_LENGTH, &infoLen);
        if (infoLen > 0)
        {
            char* infoLog = (char*)alloca(infoLen);
            glGetProgramInfoLog(result, infoLen, NULL, infoLog);
            cout << "Shader program linking error" << endl << infoLog << endl;
        }
        glDeleteProgram(result);
        return 0;
    }



    return result;
}

bool createShaderProgram1()
{
    g_shaderProgram1 = 0;

    const GLchar vsh[] =
        "#version 330\n"
        ""
        "layout(location = 0) in vec2 a_position;"
        "layout(location = 1) in vec2 a_texCoord;"
        ""
        "uniform mat4 u_mvp;"
        "uniform mat4 u_mv;"




        ""
        "out vec3 v_normal;"
        "out vec3 v_pos;"
        "out vec2 v_texCoord;"

        "float f(vec2 p) {return 0;}"
        "vec3 grad(vec2 p) {return vec3(0,0,1);}"
        ""
        "void main()"
        "{"
        "    v_texCoord=vec2((a_position.xy+0.5)/2);"
        "    v_normal =transpose(inverse(mat3(u_mv)))* normalize(grad(a_position));"
        "    float z =f(a_position);"
        "    vec4 p0=vec4(a_position.x,a_position.y,z,1);"
        "    v_pos=vec3(u_mv * p0);"
        "    gl_Position = u_mvp*  p0;"
        "}"
        ;

    const GLchar fsh[] =
        "#version 330\n"
        ""
        "in vec3 v_normal;"
        "in vec3 v_pos;"
        "in vec2 v_texCoord;"
        "uniform sampler2D u_map1;"
        "uniform sampler2D u_map2;"
        ""
        "layout(location = 0) out vec4 o_color;"
        "float rs(float n1, float n2, float cosI, float cosT) {return (n1 * cosI - n2 * cosT) / (n1 * cosI + n2 * cosT);}"
        ""
        "float rp(float n1, float n2, float cosI, float cosT) {return (n2 * cosI - n1 * cosT) / (n1 * cosT + n2 * cosI);}"
        ""
        "float ts(float n1, float n2, float cosI, float cosT) {return 2 * n1 * cosI / (n1 * cosI + n2 * cosT);}"
        ""
        "float tp(float n1, float n2, float cosI, float cosT) { return 2 * n1 * cosI / (n1 * cosT + n2 * cosI);}"
        ""
        "float thinFilmReflectance(float cos0, float lambda, float thickness, float n0, float n1, float n2) {"
        "    float PI=3.1415926535;"
        "    float d10 = (n1 > n0) ? 0 : PI;"
        "    float d12 = (n1 > n2) ? 0 : PI;"
        "    float delta = d10 + d12;"
        "    float sin1 = pow(n0 / n1, 2) * (1 - pow(cos0, 2));"
        "    if (sin1 > 1) return 1.0;"
        "    float cos1 = sqrt(1 - sin1);"
        "    float sin2 = pow(n0 / n2, 2) * (1 - pow(cos0, 2));"
        "    if (sin2 > 1) return 1.0;"
        "    float cos2 = sqrt(1 - sin2);"
        "    float alpha_s = rs(n1, n0, cos1, cos0) * rs(n1, n2, cos1, cos2);"
        "    float alpha_p = rp(n1, n0, cos1, cos0) * rp(n1, n2, cos1, cos2);"
        "    float beta_s = ts(n0, n1, cos0, cos1) * ts(n1, n2, cos1, cos2);"
        "    float beta_p = tp(n0, n1, cos0, cos1) * tp(n1, n2, cos1, cos2);"
        "    float phi = (2 * PI / lambda) * (2 * n1 * thickness * cos1) + delta;"
        "    float ts = pow(beta_s, 2) / (pow(alpha_s, 2) - 2 * alpha_s * cos(phi) + 1);"
        "    float tp = pow(beta_p, 2) / (pow(alpha_p, 2) - 2 * alpha_p * cos(phi) + 1);"
        "    float beamRatio = (n2 * cos2) / (n0 * cos0);"
        "    float t = beamRatio * (ts + tp) / 2;"
        "    return 1 - t;"
        "}"


        ""
        "void main()"
        "{"
        "   float thicknessMin = 250;"
        "   float thicknessMax = 500;"
        "   float nmedium = 1;"
        "   float nfilm = 1.5;"
        "   float ninternal = 1;"

        "   vec3 n = normalize(v_normal);"
        "   vec3 E = vec3(0,0,0);"
        "   vec3 L = vec3(0,0,5);"
        "   vec3 I = normalize(v_pos-L);"
        "   float d = max(dot(n,-I),0.1);"
        "   vec3 e = normalize(E-v_pos);"
        "   vec3 h = normalize(-I+e);"

        "   float cos0 = abs(dot(L , n));"
        "   float t = (texture2D(u_map2, v_texCoord).r+texture2D(u_map2, v_texCoord).g+texture2D(u_map2, v_texCoord).b)/ 3.0;"
        "   float thick=thicknessMin*(1.0-t)+thicknessMax*t;"
        "   float red=thinFilmReflectance(cos0, 650, thick, nmedium, nfilm, ninternal);"
        "   float green=thinFilmReflectance(cos0, 510, thick, nmedium, nfilm, ninternal);"
        "   float blue=thinFilmReflectance(cos0, 475, thick, nmedium, nfilm, ninternal);"
        "   o_color=vec4(red,green,blue,1);"
        "}"
        ;

    GLuint vertexShader, fragmentShader;

    vertexShader = createShader(vsh, GL_VERTEX_SHADER);
    fragmentShader = createShader(fsh, GL_FRAGMENT_SHADER);

    g_shaderProgram1 = createProgram(vertexShader, fragmentShader);

    g_uMVP = glGetUniformLocation(g_shaderProgram1, "u_mvp");
    g_uMV = glGetUniformLocation(g_shaderProgram1, "u_mv");
    maplocation1 = glGetUniformLocation(g_shaderProgram1, "u_map1");
    maplocation2 = glGetUniformLocation(g_shaderProgram1, "u_map2");

    glDeleteShader(vertexShader);
    glDeleteShader(fragmentShader);

    return g_shaderProgram1 != 0;
}



bool createModel1()
{

    const int elem = 2;
    const int nPoints = 100;
    const int planeSize = 4 * elem * (nPoints - 1) * (nPoints - 1);
    GLfloat* vertices = new GLfloat[planeSize];
    GLuint* indices = new GLuint[6 * (nPoints - 1) * (nPoints - 1)];
    float startPosx = -1.0f;
    float startPosy = -1.0f;
    for (int j = 0; j < (nPoints - 1); j++)
    {
        int dop = 4 * elem * (nPoints - 1) * j;
        for (int i = 0; i < (nPoints - 1) * 8; i += 8)
        {
            vertices[dop + i] = startPosx; //x
            vertices[dop + i + 1] = startPosy; //y
            vertices[dop + i + 2] = startPosx;
            vertices[dop + i + 3] = startPosy + 2.0f / (nPoints - 1);
            vertices[dop + i + 4] = startPosx + 2.0f / (nPoints - 1);
            vertices[dop + i + 5] = startPosy + 2.0f / (nPoints - 1);
            vertices[dop + i + 6] = startPosx + 2.0f / (nPoints - 1);
            vertices[dop + i + 7] = startPosy;
            startPosx = startPosx + 2.0f / (nPoints - 1);
            /* for (int k = 0; k < 12; k++)
             {
                 cout << ' ' << vertices[dop+i + k];
                 if ((k+1)%3==0 &&k!=0 )
                     std::cout << endl;

             }
             cout << endl;
             cout << "-----------------------" << endl;*/
        }
        startPosy = startPosy + 2.0f / (nPoints - 1);
        startPosx = -1.0f;
    }

    int cur = 0;
    for (int i = 0; i < (nPoints - 1) * (nPoints - 1) * 6; i += 6)
    {
        indices[i] = cur;
        indices[i + 1] = cur + 1;
        indices[i + 2] = cur + 2;
        indices[i + 3] = cur + 2;
        indices[i + 4] = cur + 3;
        indices[i + 5] = cur;
        cur += 4;
    }
    /*cout << "Группикровка: " << endl;
    for (int k = 0; k < (nPoints - 1) * (nPoints - 1) * 6; k++)
    {
        cout << ' ' << indices[k];
        if ((k + 1) % 6 == 0 && k != 0)
            std::cout << endl;

    }
    cout << endl;
    cout << "-----------------------" << endl;*/


    glGenVertexArrays(1, &g_model1.vao);
    glBindVertexArray(g_model1.vao);

    glGenBuffers(1, &g_model1.vbo);
    glBindBuffer(GL_ARRAY_BUFFER, g_model1.vbo);
    glBufferData(GL_ARRAY_BUFFER, (nPoints - 1) * (nPoints - 1) * 4 * elem * sizeof(GLfloat), vertices, GL_STATIC_DRAW);

    glGenBuffers(1, &g_model1.ibo);
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, g_model1.ibo);
    glBufferData(GL_ELEMENT_ARRAY_BUFFER, (nPoints - 1) * (nPoints - 1) * 6 * sizeof(GLuint), indices, GL_STATIC_DRAW);

    g_model1.indexCount = (nPoints - 1) * (nPoints - 1) * 6;

    glEnableVertexAttribArray(0);
    glVertexAttribPointer(0, 2, GL_FLOAT, GL_FALSE, 2 * sizeof(GLfloat), (const GLvoid*)0);
    glEnableVertexAttribArray(1);
    glVertexAttribPointer(1, 2, GL_FLOAT, GL_FALSE, 3 * sizeof(GLfloat), (const GLvoid*)(0 * sizeof(GLfloat)));

    delete[] vertices;
    delete[] indices;

    return g_model1.vbo != 0 && g_model1.ibo != 0 && g_model1.vao != 0;
}



bool createTexture()
{
    unsigned char* image;

    image = stbi_load("SsW016.jpg", &texture2_width, &texture2_height, 0, 0);

    glGenTextures(1, &texture2_ID);
    glBindTexture(GL_TEXTURE_2D, texture2_ID);
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, texture2_width, texture2_height, 0, GL_RGB, GL_UNSIGNED_BYTE, image);
    delete[] image;
    //glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
    //glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);







    return true;
}

bool init()
{
    // Set initial color of color buffer to white.
    glClearColor(1.0f, 0.0f, 0.0f, 1.0f);

    glEnable(GL_DEPTH_TEST);
    //glEnable(GL_LINE_SMOOTH);



    return createShaderProgram1() && createModel1() && createTexture();
}

void reshape(GLFWwindow* window, int width, int height)
{
    glViewport(0, 0, width, height);
}




void draw()
{
    glfwGetWindowSize(g_window, &width, &height);
    // Clear color buffer.
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    glBindTexture(GL_TEXTURE_2D, texture1_ID);
    glActiveTexture(GL_TEXTURE0);
    glUniform1i(maplocation1, 0);
    glActiveTexture(GL_TEXTURE1);
    glBindTexture(GL_TEXTURE_2D, texture2_ID);
    glUniform1i(maplocation2, 1);

    glUseProgram(g_shaderProgram1);

    glBindVertexArray(g_model1.vao);







    glm::vec3 front;
    front.x = cos(glm::radians(pitch_var)) * cos(glm::radians(yaw_var));
    front.z = sin(glm::radians(pitch_var));
    front.y = cos(glm::radians(pitch_var)) * sin(glm::radians(yaw_var));


    glm::vec3 rot = glm::normalize(front);
    glm::mat4 projection(1.f);
    projection = glm::perspective(45.0f, (GLfloat)4.0f / (GLfloat)3.0f, 0.1f, 100.0f);
    glm::mat4 model(1.f);
    glm::mat4 view(1.f);
    model = glm::scale(model, vec3(0.5f, 0.5f, 0.5f));
    model = glm::rotate(model, 3.0f, rot);
    view = glm::translate(view, glm::vec3(0.0f, -0.0f, -1.f));

    glm::mat4 mvp = projection * view * model;
    glm::mat4 mv = view * model;

    glUniformMatrix4fv(g_uMVP, 1, GL_FALSE, glm::value_ptr(mvp));
    glUniformMatrix4fv(g_uMV, 1, GL_FALSE, glm::value_ptr(mv));


    if (!makeBody)
        glDrawElements(GL_TRIANGLES, g_model1.indexCount, GL_UNSIGNED_INT, NULL);




}



void cleanup()
{
    if (g_shaderProgram1 != 0)
        glDeleteProgram(g_shaderProgram1);
    if (g_model1.vbo != 0)
        glDeleteBuffers(1, &g_model1.vbo);
    if (g_model1.ibo != 0)
        glDeleteBuffers(1, &g_model1.ibo);
    if (g_model1.vao != 0)
        glDeleteVertexArrays(1, &g_model1.vao);


    glDeleteTextures(1, &texture1_ID);
    glDeleteTextures(1, &texture2_ID);
}

bool initOpenGL()
{
    // Initialize GLFW functions.
    if (!glfwInit())
    {
        cout << "Failed to initialize GLFW" << endl;
        return false;
    }

    // Request OpenGL 3.3 without obsoleted functions.
    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
    glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE);
    glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);

    // Create window.
    g_window = glfwCreateWindow(800, 600, "OpenGL Test", NULL, NULL);
    if (g_window == NULL)
    {
        cout << "Failed to open GLFW window" << endl;
        glfwTerminate();
        return false;
    }

    // Initialize OpenGL context with.
    glfwMakeContextCurrent(g_window);

    // Set internal GLEW variable to activate OpenGL core profile.
    glewExperimental = true;

    // Initialize GLEW functions.
    if (glewInit() != GLEW_OK)
    {
        cout << "Failed to initialize GLEW" << endl;
        return false;
    }

    // Ensure we can capture the escape key being pressed.
    glfwSetInputMode(g_window, GLFW_STICKY_KEYS, GL_TRUE);

    // Set callback for framebuffer resizing event.
    glfwSetFramebufferSizeCallback(g_window, reshape);

    return true;
}

void tearDownOpenGL()
{
    // Terminate GLFW.
    glfwTerminate();
}


static void cursor_position_callback(GLFWwindow* window, double xpos, double ypos)
{
    GLfloat xoffset = xpos - lastX;
    GLfloat yoffset = lastY - ypos;
    lastX = xpos;
    lastY = ypos;

    GLfloat sensitivity = 0.05f;
    xoffset *= sensitivity;
    yoffset *= sensitivity;
    yaw_var += xoffset;
    pitch_var += yoffset;

    if (pitch_var > 179.0f)
        pitch_var = 179.0f;
    if (pitch_var < -179.0f)
        pitch_var = -179.0f;
    //cout << pitch_var << endl;
    cout << xpos << ' ' << ypos << endl;
}




void mouse_button_callback(GLFWwindow* window, int button, int action, int mods)
{
    if (action == 1) {
        if (lastX >= 130 && lastX <= 670 && lastY >= 30 && lastY <= 570) {

            pointsArray[cur_i] = PointCoords((lastX - 130.0 - 270.0) / 270.0, ((600 - lastY) - 30 - 270) / 270.0);
            cur_i++;

        }
    }
}

void key_callback(GLFWwindow* window, int key, int scancode, int action, int mods)
{



}

int main()
{
    for (int i = 0; i < 10; i++)
    {
        pointsArray[i] = PointCoords(-1, -1);
    }
    // Initialize OpenGL
    if (!initOpenGL())
        return -1;

    // Initialize graphical resources.
    bool isOk = init();
    glfwSetInputMode(g_window, GLFW_CURSOR, GLFW_CURSOR_DISABLED);
    glfwSetKeyCallback(g_window, key_callback);
    glfwSetMouseButtonCallback(g_window, mouse_button_callback);
    glfwSetCursorPosCallback(g_window, cursor_position_callback);
    if (isOk)
    {
        // Main loop until window closed or escape pressed.
        while (glfwGetKey(g_window, GLFW_KEY_ESCAPE) != GLFW_PRESS && glfwWindowShouldClose(g_window) == 0)
        {
            // Draw scene.
            draw();
            // Swap buffers.
            glfwSwapBuffers(g_window);
            // Poll window events.
            glfwPollEvents();
        }
    }

    // Cleanup graphical resources.
    cleanup();

    // Tear down OpenGL.
    tearDownOpenGL();

    return isOk ? 0 : -1;
}