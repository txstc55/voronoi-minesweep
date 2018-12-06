#ifndef SHADER_H
#define SHADER_H

#include <iostream>
#include <string>
#include <vector>
#include <Eigen/Core>
#include <algorithm>
#include <string>
#include <Eigen/Geometry>
#include <fstream>
#include <tuple>
#include <map>
#include <math.h>
#include <sstream>
#include <iostream>
#include <imgui/imgui.h>
#include <imgui/imgui_impl_opengl3.h>
#include <imgui/imgui_impl_glfw.h>
#include <queue>

#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp>
#include <set>

#include <ft2build.h>
#include FT_FREETYPE_H

#ifdef _WIN32
#  include <windows.h>
#  undef max
#  undef min
#  undef DrawText
#endif

#ifndef __APPLE__
#  define GLEW_STATIC
#  include <GL/glew.h>
#endif

#ifdef __APPLE__
#   include <OpenGL/gl3.h>
#   define __gl_h_ /* Prevent inclusion of the old gl.h */
#else
#   ifdef _WIN32
#       include <windows.h>
#   endif
#   include <GL/gl.h>
#endif


#define PI 3.14159265



class Shader
{
public:
    GLuint Program;
    // Constructor generates the shader on the fly
    Shader( const GLchar *vertexPath, const GLchar *fragmentPath )
    {
        // 1. Retrieve the vertex/fragment source code from filePath
        std::string vertexCode;
        std::string fragmentCode;
        std::ifstream vShaderFile;
        std::ifstream fShaderFile;
        // ensures ifstream objects can throw exceptions:
        vShaderFile.exceptions ( std::ifstream::badbit );
        fShaderFile.exceptions ( std::ifstream::badbit );
        try
        {
            // Open files
            vShaderFile.open( vertexPath );
            fShaderFile.open( fragmentPath );
            std::stringstream vShaderStream, fShaderStream;
            // Read file's buffer contents into streams
            vShaderStream << vShaderFile.rdbuf( );
            fShaderStream << fShaderFile.rdbuf( );
            // close file handlers
            vShaderFile.close( );
            fShaderFile.close( );
            // Convert stream into string
            vertexCode = vShaderStream.str( );
            fragmentCode = fShaderStream.str( );
        }
        catch ( std::ifstream::failure e )
        {
            std::cout << "ERROR::SHADER::FILE_NOT_SUCCESFULLY_READ" << std::endl;
        }
        const GLchar *vShaderCode = vertexCode.c_str( );
        const GLchar *fShaderCode = fragmentCode.c_str( );
        // 2. Compile shaders
        GLuint vertex, fragment;
        GLint success;
        GLchar infoLog[512];
        // Vertex Shader
        vertex = glCreateShader( GL_VERTEX_SHADER );
        glShaderSource( vertex, 1, &vShaderCode, NULL );
        glCompileShader( vertex );
        // Print compile errors if any
        glGetShaderiv( vertex, GL_COMPILE_STATUS, &success );
        if ( !success )
        {
            glGetShaderInfoLog( vertex, 512, NULL, infoLog );
            std::cout << "ERROR::SHADER::VERTEX::COMPILATION_FAILED\n" << infoLog << std::endl;
        }
        // Fragment Shader
        fragment = glCreateShader( GL_FRAGMENT_SHADER );
        glShaderSource( fragment, 1, &fShaderCode, NULL );
        glCompileShader( fragment );
        // Print compile errors if any
        glGetShaderiv( fragment, GL_COMPILE_STATUS, &success );
        if ( !success )
        {
            glGetShaderInfoLog( fragment, 512, NULL, infoLog );
            std::cout << "ERROR::SHADER::FRAGMENT::COMPILATION_FAILED\n" << infoLog << std::endl;
        }
        // Shader Program
        this->Program = glCreateProgram( );
        glAttachShader( this->Program, vertex );
        glAttachShader( this->Program, fragment );
        glLinkProgram( this->Program );
        // Print linking errors if any
        glGetProgramiv( this->Program, GL_LINK_STATUS, &success );
        if (!success)
        {
            glGetProgramInfoLog( this->Program, 512, NULL, infoLog );
            std::cout << "ERROR::SHADER::PROGRAM::LINKING_FAILED\n" << infoLog << std::endl;
        }
        // Delete the shaders as they're linked into our program now and no longer necessery
        glDeleteShader( vertex );
        glDeleteShader( fragment );
        
    }
    // Uses the current shader
    void Use( )
    {
        glUseProgram( this->Program );
    }
};












class VertexArrayObject
{
public:
    unsigned int id;

    VertexArrayObject() : id(0) {}

    // Create a new VAO
    void init();

    // Select this VAO for subsequent draw calls
    void bind();

    // Release the id
    void free();
};

class VertexBufferObject
{
public:
    typedef unsigned int GLuint;
    typedef int GLint;

    GLuint id;
    GLuint rows;
    GLuint cols;

    VertexBufferObject() : id(0), rows(0), cols(0) {}

    // Create a new empty VBO
    void init();

    // Updates the VBO with a matrix M
    void update(const Eigen::MatrixXf& M);

    // Select this VBO for subsequent draw calls
    void bind();

    // Release the id
    void free();
};

// This class wraps an OpenGL program composed of two shaders
class Program
{
public:
  typedef unsigned int GLuint;
  typedef int GLint;

  GLuint vertex_shader;
  GLuint fragment_shader;
  GLuint program_shader;

  Program() : vertex_shader(0), fragment_shader(0), program_shader(0) { }

  // Create a new shader from the specified source strings
  bool init(const std::string &vertex_shader_string,
  const std::string &fragment_shader_string,
  const std::string &fragment_data_name);

  // Select this shader for subsequent draw calls
  void bind();

  // Release all OpenGL objects
  void free();

  // Return the OpenGL handle of a named shader attribute (-1 if it does not exist)
  GLint attrib(const std::string &name) const;

  // Return the OpenGL handle of a uniform attribute (-1 if it does not exist)
  GLint uniform(const std::string &name) const;

  // Bind a per-vertex array attribute
  GLint bindVertexAttribArray(const std::string &name, VertexBufferObject& VBO) const;

  GLuint create_shader_helper(GLint type, const std::string &shader_string);

};

// From: https://blog.nobel-joergensen.com/2013/01/29/debugging-opengl-using-glgeterror/
void _check_gl_error(const char *file, int line);

///
/// Usage
/// [... some opengl calls]
/// glCheckError();
///



// return the matrix based on center, rotation and scaling
Eigen::Matrix4Xf rotation_scale_matrix(float rz, float rx, float ry, float s, Eigen::Vector3f o);

Eigen::Vector3f return_look_up_vector(Eigen::Vector3f camera_p, Eigen::Vector3f end_p, float rx, float ry);


// Eigen::Matrix4Xf camera_matrix(Eigen::Vector3f e, Eigen::Vector3f g);

// return the normal of the face
Eigen::Vector3f face_normal(Eigen::Vector3f a, Eigen::Vector3f b, Eigen::Vector3f c);

Eigen::Matrix3Xf return_rotation_matrix(float rz, float rx, float ry);

Eigen::Matrix4Xf return_camera_matrix(Eigen::Vector3f camera_p, Eigen::Vector3f end_p, float rx, float ry);


Eigen::Matrix4Xf return_orth_proj_matrix(float s, Eigen::Vector3f e);

Eigen::Matrix4Xf return_pers_proj_matrix(float s, Eigen::Vector3f e);

// check intersection, and record the point
bool RayIntersectsTriangle(Eigen::Vector3f rayOrigin, 
                           Eigen::Vector3f rayVector, 
                           Eigen::Vector3f vertex0,
                           Eigen::Vector3f vertex1,
                           Eigen::Vector3f vertex2,
                           float& t);

bool RayIntersectBox(Eigen::Vector3f center, float rot_z, float rot_x, float rot_y, float scale, Eigen::Matrix3Xf cube, Eigen::Vector3f rayOrigin, Eigen::Vector3f rayVector, float& t);


int orientation(Eigen::Vector3f p1, Eigen::Vector3f p2, Eigen::Vector3f p3);

// given a color, get the complementary color
Eigen::Vector3f complementary_color(Eigen::Vector3f c);


void RenderText(Shader &shader, std::string text, GLfloat x, GLfloat y, GLfloat scale, glm::vec3 color);


// load number meshes
void load_number(std::string file, Eigen::Matrix3Xf& number_edges, float& number_width);


// load all number meshes
void load_all_number(std::vector<Eigen::Matrix3Xf>& all_numbers_edges, std::vector<float>& all_numbers_width);


// get the square that centers at the center of the polygon
// which has maximum area, and does not go out of the polygon
Eigen::Matrix3Xf GetMaximumSquare(Eigen::Matrix3Xf bounding_box, Eigen::Vector3f center);

// get the centers and the size of numbers
// given a square box
void getNumberCenterSize(Eigen::Matrix3Xf square, 
             int number, 
             Eigen::Vector3f square_center, 
             std::vector<Eigen::Vector3f>& number_centers, 
             std::vector<float>& each_number_size, 
             std::vector<float> all_number_width);



// check if two line segments intersects
// the two lines are from p1 to p2, and from p3 to p4
// if intersects, it will be in form of p1 + u * (p3 - p1)
bool LineIntersection(Eigen::Vector3f p1, Eigen::Vector3f p2, Eigen::Vector3f p3, Eigen::Vector3f p4, float& u);


// load voronoi points and all other maps
void load_voronoi(std::string file, 
                  int number_of_mines, 
                  std::vector<Eigen::Matrix3Xf>& cells, 
                  std::vector<Eigen::Vector3f>& cell_centers,
                  std::vector<Eigen::Matrix3Xf>& cell_center_squares,
                  std::map<int, std::set<int> >& neighbors, 
                  std::set<int>& mine_cells, 
                  std::map<int, int>& neighbor_mines,
                  std::map<int, std::set<int> >& square_cell_maps);



bool PointInPolygon(Eigen::Matrix3Xf cell, Eigen::Vector3f center, float x, float y);

void load_shape(std::string file,
               Eigen::Matrix3Xf& shape_edges);


std::set<int> open_cells(int ind, std::map<int, int> neighbor_mines, std::set<int> mine_cells, std::map<int, std::set<int> > neighbors, std::set<int> flagged_cells);

float RandomFloat(float a, float b);

#define check_gl_error() _check_gl_error(__FILE__,__LINE__)
#endif


