#include "Helpers.h"

#include <iostream>
#include <fstream>
#define PI 3.14159265


void VertexArrayObject::init()
{
	glGenVertexArrays(1, &id);
	check_gl_error();
}

void VertexArrayObject::bind()
{
	glBindVertexArray(id);
	check_gl_error();
}

void VertexArrayObject::free()
{
	glDeleteVertexArrays(1, &id);
	check_gl_error();
}

void VertexBufferObject::init()
{
	glGenBuffers(1,&id);
	check_gl_error();
}

void VertexBufferObject::bind()
{
	glBindBuffer(GL_ARRAY_BUFFER,id);
	check_gl_error();
}

void VertexBufferObject::free()
{
	glDeleteBuffers(1,&id);
	check_gl_error();
}

void VertexBufferObject::update(const Eigen::MatrixXf& M)
{
	assert(id != 0);
	glBindBuffer(GL_ARRAY_BUFFER, id);
	glBufferData(GL_ARRAY_BUFFER, sizeof(float)*M.size(), M.data(), GL_DYNAMIC_DRAW);
	rows = M.rows();
	cols = M.cols();
	check_gl_error();
}

bool Program::init(
	const std::string &vertex_shader_string,
	const std::string &fragment_shader_string,
	const std::string &fragment_data_name)
{
	using namespace std;
	vertex_shader = create_shader_helper(GL_VERTEX_SHADER, vertex_shader_string);
	fragment_shader = create_shader_helper(GL_FRAGMENT_SHADER, fragment_shader_string);

	if (!vertex_shader || !fragment_shader)
		return false;

	program_shader = glCreateProgram();

	glAttachShader(program_shader, vertex_shader);
	glAttachShader(program_shader, fragment_shader);

	glBindFragDataLocation(program_shader, 0, fragment_data_name.c_str());
	glLinkProgram(program_shader);

	GLint status;
	glGetProgramiv(program_shader, GL_LINK_STATUS, &status);

	if (status != GL_TRUE)
	{
		char buffer[512];
		glGetProgramInfoLog(program_shader, 512, NULL, buffer);
		cerr << "Linker error: " << endl << buffer << endl;
		program_shader = 0;
		return false;
	}

	check_gl_error();
	return true;
}

void Program::bind()
{
	glUseProgram(program_shader);
	check_gl_error();
}

GLint Program::attrib(const std::string &name) const
{
	return glGetAttribLocation(program_shader, name.c_str());
}

GLint Program::uniform(const std::string &name) const
{
	return glGetUniformLocation(program_shader, name.c_str());
}

GLint Program::bindVertexAttribArray(
				const std::string &name, VertexBufferObject& VBO) const
{
	GLint id = attrib(name);
	if (id < 0)
		return id;
	if (VBO.id == 0)
	{
		glDisableVertexAttribArray(id);
		return id;
	}
	VBO.bind();
	glEnableVertexAttribArray(id);
	
	// check_gl_error();

	glVertexAttribPointer(id, VBO.rows, GL_FLOAT, GL_FALSE, 0, 0);

	check_gl_error();
	

	return id;
}

void Program::free()
{
	if (program_shader)
	{
		glDeleteProgram(program_shader);
		program_shader = 0;
	}
	if (vertex_shader)
	{
		glDeleteShader(vertex_shader);
		vertex_shader = 0;
	}
	if (fragment_shader)
	{
		glDeleteShader(fragment_shader);
		fragment_shader = 0;
	}
	check_gl_error();
}

GLuint Program::create_shader_helper(GLint type, const std::string &shader_string)
{
	using namespace std;
	if (shader_string.empty())
		return (GLuint) 0;

	GLuint id = glCreateShader(type);
	const char *shader_string_const = shader_string.c_str();
	glShaderSource(id, 1, &shader_string_const, NULL);
	glCompileShader(id);

	GLint status;
	glGetShaderiv(id, GL_COMPILE_STATUS, &status);

	if (status != GL_TRUE)
	{
		char buffer[512];
		if (type == GL_VERTEX_SHADER)
			cerr << "Vertex shader:" << endl;
		else if (type == GL_FRAGMENT_SHADER)
			cerr << "Fragment shader:" << endl;
		else if (type == GL_GEOMETRY_SHADER)
			cerr << "Geometry shader:" << endl;
		cerr << shader_string << endl << endl;
		glGetShaderInfoLog(id, 512, NULL, buffer);
		cerr << "Error: " << endl << buffer << endl;
		return (GLuint) 0;
	}
	check_gl_error();

	return id;
}

void _check_gl_error(const char *file, int line)
{
	GLenum err (glGetError());

	while(err!=GL_NO_ERROR)
	{
		std::string error;

		switch(err)
		{
			case GL_INVALID_OPERATION:      error="INVALID_OPERATION";      break;
			case GL_INVALID_ENUM:           error="INVALID_ENUM";           break;
			case GL_INVALID_VALUE:          error="INVALID_VALUE";          break;
			case GL_OUT_OF_MEMORY:          error="OUT_OF_MEMORY";          break;
			case GL_INVALID_FRAMEBUFFER_OPERATION:  error="INVALID_FRAMEBUFFER_OPERATION";  break;
		}

		std::cerr << "GL_" << error.c_str() << " - " << file << ":" << line << std::endl;
		err = glGetError();
	}
}


// returns a matrix that calculate the position of the object in the world based on
// its center, its rotation (around z axis), and scaling
Eigen::Matrix4Xf rotation_scale_matrix(float rz, float rx, float ry, float s, Eigen::Vector3f o){
		rz = PI*rz/180;
		rx = PI*rx/180;
		ry = PI*ry/180;


		Eigen::Matrix3Xf rotz(3, 3);
		rotz<<
		cos(rz), sin(rz), 0,
		-sin(rz), cos(rz), 0,
		0, 0, 1;

		Eigen::Matrix3Xf roty(3, 3);
		roty<<
		cos(ry), 0, -sin(ry),
		0, 1, 0,
		sin(ry), 0, cos(ry);

		Eigen::Matrix3Xf rotx(3, 3);
		rotx<<
		1, 0, 0, 
		0, cos(rx), sin(rx), 
		0, -sin(rx), cos(rx);

		Eigen::Matrix3Xf rot(3, 3);
		rot = rotz*roty*rotx;

		Eigen::Matrix3Xf scale(3, 3);
		scale<<
		s, 0, 0,
		0, s, 0,
		0, 0, s;

		Eigen::Matrix3Xf rs(3, 3);
		rs = rot*scale;




		Eigen::Matrix4Xf trans(4, 4);
		trans.col(0) = Eigen::Vector4f(rs(0, 0), rs(1, 0), rs(2, 0), 0);
		trans.col(1) = Eigen::Vector4f(rs(0, 1), rs(1, 1), rs(2, 1), 0);
		trans.col(2) = Eigen::Vector4f(rs(0, 2), rs(1, 2), rs(2, 2), 0);
		trans.col(3) = Eigen::Vector4f(o(0), o(1), o(2), 1);

		return trans;
}


int have_same_sign(float a, float b, float c){
		if (a<=0 && b<=0 &&c<=0){
				return 1;
		}else if (a>=0 && b>=0 &&c>=0){
				return -1;
		}else{
				return 0;
		}
}





Eigen::Matrix3Xf return_rotation_matrix(float rz, float rx, float ry){
		rz = PI*rz/180;
		rx = PI*rx/180;
		ry = PI*ry/180;


		Eigen::Matrix3Xf rotz(3, 3);
		rotz<<
		cos(rz), sin(rz), 0,
		-sin(rz), cos(rz), 0,
		0, 0, 1;

		Eigen::Matrix3Xf roty(3, 3);
		roty<<
		cos(ry), 0, -sin(ry),
		0, 1, 0,
		sin(ry), 0, cos(ry);

		Eigen::Matrix3Xf rotx(3, 3);
		rotx<<
		1, 0, 0, 
		0, cos(rx), sin(rx), 
		0, -sin(rx), cos(rx);

		Eigen::Matrix3Xf rot(3, 3);
		rot = rotz*roty*rotx;
		return rot;

}

// calculate a face's normal based on the the three vertex and its center
// the idea is: we will first calculate the normal in the normal way
// then we need to check if it is facing outward
// to do this, we can check the sign of three vertex minus the center
// the axis that matters is the dimension that has the same sign
Eigen::Vector3f face_normal(Eigen::Vector3f a, Eigen::Vector3f b, Eigen::Vector3f c){
		Eigen::Vector3f edge1 = b-a;
		Eigen::Vector3f edge2 = c-a;
		Eigen::Vector3f n = edge1.cross(edge2).normalized();

		return n;

}


bool RayIntersectsTriangle(Eigen::Vector3f rayOrigin, 
													 Eigen::Vector3f rayVector, 
													 Eigen::Vector3f vertex0,
													 Eigen::Vector3f vertex1,
													 Eigen::Vector3f vertex2,
													 float& t){
		const double EPSILON = 0.0000001; 
		Eigen::Vector3f edge1, edge2, h, s, q;
		double a,f,u,v;
		edge1 = vertex1 - vertex0;
		edge2 = vertex2 - vertex0;
		h = rayVector.cross(edge2);
		a = edge1.dot(h);
		if (a > -EPSILON && a < EPSILON)
				return false;
		f = 1/a;
		s = rayOrigin - vertex0;
		u = f * (s.dot(h));
		if (u < 0.0 || u > 1.0)
				return false;
		q = s.cross(edge1);
		v = f * rayVector.dot(q);
		if (v < 0.0 || u + v > 1.0)
				return false;
		// At this stage we can compute t to find out where the intersection point is on the line.
		t = f * edge2.dot(q);
		if (t > EPSILON) // ray intersection
		{
				return true;
		}
		else // This means that there is a line intersection but not a ray intersection.
				return false;
}



Eigen::Vector3f return_look_up_vector(Eigen::Vector3f camera_p, Eigen::Vector3f end_p, float rx, float ry){
		Eigen::Matrix3Xf rot = return_rotation_matrix(0, rx, ry);
		Eigen::Vector3f camera_p_rotated = rot*camera_p;
		Eigen::Vector3f end_p_rotated = rot*end_p;
		return (end_p_rotated - camera_p_rotated).normalized();

}


Eigen::Matrix4Xf return_camera_matrix(Eigen::Vector3f e, Eigen::Vector3f end_p, float rx, float ry){
	Eigen::Vector3f t = return_look_up_vector(e, end_p, rx, ry);
	Eigen::Matrix3Xf rot = return_rotation_matrix(0.0, rx, ry);
	e = rot*e;
	end_p = rot*end_p;

	
	Eigen::Vector3f g = (Eigen::Vector3f(0.0, 0.0, 0.0) - e).normalized();
	Eigen::Vector3f w = - (g.normalized());
	Eigen::Vector3f u = (t.cross(w)).normalized();
	Eigen::Vector3f v = w.cross(u);


	// std::cout<<"view_up vector: "<<std::endl;
	// std::cout<<t<<std::endl;

	// std::cout<<"eye position: "<<std::endl;
	// std::cout<<e<<std::endl;

	// std::cout<<"gaze direction: "<<std::endl;
	// std::cout<<g<<std::endl;


	Eigen::Matrix4Xf camera_matrix(4, 4);

	camera_matrix.col(0) = Eigen::Vector4f(u(0), u(1), u(2), 0.0);
	camera_matrix.col(1) = Eigen::Vector4f(v(0), v(1), v(2), 0.0);
	camera_matrix.col(2) = Eigen::Vector4f(w(0), w(1), w(2), 0.0);
	camera_matrix.col(3) = Eigen::Vector4f(e(0), e(1), e(2), 1.0);

	return camera_matrix.inverse();



}


Eigen::Matrix4Xf return_orth_proj_matrix(float s, Eigen::Vector3f e){
	Eigen::Vector3f eye_mapped_to_z = Eigen::Vector3f(0, 0, sqrt(pow(e(0), 2)+ pow(e(1), 2)+ pow(e(2), 2)));
	Eigen::Vector3f bottom_left = Eigen::Vector3f(-1, -1, 3);
	Eigen::Vector3f top_right = Eigen::Vector3f(1, 1, -1) - eye_mapped_to_z;

 


	bottom_left = bottom_left*s - Eigen::Vector3f(0, 0, 2);
	top_right = top_right*s - Eigen::Vector3f(0, 0, 2);

	float r = top_right(0);
	float t = top_right(1);
	float f = top_right(2);
	float l = bottom_left(0);
	float b = bottom_left(1);
	float n = bottom_left(2);

	Eigen::Matrix4Xf proj(4, 4);
	proj<<
	2/(r-l), 0, 0, -(r+l)/(r-l),
	0, 2/(t-b), 0, -(t+b)/(t-b),
	0, 0, -2/(n-f), (n+f)/(n-f),
	0, 0, 0, 1;
	return proj;
}



Eigen::Matrix4Xf return_pers_proj_matrix(float s, Eigen::Vector3f e){

	Eigen::Vector3f eye_mapped_to_z = Eigen::Vector3f(0, 0, sqrt(pow(e(0), 2)+ pow(e(1), 2)+ pow(e(2), 2)));  
	Eigen::Vector3f bottom_left = Eigen::Vector3f(-1, -1, 1 );
	Eigen::Vector3f top_right = Eigen::Vector3f(1, 1, -1) - eye_mapped_to_z;



	bottom_left = bottom_left*s - Eigen::Vector3f(0, 0, 2);
	top_right = top_right*s - Eigen::Vector3f(0, 0, 2);

	float r = top_right(0);
	float t = top_right(1);
	float f = top_right(2);
	float l = bottom_left(0);
	float b = bottom_left(1);
	float n = bottom_left(2);

	Eigen::Matrix4Xf proj(4, 4);
	proj<<
	2*abs(n)/(r-l), 0, (r+l)/(r - l), 0,
	0, 2*abs(n)/(t-b), (b+t)/(t - b), 0,
	0, 0, (abs(f)+abs(n))/(abs(n)-abs(f)), 2*abs(f)*abs(n)/(abs(n) - abs(f)),
	0, 0, -1, 0;
	return proj;
}





// check if ray intersects the bounding box
bool RayIntersectBox(Eigen::Vector3f center, float rz, float rx, float ry, float s, Eigen::Matrix3Xf cube, Eigen::Vector3f rayOrigin, Eigen::Vector3f rayVector, float& min_t){
	Eigen::Matrix4Xf rot_scale(4, 4);
	float t = std::numeric_limits<float>::max();
	bool found_intersect = false;
	rot_scale = rotation_scale_matrix(rz, rx, ry, s, center);
	for (int i=0; i<cube.cols(); i+=3){
		Eigen::Vector4f pa, pb, pc;
		pa = Eigen::Vector4f(cube.col(i)(0), cube.col(i)(1), cube.col(i)(2), 1);
		pb = Eigen::Vector4f(cube.col(i+1)(0), cube.col(i+1)(1), cube.col(i+1)(2), 1);
		pc = Eigen::Vector4f(cube.col(i+2)(0), cube.col(i+2)(1), cube.col(i+2)(2), 1);
		pa = rot_scale*pa;
		pb = rot_scale*pb;
		pc = rot_scale*pc;
		Eigen::Vector3f a, b, c;
		a = Eigen::Vector3f(pa(0), pa(1), pa(2));
		b = Eigen::Vector3f(pb(0), pb(1), pb(2));
		c = Eigen::Vector3f(pc(0), pc(1), pc(2));
		if (RayIntersectsTriangle(rayOrigin, 
													 rayVector, 
													 a,
													 b,
													 c,
													 t)){
			if (t<min_t){
				min_t = t;
				found_intersect = true;
			}
		}
	}
	return found_intersect;
}

int orientation(Eigen::Vector3f p1, Eigen::Vector3f p2, Eigen::Vector3f p3, Eigen::Vector3f v) 
{ 
	Eigen::Vector3f N = (p2-p1).cross(p3 - p1);
	float w = N.dot(p1 - v);

	if (w == 0){
		return 0;
	}
	if (w>0){
		return 1;
	}
	if (w<0){
		return 2;
	}

    // float val = (p2(1) - p1(1)) * (p3(0) - p2(0)) - 
    //           (p2(0) - p1(0)) * (p3(1) - p2(1)); 
  
    // if (val == 0) return 0;  // colinear 
  
    // return (val > 0)? 1: 2; // clock or counterclock wise 
} 


Eigen::Vector3f complementary_color(Eigen::Vector3f c){
	float r, g, b;
	r = c(0);
	g = c(1);
	b = c(2);

	float max_c = std::max(b, std::max(r, g));
	float min_c = std::min(b, std::min(r, g));
	Eigen::Vector3f comp = Eigen::Vector3f((max_c == r || min_c == r)?1.0 - r:r, (max_c == g || min_c == g)?1.0 - g:g, (max_c == b || min_c == b)?1.0 - b:b);


	return comp;
}





// load the number off object, and get their width
void load_number(std::string file, Eigen::Matrix3Xf& number_edges, float& number_width){
	Eigen::Matrix3Xf number_shape_points;
	std::ifstream shape_in(file);
	std::string dummyLine;
	getline(shape_in, dummyLine);
	int num_vert, num_face, num_edge;
	shape_in>>num_vert>>num_face>>num_edge;

	// to get the bounding box of a number
	float x_max = - std::numeric_limits<float>::infinity();
	float x_min = std::numeric_limits<float>::infinity();
	float y_max = - std::numeric_limits<float>::infinity();
	float y_min = std::numeric_limits<float>::infinity();
	float x_center = 0;
	float y_center = 0;
	float z_center = -0.1;

	for (int i=0; i<num_vert; i++){
		float x, y, z;
		shape_in>>x>>y>>z;
		Eigen::Vector3f point(x, y, z);
		number_shape_points.conservativeResize(3, number_shape_points.cols()+1);
		number_shape_points.col(number_shape_points.cols()-1) = point;

		if (x>x_max){
				x_max = x;
		}
		if (x<x_min){
				x_min = x;
		}
		if (y>y_max){
				y_max = y;
		}
		if (y<y_min){
				y_min = y;
		}
	}

	x_center = (x_max + x_min)/2;
	y_center = (y_max + y_min)/2;
	Eigen::Vector3f shape_center(x_center, y_center, z_center);

	float y_max_scale = std::max(abs(y_max - y_center), abs(y_center - y_min));
	float x_max_scale = std::max(abs(x_max - x_center), abs(x_center - x_min));
	float max_scale = std::max(y_max_scale, x_max_scale);

	number_width = ((x_max - x_min)/max_scale);


	// now each point should be in a 2 by 2 box
	for (int i=0; i<number_shape_points.cols(); i++){
		Eigen::Vector3f point = number_shape_points.col(i);
		point = (point - shape_center)/max_scale;
		number_shape_points.col(i) = point;
		// number_shape_points.col(i)(2) = 1.0;
		// std::cout<<number_shape_points.col(i)(0)<<" "<<number_shape_points.col(i)(1)<<" "<<number_shape_points.col(i)(2)<<std::endl;
	}

	for (int i=0; i<num_face; i++){
		int n_edge;
		int x, y, z;
		shape_in>>n_edge>>x>>y>>z;
		Eigen::Vector3f px, py, pz;



		number_edges.conservativeResize(3, number_edges.cols()+1);
		px = number_shape_points.col(x);
		number_edges.col(number_edges.cols()-1) = px;

		number_edges.conservativeResize(3, number_edges.cols()+1);
		py = number_shape_points.col(y);
		number_edges.col(number_edges.cols()-1) = py;

		number_edges.conservativeResize(3, number_edges.cols()+1);
		pz = number_shape_points.col(z);
		number_edges.col(number_edges.cols()-1) = pz;
	}
}


void load_all_number(std::vector<Eigen::Matrix3Xf>& all_numbers_edges, std::vector<float>& all_numbers_width){
	for (int i=0; i<10; i++){
		Eigen::Matrix3Xf number_edge;
		float number_width;
		std::string file = "../numbers_off/"+std::to_string(i)+".off";
		load_number(file, number_edge, number_width);
		all_numbers_edges.push_back(number_edge);
		all_numbers_width.push_back(number_width);
	}
}


bool LineIntersection(Eigen::Vector3f p1, Eigen::Vector3f p2, Eigen::Vector3f p3, Eigen::Vector3f p4, float& u){
	Eigen::Vector3f da = p2 - p1;
	Eigen::Vector3f db = p4 - p3;
	Eigen::Vector3f dc = p3 - p1;

	if (abs(dc.dot(da.cross(db))) !=0){
		return false;
	}

	u = dc.cross(db).dot(da.cross(db)) / pow(da.cross(db).norm(), 2);
	// float base = (p4(1) - p3(1))*(p2(0) - p1(0)) - (p4(0) - p3(0))*(p2(1) - p1(1));
	// if (base == 0){
	// 	return false;
	// }

	// float top = (p4(0) - p3(0))*(p1(1) - p3(1)) - (p4(1) - p3(1))*(p1(0) - p3(0));
	// u = top/base;
	if (0<=u && u<=1){
		return true;
	}else{
		return false;
	}
	
}




Eigen::Matrix3Xf GetMaximumSquare(Eigen::Matrix3Xf bounding_box, Eigen::Vector3f center){
	Eigen::Matrix3Xf tmp;

	float x_max = - std::numeric_limits<float>::infinity();
	float x_min = std::numeric_limits<float>::infinity();
	float y_max = - std::numeric_limits<float>::infinity();
	float y_min = std::numeric_limits<float>::infinity();

	for (int i=0; i<bounding_box.cols(); i++){
		Eigen::Vector3f bounding_point = bounding_box.col(i);
		float x = bounding_point(0);
		float y = bounding_point(1);
		if (x>x_max)
			x_max = x;
		if (x<x_min)
			x_min = x;
		if (y>y_max)
			y_max = y;
		if (y<y_min)
			y_min = y;

	}

	float center_x = center(0);
	float center_y = center(1);

	float x_max_range = std::max(x_max - center_x, center_x - x_min);
	float y_max_range = std::max(y_max - center_y, center_y - y_min);

	float max_half_length = std::max(y_max_range, x_max_range);


	// define the four corner of the bounding square of the cell
	Eigen::Vector3f top_right, top_left, bot_right, bot_left;
	top_right = center + Eigen::Vector3f(sqrt(2) * max_half_length, sqrt(2) * max_half_length, 0);
	top_left = center + Eigen::Vector3f(-sqrt(2) * max_half_length, sqrt(2) * max_half_length, 0);
	bot_right = center + Eigen::Vector3f(sqrt(2) * max_half_length, -sqrt(2) * max_half_length, 0);
	bot_left = center + Eigen::Vector3f(-sqrt(2) * max_half_length, -sqrt(2) * max_half_length, 0);

	float u_min = 1.01;
	for (int i=0; i<bounding_box.cols(); i++){
		float u = 0;
		if (LineIntersection(center, top_right, bounding_box.col(i), bounding_box.col((i+1) % (int)(bounding_box.cols())), u) ){
			if (u<u_min){
				u_min = u;
			}
		}
		if (LineIntersection(center, bot_right, bounding_box.col(i), bounding_box.col((i+1) % (int)(bounding_box.cols())), u) ){
			if (u<u_min){
				u_min = u;
			}
		}
		if (LineIntersection(center, top_left, bounding_box.col(i), bounding_box.col((i+1) % (int)(bounding_box.cols())), u) ){
			if (u<u_min){
				u_min = u;
			}
		}
		if (LineIntersection(center, bot_left, bounding_box.col(i), bounding_box.col((i+1) % (int)(bounding_box.cols())), u) ){
			if (u<u_min){
				u_min = u;
			}
		}
	}

	tmp.resize(3, 4);
	top_right = center + u_min * Eigen::Vector3f(sqrt(2) * max_half_length, sqrt(2) * max_half_length, bounding_box.col(0)(2)+0.01);
	top_left = center + u_min * Eigen::Vector3f(-sqrt(2) * max_half_length, sqrt(2) * max_half_length, bounding_box.col(0)(2)+0.01);
	bot_right = center + u_min * Eigen::Vector3f(sqrt(2) * max_half_length, -sqrt(2) * max_half_length, bounding_box.col(0)(2)+0.01);
	bot_left = center + u_min * Eigen::Vector3f(-sqrt(2) * max_half_length, -sqrt(2) * max_half_length, bounding_box.col(0)(2)+0.01);

	// top_right(2) = -0.5;
	// top_left(2) = -0.5;
	// bot_right(2) = -0.5;
	// bot_left(2) = -0.5;

	tmp.col(0) = top_right;
	tmp.col(1) = top_left;
	tmp.col(2) = bot_left;
	tmp.col(3) = bot_right;

	return tmp;
}



void getNumberCenterSize(Eigen::Matrix3Xf square, 
						 int number, 
						 Eigen::Vector3f square_center, 
						 std::vector<Eigen::Vector3f>& number_centers, 
						 std::vector<float>& each_number_size, 
						 std::vector<float> all_number_width){

	Eigen::Vector3f p1 = square.col(0);
	Eigen::Vector3f p2 = square.col(1);
	
	float square_width = std::max(abs(p1(0) - p2(0)), abs(p1(1) - p2(1)));
	if (number < 10){
		number_centers.push_back(square_center);
		each_number_size.push_back(square_width/2);
	}else{
		// std::cout<<number<<std::endl;
		int number_length = std::to_string(number).length();
		float number_width = 0;
		for (int i=0; i<number_length; i++){
			number_width += all_number_width[(int)(std::to_string(number)[i]) - 48] + 2.0/(number_length + 1 + 5);
			// std::cout<<(std::to_string(number)[i])<<" " <<all_number_width[(int)(std::to_string(number)[i]) - 48]<<std::endl;
		}


		if (number_width <= 2){

			Eigen::Vector3f number_start_boundary = square_center;
			number_start_boundary(0) -= number_width * (square_width/2) / 2;



			for (int i=0; i<number_length; i++){
				number_start_boundary(0) += (2.0/(number_length + 5)) * (square_width / 2) / 2;
				number_start_boundary(0) += all_number_width[(int)(std::to_string(number)[i]) - 48]/4 * square_width;
				number_centers.push_back(number_start_boundary);
				number_start_boundary(0) += all_number_width[(int)(std::to_string(number)[i]) - 48] * (square_width / 2) / 2;
				each_number_size.push_back(square_width/2);
			}

		}else{

			Eigen::Vector3f number_start_boundary = square_center;
			number_start_boundary(0) -= square_width/2;
			for (int i=0; i<number_length; i++){
				number_start_boundary(0) += (2.0/(number_length + 5)) * (square_width / number_width) / 2;
				number_start_boundary(0) += all_number_width[(int)(std::to_string(number)[i]) - 48] * (square_width / number_width) / 2;
				number_centers.push_back(number_start_boundary);
				number_start_boundary(0) += all_number_width[(int)(std::to_string(number)[i]) - 48] * (square_width / number_width) / 2;
				each_number_size.push_back(square_width/number_width);
			}
		}
		// std::cout<<"333333333333333333333333333333333333333333333333333333333333333333"<<std::endl;
	}

}





// load voronoi points and all the other maps
void load_voronoi(std::string file, 
				  int number_of_mines, 
				  std::vector<Eigen::Matrix3Xf>& cells, 
				  std::vector<Eigen::Vector3f>& cell_centers,
				  std::vector<Eigen::Matrix3Xf>& cell_center_squares,
				  std::map<int, std::set<int> >& neighbors, 
				  std::set<int>& mine_cells, 
				  std::map<int, int>& neighbor_mines,
				  std::map<int, std::set<int> >& square_cell_maps){

	Eigen::Matrix3Xf voronoi_points;
	std::ifstream shape_in(file);
	int num_vert, num_face;
	shape_in>>num_vert>>num_face;
	// std::cout<<num_face<<std::endl;




	// define the square cells
	// the number of square cells is square_cell_partition**2
	int square_cell_partition = (int)floor(cbrt(num_face));
	float square_cell_length = 2.0/(1.0*square_cell_partition);


	// initialize the square cell map
	for (int i=0; i<pow(square_cell_partition, 2); i++){
		std::set<int> tmp;
		square_cell_maps.insert(std::make_pair(i, tmp));
	}
	// how the cell numbers look like:
	// | 6 | 7 | 8 |
	// | 3 | 4 | 5 |
	// | 0 | 1 | 2 |
	// bottom right is negative, top right is positive

	// get the points
	// the original points are from 0, 0 to 1, 1
	// we map it to from -1, -1 to 1, 1


	// map, key is the vertex index, value is a set, which contains the cell index
	std::map<int, std::set<int> > point_to_cell;



	for (int i=0; i<num_vert; i++){
		double x, y;
		shape_in>>x>>y;
		std::set<int> tmp;
		// there can be floating points
		if (y>1){
			y = 1.0;
		}
		if (y<0){
			y = 0.0;
		}
		if (x<0){
			x = 0.0;
		}
		if (x>1){
			x = 1.0;
		}
		Eigen::Vector3f point(x*2-1, y*2-1, 0.0);
		voronoi_points.conservativeResize(3, voronoi_points.cols()+1);
		voronoi_points.col(voronoi_points.cols()-1) = point;
		point_to_cell.insert(std::make_pair(i, tmp));
	}

	// somehow we need this line
	// my guess is there is a line breaker
	// that we need to skip
	std::string dummy_line;
	getline(shape_in, dummy_line);

	// load each cells
	for (int i=0; i<num_face; i++){
		Eigen::Matrix3Xf cell_boundaries;
		std::string s;
		getline(shape_in, s);

		Eigen::Vector3f cell_center = Eigen::Vector3f(0.0, 0.0, 0.0);
		int n = 0;

		float x_max = - std::numeric_limits<float>::infinity();
		float x_min = std::numeric_limits<float>::infinity();
		float y_max = - std::numeric_limits<float>::infinity();
		float y_min = std::numeric_limits<float>::infinity();

		std::istringstream iss(s);
		for(std::string s; iss >> s; ){
			int b = std::stoi(s);
			Eigen::Vector3f b_point = voronoi_points.col(b);

			float x = b_point(0)+1.0;
			float y = b_point(1)+1.0;

			int row_num = (int)floor(y/square_cell_length);
			if (row_num >= square_cell_partition){
				row_num = square_cell_partition - 1;
			}

			int col_num = (int)floor(x/square_cell_length);
			if (col_num >= square_cell_partition){
				col_num = square_cell_partition - 1;
			}

			int square_cell_index = square_cell_partition * row_num + col_num;



			// for each square cell
			// if one of the boundary points is in that square cell
			// we determine this cell is in that square cell
			square_cell_maps[square_cell_index].insert(i);

			cell_boundaries.conservativeResize(3, cell_boundaries.cols()+1);
			cell_boundaries.col(cell_boundaries.cols()-1) = b_point;

			point_to_cell[b].insert(i);

			cell_center += voronoi_points.col(b);
			n++;
		}
		int offset = 0;
		int order = 0;

		while (order==0){
			order = orientation(cell_boundaries.col(offset), cell_boundaries.col(offset+1), cell_boundaries.col(offset+2), Eigen::Vector3f(0, 0, 10));
			offset++;
		}


		// you have to make the cell points conter colckwise
		if (order == 1){
			Eigen::Matrix3Xf flipped;
			flipped.resize(3, cell_boundaries.cols());
			for (int i=0; i<cell_boundaries.cols(); i++){
				flipped.col(i) = cell_boundaries.col(cell_boundaries.cols()-1-i);
			}
			cell_boundaries = flipped;
		}
		cell_center = cell_center/n;
		// cell_center(0) = (x_max + x_min)/2;
		// cell_center(1) = (y_max + y_min)/2;
		cell_centers.push_back(cell_center);
		cells.push_back(cell_boundaries);
		std::set<int> tmp;
		neighbors.insert(std::make_pair(i, tmp));

		cell_center_squares.push_back(GetMaximumSquare(cell_boundaries, cell_center));
	}

	// generate mine cells
	while (mine_cells.size()<number_of_mines){
		int random_cell_number = rand()%num_face;
		mine_cells.insert(random_cell_number);
	}

	// generate a set of int that represents the cells that has mines
	for (auto it = point_to_cell.begin(); it!=point_to_cell.end(); ++it){
		int point_index = it->first;
		std::set<int> cell_index_set = it->second;

		for (auto it2 = cell_index_set.begin(); it2!=cell_index_set.end(); ++it2){
			int cell_ind = *it2;
			neighbors[cell_ind].insert(cell_index_set.begin(), cell_index_set.end());
		}
	}

	for (int i=0; i<num_face; i++){
		neighbors[i].erase(i);
		neighbor_mines.insert(std::make_pair(i, 0));
		for (auto it = neighbors[i].begin(); it!=neighbors[i].end(); ++it){
			if (mine_cells.find(*it)!=mine_cells.end()){
				neighbor_mines[i]+=1;
			}
		}
	}
}




bool PointInPolygon(Eigen::Matrix3Xf cell, Eigen::Vector3f center, float x, float y){

	center(2) = 0;
	Eigen::Vector3f p2 = Eigen::Vector3f(x, y, 0);

	for (int i=0; i<cell.cols(); i++){
		Eigen::Vector3f p3 = cell.col(i);
		Eigen::Vector3f p4 = cell.col((i+1)%cell.cols());
		p3(2) = 0;
		p4(2) = 0;

		float u_tmp = 0;

		// if two line intersects, return false
		if (LineIntersection(center, p2, p3, p4, u_tmp)){
			return false;
		}
	}
	return true;
}





float RandomFloat(float a, float b) {
    float random = ((float) rand()) / (float) RAND_MAX;
    float diff = b - a;
    float r = random * diff;
    return a + r;
}

void load_shape(std::string file,
               Eigen::Matrix3Xf& shape_edges){


	Eigen::Matrix3Xf shape_points;
	std::ifstream shape_in(file);
	std::string dummyLine;
	getline(shape_in, dummyLine);

	int num_vert, num_face, num_edge;

	float x_max = - std::numeric_limits<float>::infinity();
	float x_min = std::numeric_limits<float>::infinity();
	float y_max = - std::numeric_limits<float>::infinity();
	float y_min = std::numeric_limits<float>::infinity();
	float z_max = - std::numeric_limits<float>::infinity();
	float z_min = std::numeric_limits<float>::infinity();

	float x_center = 0;
	float y_center = 0;
	float z_center = 0;


	shape_in>>num_vert>>num_face>>num_edge;
	for (unsigned i=0; i<num_vert; i++){
		// Vector3d point;
		double x, y, z;
		shape_in>>x>>y>>z;

		if (x>x_max){
				x_max = x;
		}
		if (x<x_min){
				x_min = x;
		}
		if (y>y_max){
				y_max = y;
		}
		if (y<y_min){
				y_min = y;
		}
		if (z>z_max){
				z_max = z;
		}
		if (z<z_min){
				z_min = z;
		}
		x_center+=x;
		y_center+=y;
		z_center+=z;

		Eigen::Vector3f point(x, y, z);
		shape_points.conservativeResize(3, shape_points.cols()+1);
		shape_points.col(shape_points.cols()-1) = point;
	}
	x_center = x_center/(1.0*num_vert);
	y_center = y_center/(1.0*num_vert);
	z_center = z_center/(1.0*num_vert);
	Eigen::Vector3f shape_center(x_center, y_center, z_center);

	float z_max_scale = std::max(abs(z_max - z_center), abs(z_center - z_min));
	float y_max_scale = std::max(abs(y_max - y_center), abs(y_center - y_min));
	float x_max_scale = std::max(abs(x_max - x_center), abs(x_center - x_min));

	float max_scale = std::max(z_max_scale, std::max(y_max_scale, x_max_scale));

	for (int i=0; i<shape_points.cols(); i++){
		Eigen::Vector3f point = shape_points.col(i);
		point = (point - shape_center)/max_scale;
		shape_points.col(i) = point/2;
	}

	for (unsigned i=0; i<num_face; i++){
		int n_edge;
		int x, y, z;
		shape_in>>n_edge>>x>>y>>z;

		Eigen::Vector3f px, py, pz;
		shape_edges.conservativeResize(3, shape_edges.cols()+1);
		px = shape_points.col(x);

		shape_edges.col(shape_edges.cols()-1) = px;
		// std::cout<<tmp(0)<<" "<<tmp(1)<<" "<<tmp(2)<<std::endl;

		shape_edges.conservativeResize(3, shape_edges.cols()+1);
		py = shape_points.col(y);

		shape_edges.col(shape_edges.cols()-1) = py;
		// std::cout<<tmp(0)<<" "<<tmp(1)<<" "<<tmp(2)<<std::endl;

		shape_edges.conservativeResize(3, shape_edges.cols()+1);
		pz = shape_points.col(z);

		shape_edges.col(shape_edges.cols()-1) = pz;
	}

}


std::set<int> open_cells(int ind, std::map<int, int> neighbor_mines, std::set<int> mine_cells, std::map<int, std::set<int> > neighbors, std::set<int> flagged_cells){
	std::queue<int> cell_queue;
	cell_queue.push(ind);

	std::set<int> should_open;
	should_open.insert(ind);

	if (neighbor_mines[ind] !=0){
		return should_open;
	}

	while (cell_queue.size()!=0){
		int c = cell_queue.front();
		cell_queue.pop();
		if (mine_cells.find(c) == mine_cells.end() && flagged_cells.find(c) == flagged_cells.end() ){
			std::set<int> c_neighbor = neighbors[c];
			for (auto it = c_neighbor.begin(); it!=c_neighbor.end(); ++it){
				if (should_open.find(*it) == should_open.end()){
					if (neighbor_mines[*it]!=0){
						should_open.insert(*it);
					}else{
						cell_queue.push(*it);
					}
					
				}
			}
			should_open.insert(c);
		}
	}
	return should_open;
}


Eigen::Matrix3Xf populate_triangle(Eigen::Matrix3Xf t){
	float thresh = 0.01;
	Eigen::Vector3f a = t.col(0);
	Eigen::Vector3f b = t.col(1);
	Eigen::Vector3f c = t.col(2);
	if ((b-a).norm()<thresh && (c-b).norm()<thresh && (a-c).norm()<thresh){
		return t;
	}
	double af, bf;
	Eigen::Vector3f ab = (a+b)/2;
	Eigen::Vector3f bc = (b+c)/2;
	Eigen::Vector3f ca = (c+a)/2;

	ab = ab/ab.norm();
	bc = bc/bc.norm();
	ca = ca/ca.norm();




	Eigen::Matrix3Xf result_final, result1, result2, result3, result4, t1, t2, t3, t4;
	t1.resize(3, 3);
	t2.resize(3, 3);
	t3.resize(3, 3);
	t4.resize(3, 3);

	t1<<a, ab, ca;
	t2<<ab, b, bc;
	t3<<ca, bc, c;
	t4<<ab, bc, ca;

	result1 = populate_triangle(t1);
	result2 = populate_triangle(t2);
	result3 = populate_triangle(t3);
	result4 = populate_triangle(t4);

	result_final.resize(3, result1.cols()+result2.cols()+result3.cols()+result4.cols());
	result_final<<result1, result2, result3, result4;
	return result_final;
}


void get_polar(Eigen::Vector3f p, double& a, double& b){
	a = atan2(p(0), p(2)); // phi
	if (a<0){
		a += 2*PI;
	}
	b = atan2(std::hypot(p(0), p(2)), p(1)); // theta
	if (b<0){
		b += 2*PI;
	}
	

}

Eigen::Vector3f polar_to_sphere(Eigen::Vector3f p, double a, double b){
	Eigen::Vector3f tmp;
	tmp(2) = sin(b)*cos(a);
	tmp(0) = sin(b)*sin(a);
	tmp(1) = cos(b);

	// if (p(0)<0)
	// 	tmp(0)*=-1;
	// if (p(1)<0)
	// 	tmp(1)*=-1;
	// if (p(2)<0)
	// 	tmp(2)*=-1;

	return tmp;
}

Eigen::Matrix3Xf GetMaximumSphericalSquare(Eigen::Matrix3Xf bounding_box, Eigen::Vector3f center){
	float min_d = 3;
	for (int i=0; i<bounding_box.cols(); i++){
		Eigen::Vector3f p1 = bounding_box.col(i);
		Eigen::Vector3f p2 = bounding_box.col((i+1)%bounding_box.cols());
		Eigen::Vector3f l1 = p2 - p1;
		Eigen::Vector3f l2 = center - p1;
		float d = (l1.cross(l2)).norm()/l1.norm();
		if (d<min_d){
			min_d = d;
		}
	}
	Eigen::Matrix3Xf square;
	square.resize(3, 4);
	square.col(0) = Eigen::Vector3f(0, 0, 0) + min_d*Eigen::Vector3f(1, 1, 0);
	square.col(1) = Eigen::Vector3f(0, 0, 0) + min_d*Eigen::Vector3f(-1, 1, 0);
	square.col(2) = Eigen::Vector3f(0, 0, 0) + min_d*Eigen::Vector3f(-1, -1, 0);
	square.col(3) = Eigen::Vector3f(0, 0, 0) + min_d*Eigen::Vector3f(1, -1, 0);
	return square;
}




void load_spherical_voronoi(std::string file,
                            int number_of_mines,
                            std::vector<Eigen::Matrix3Xf>& spherical_cells, 
                            std::vector<Eigen::Matrix3Xf>& spherical_cell_triangles,
                            std::vector<Eigen::Vector3f>& spherical_cell_centers,
                            std::vector<Eigen::Vector3f>& spherical_cell_centers_unnormed,
                            std::vector<Eigen::Matrix3Xf>& sphere_center_squares,
                            std::map<int, std::set<int> >& neighbors,
                            std::set<int>& spherical_mine_cells,
                            std::map<int, int>& spherical_neighbor_mines,
                            std::map<int, std::set<int> >& spherical_layer_maps){
	Eigen::Matrix3Xf spherical_cell_points;
	std::ifstream shape_in(file);
	int num_vert, num_face;
	shape_in>>num_vert>>num_face;
	std::map<int, std::set<int> > point_to_cell;

	for (int i=0; i<num_vert; i++){
		float x, y, z;
		shape_in>>x>>y>>z;
		std::set<int> tmp;
		if (y>1){
			y = 1.0;
		}
		if (y<-1){
			y = -1.0;
		}
		if (x<-1){
			x = -1.0;
		}
		if (x>1){
			x = 1.0;
		}
		if (z>1.0){
			z = 1.0;
		}
		if (z<-1.0){
			z = -1.0;
		}

		Eigen::Vector3f point(x, y, z);
		spherical_cell_points.conservativeResize(3, spherical_cell_points.cols()+1);
		spherical_cell_points.col(spherical_cell_points.cols()-1) = point;
		point_to_cell.insert(std::make_pair(i, tmp));
	}

	std::string dummy_line;
	getline(shape_in, dummy_line);

	for (int i=0; i<num_face; i++){
		Eigen::Matrix3Xf cell_boundaries;
		std::string s;
		getline(shape_in, s);
		Eigen::Vector3f cell_center = Eigen::Vector3f(0.0, 0.0, 0.0);
		int n=0;

		std::istringstream iss(s);

		for(std::string s; iss >> s; ){
			int b = std::stoi(s);
			// std::cout<<b<<" ";
			Eigen::Vector3f b_point = spherical_cell_points.col(b);
			cell_boundaries.conservativeResize(3, cell_boundaries.cols()+1);
			cell_boundaries.col(cell_boundaries.cols()-1) = b_point;
			point_to_cell[b].insert(i);
			cell_center+=b_point;
			n++;
		}
		// std::cout<<std::endl;


		// make them conter clockwise
		int offset = 0;
		int order = 0;
		while (order == 0){
			order = orientation(cell_boundaries.col(offset), cell_boundaries.col(offset+1), cell_boundaries.col(offset+2), Eigen::Vector3f(0.0, 0.0, 0.0));
			offset++;
		}
		if (order == 2){
			Eigen::Matrix3Xf flipped;
			flipped.resize(3, cell_boundaries.cols());
			for (int i=0; i<cell_boundaries.cols(); i++){
				flipped.col(i) = cell_boundaries.col(cell_boundaries.cols()-1-i);
			}
			cell_boundaries = flipped;
		}

		cell_center = cell_center / n;

		Eigen::Vector3f cell_center_unnormed = cell_center;

		cell_center = cell_center / cell_center.norm();

		spherical_cell_centers.push_back(cell_center);
		// spherical_cell_centers_unnormed.push_back(cell_center_unnormed);
		std::set<int> tmp;
		neighbors.insert(std::make_pair(i, tmp));
		spherical_cells.push_back(cell_boundaries);

		Eigen::Matrix3Xf cell_triangles;
		for (int k=0; k<cell_boundaries.cols(); k++){
			Eigen::Matrix3Xf tmp;
			tmp.resize(3, 3);
			tmp<<cell_center, cell_boundaries.col(k), cell_boundaries.col((k+1)%cell_boundaries.cols());
			Eigen::Matrix3Xf triangles = populate_triangle(tmp);

			Eigen::Matrix3Xf tmp_result;
			tmp_result.resize(3, cell_triangles.cols()+triangles.cols());
			tmp_result<<cell_triangles, triangles;
			cell_triangles = tmp_result;

		}
		spherical_cell_triangles.push_back(cell_triangles);


		Eigen::Matrix3Xf square;
		square = GetMaximumSphericalSquare(cell_boundaries, cell_center_unnormed);
		sphere_center_squares.push_back(square);


	}

	while(spherical_mine_cells.size()<number_of_mines){
		int random_cell_number = rand()%num_face;
		spherical_mine_cells.insert(random_cell_number);
	}

	for (auto it = point_to_cell.begin(); it!=point_to_cell.end(); ++it){
		int point_index = it->first;
		std::set<int> cell_index_set = it->second;

		for (auto it2 = cell_index_set.begin(); it2!=cell_index_set.end(); ++it2){
			int cell_ind = *it2;
			neighbors[cell_ind].insert(cell_index_set.begin(), cell_index_set.end());
		}
	}

	for (int i=0; i<num_face; i++){
		neighbors[i].erase(i);
		spherical_neighbor_mines.insert(std::make_pair(i, 0));
		for (auto it = neighbors[i].begin(); it!=neighbors[i].end(); ++it){
			if (spherical_mine_cells.find(*it)!=spherical_mine_cells.end()){
				spherical_neighbor_mines[i]+=1;
			}
		}
	}

	for (int i=0; i<num_face; i++){
		float x, y, z;
		shape_in>>x>>y>>z;
		// std::cout<<"x: "<<x<<" y: "<<y<<" z: "<<z<<std::endl;
		Eigen::Vector3f tmp = Eigen::Vector3f(x, y, z);
		spherical_cell_centers_unnormed.push_back(tmp/tmp.norm());

	}








}