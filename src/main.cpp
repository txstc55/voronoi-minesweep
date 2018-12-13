// Std. Includes
#include <iostream>
#include <map>
#include <string>
// GLEW
#define GLEW_STATIC
#include <GL/glew.h>
// GLFW
#include <GLFW/glfw3.h>
// GLM
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp>
// FreeType
#include <ft2build.h>
#include FT_FREETYPE_H
// GL includes
#include "Helpers.h"

// Timer
#include <chrono>

#define PI 3.14159265


Eigen::Vector3f background = Eigen::Vector3f(211.0/255.0, 219.0/255.0, 229.0/255.0);

VertexBufferObject vbo_number;

Eigen::Vector3f camera_position = Eigen::Vector3f(0.0, 0.0, 4.0);
Eigen::Vector3f camera_look_up_endpoint = Eigen::Vector3f(0.0, 1.0, 4.0);


Eigen::Matrix4f view(4,4);
Eigen::Matrix4Xf camera_matrix(4, 4);
Eigen::Matrix4Xf projection_matrix(4, 4);
float camera_scale = 1.0;


// square mine grid elements
int num_of_cells;
int number_of_mines;
std::map<int, std::set<int> > square_cell_maps;
std::vector<Eigen::Matrix3Xf> cells;
std::vector<Eigen::Vector3f> cell_centers;
int square_cell_partition;
std::vector<Eigen::Matrix3Xf> cell_center_squares;
std::map<int, std::set<int> > neighbors;
std::set<int> mine_cells;
std::map<int, int> neighbor_mines;


std::vector<Eigen::Vector3f> cell_colors;
std::vector<Eigen::Vector3f> cell_complementary_colors;


// sphere mine variables

std::vector<Eigen::Matrix3Xf> spherical_cells; 
std::vector<Eigen::Vector3f> spherical_cell_centers;
std::vector<Eigen::Vector3f> spherical_cell_centers_unnormed;
std::vector<Eigen::Matrix3Xf> spherical_cell_triangles;
std::vector<Eigen::Matrix3Xf> sphere_center_squares;
std::set<int> spherical_mine_cells;
std::map<int, int> spherical_neighbor_mines;
std::map<int, std::set<int> > spherical_layer_maps;





// game variables
int left_clicked_cell = -1;
int right_clicked_cell = -1;
int left_button_down = 0;
int right_button_down = 0;
std::set<int> opened_cell;
bool game_over = false;
std::set<int> flagged_cells;
int cursor_position_cell = -1;
bool win_game = false;
auto t_now= std::chrono::high_resolution_clock::now();
auto t_start = std::chrono::high_resolution_clock::now();
int finish_time;
char choice;



// nu,ber variables
std::vector<Eigen::Matrix3Xf> all_numbers_edges;
std::vector<float> all_numbers_width;
std::vector<std::vector<Eigen::Vector3f> > number_centers_for_each_cell;
std::vector<std::vector<float> > each_number_size_for_each_cell;


// predefined commands
std::string command = "python ../src/voronoi.py ";	
std::string file = "../src/voronoi_data";

std::string sphere_command = "python ../src/spherical_voronoi.py ";
std::string sphere_file = "../src/spherical_voronoi_data";



//global settings
float ry = 0;
float rx = 0;
float x_position_last = 0;
float y_position_last = 0;
Eigen::Matrix4Xf rotation;



// Properties
const GLuint WIDTH = 1000, HEIGHT = 600;

/// Holds all state information relevant to a character as loaded using FreeType
struct Character {
	GLuint TextureID;   // ID handle of the glyph texture
	glm::ivec2 Size;    // Size of glyph
	glm::ivec2 Bearing;  // Offset from baseline to left/top of glyph
	GLuint Advance;    // Horizontal offset to advance to next glyph
};

std::map<GLchar, Character> Characters;
GLuint VAO, VBO;

void RenderText(Shader &shader, std::string text, GLfloat x, GLfloat y, GLfloat scale, glm::vec3 color);

void framebuffer_size_callback(GLFWwindow* window, int width, int height)
{
	glViewport(0, 0, width, height);
}

void empty_mouse_button_callback(GLFWwindow* window, int button, int action, int mods){}
void empty_mouse_cursor_position_callback(GLFWwindow* window, double xpos, double ypos){}

void mouse_button_callback(GLFWwindow* window, int button, int action, int mods){

	float square_cell_length = 2.0/(1.0*square_cell_partition);


	if (button == GLFW_MOUSE_BUTTON_LEFT && action == GLFW_PRESS){
		left_button_down = 1;
		double xpos, ypos;
		glfwGetCursorPos(window, &xpos, &ypos);

		// Get the size of the window
		int width, height;
		glfwGetWindowSize(window, &width, &height);

		// Convert screen position to world coordinates

		Eigen::Vector4f p_screen(xpos,height-1-ypos,0,1);
		Eigen::Vector4f p_canonical((p_screen[0]/width)*2-1,(p_screen[1]/height)*2-1, -2,1);
		Eigen::Vector4f p_world = (view* projection_matrix* camera_matrix).inverse()*p_canonical;
		double xworld = p_world[0];
		double yworld = p_world[1];

		if (xworld>=-1 && xworld<=1 && yworld>=-1 && yworld<=1){
			int row_num = (int)floor((yworld+1)/square_cell_length);
			if (row_num >= square_cell_partition){
				row_num = square_cell_partition - 1;
			}

			int col_num = (int)floor((xworld+1)/square_cell_length);
			if (col_num >= square_cell_partition){
				col_num = square_cell_partition - 1;
			}

			int square_cell_index = square_cell_partition * row_num + col_num;

			// std::cout<<"xworld: "<<xworld<<", yworld: "<<yworld<<std::endl;
			// std::cout<<"square cell number: "<<square_cell_index<<std::endl;

			std::set<int> square_cell_contained_cells= square_cell_maps[square_cell_index];

			for (auto it = square_cell_contained_cells.begin(); it!=square_cell_contained_cells.end(); ++it){
				int ind = (*it);

				if (PointInPolygon(cells[ind], cell_centers[ind], xworld, yworld)){
					std::cout<<"Left clicked on cell: "<<ind<<std::endl;
					left_clicked_cell = ind;
					break;
				}
			}
		}else{
			left_clicked_cell = -1;
		}
	}


	if (button == GLFW_MOUSE_BUTTON_LEFT && action == GLFW_RELEASE){
		left_button_down = 0;
		double xpos, ypos;
		glfwGetCursorPos(window, &xpos, &ypos);

		// Get the size of the window
		int width, height;
		glfwGetWindowSize(window, &width, &height);

		// Convert screen position to world coordinates

		Eigen::Vector4f p_screen(xpos,height-1-ypos,0,1);
		Eigen::Vector4f p_canonical((p_screen[0]/width)*2-1,(p_screen[1]/height)*2-1, -2,1);
		Eigen::Vector4f p_world = (view* projection_matrix* camera_matrix).inverse()*p_canonical;
		double xworld = p_world[0];
		double yworld = p_world[1];

		if (xworld>=-1 && xworld<=1 && yworld>=-1 && yworld<=1){
			int row_num = (int)floor((yworld+1)/square_cell_length);
			if (row_num >= square_cell_partition){
				row_num = square_cell_partition - 1;
			}

			int col_num = (int)floor((xworld+1)/square_cell_length);
			if (col_num >= square_cell_partition){
				col_num = square_cell_partition - 1;
			}

			int square_cell_index = square_cell_partition * row_num + col_num;

			// std::cout<<"xworld: "<<xworld<<", yworld: "<<yworld<<std::endl;
			// std::cout<<"square cell number: "<<square_cell_index<<std::endl;

			std::set<int> square_cell_contained_cells= square_cell_maps[square_cell_index];

			for (auto it = square_cell_contained_cells.begin(); it!=square_cell_contained_cells.end(); ++it){
				int ind = (*it);

				if (PointInPolygon(cells[ind], cell_centers[ind], xworld, yworld)){
					std::cout<<"Left released on cell: "<<ind<<std::endl;
					if (ind != left_clicked_cell){
						left_clicked_cell = -1;
						break;
					}


					if (flagged_cells.find(ind) == flagged_cells.end()){

						if (mine_cells.find(ind)!=mine_cells.end() ){
							game_over = true;
							win_game = false;
							finish_time = (int)std::chrono::duration_cast<std::chrono::duration<float>>(t_now - t_start).count();

						}else{
							std::cout<<"Not flagged"<<std::endl;
							std::set<int> should_open = open_cells(ind, neighbor_mines, mine_cells, neighbors, flagged_cells);
							opened_cell.insert(should_open.begin(), should_open.end());
						}
					}
					break;
				}
			}
		}
	}


	if (button == GLFW_MOUSE_BUTTON_RIGHT && action == GLFW_PRESS){
		double xpos, ypos;
		right_button_down = 1;
		glfwGetCursorPos(window, &xpos, &ypos);

		// Get the size of the window
		int width, height;
		glfwGetWindowSize(window, &width, &height);

		// Convert screen position to world coordinates

		Eigen::Vector4f p_screen(xpos,height-1-ypos,0,1);
		Eigen::Vector4f p_canonical((p_screen[0]/width)*2-1,(p_screen[1]/height)*2-1, -2,1);
		Eigen::Vector4f p_world = (view* projection_matrix* camera_matrix).inverse()*p_canonical;
		double xworld = p_world[0];
		double yworld = p_world[1];

		if (xworld>=-1 && xworld<=1 && yworld>=-1 && yworld<=1){
			int row_num = (int)floor((yworld+1)/square_cell_length);
			if (row_num >= square_cell_partition){
				row_num = square_cell_partition - 1;
			}

			int col_num = (int)floor((xworld+1)/square_cell_length);
			if (col_num >= square_cell_partition){
				col_num = square_cell_partition - 1;
			}

			int square_cell_index = square_cell_partition * row_num + col_num;

			// std::cout<<"xworld: "<<xworld<<", yworld: "<<yworld<<std::endl;
			// std::cout<<"square cell number: "<<square_cell_index<<std::endl;

			std::set<int> square_cell_contained_cells= square_cell_maps[square_cell_index];

			for (auto it = square_cell_contained_cells.begin(); it!=square_cell_contained_cells.end(); ++it){
				int ind = (*it);
				right_clicked_cell = ind;
				if (PointInPolygon(cells[ind], cell_centers[ind], xworld, yworld)){
					std::cout<<"Right clicked on cell: "<<ind<<std::endl;
					break;
				}
			}
		}
	}


	if (button == GLFW_MOUSE_BUTTON_RIGHT && action == GLFW_RELEASE){
		double xpos, ypos;
		right_button_down = 0;
		glfwGetCursorPos(window, &xpos, &ypos);

		// Get the size of the window
		int width, height;
		glfwGetWindowSize(window, &width, &height);

		// Convert screen position to world coordinates

		Eigen::Vector4f p_screen(xpos,height-1-ypos,0,1);
		Eigen::Vector4f p_canonical((p_screen[0]/width)*2-1,(p_screen[1]/height)*2-1, -2,1);
		Eigen::Vector4f p_world = (view* projection_matrix* camera_matrix).inverse()*p_canonical;
		double xworld = p_world[0];
		double yworld = p_world[1];

		if (xworld>=-1 && xworld<=1 && yworld>=-1 && yworld<=1){
			int row_num = (int)floor((yworld+1)/square_cell_length);
			if (row_num >= square_cell_partition){
				row_num = square_cell_partition - 1;
			}

			int col_num = (int)floor((xworld+1)/square_cell_length);
			if (col_num >= square_cell_partition){
				col_num = square_cell_partition - 1;
			}

			int square_cell_index = square_cell_partition * row_num + col_num;

			// std::cout<<"xworld: "<<xworld<<", yworld: "<<yworld<<std::endl;
			// std::cout<<"square cell number: "<<square_cell_index<<std::endl;

			std::set<int> square_cell_contained_cells= square_cell_maps[square_cell_index];

			for (auto it = square_cell_contained_cells.begin(); it!=square_cell_contained_cells.end(); ++it){
				int ind = (*it);

				if (PointInPolygon(cells[ind], cell_centers[ind], xworld, yworld)){
					std::cout<<"Right released on cell: "<<ind<<std::endl;
					if (ind != right_clicked_cell){
						right_clicked_cell = -1;
						break;
					}

					if (flagged_cells.find(ind)==flagged_cells.end()){
						if (flagged_cells.size()<number_of_mines && opened_cell.find(ind) == opened_cell.end()){
							flagged_cells.insert(ind);
						}
						

					}else{
						flagged_cells.erase(ind);
					}

					break;
				}
			}
		}
	}

	// if (action == GLFW_PRESS && ((button == GLFW_MOUSE_BUTTON_RIGHT && left_button_down == 1) || (button == GLFW_MOUSE_BUTTON_LEFT && right_button_down == 1))){
	// 	std::cout<<"Both key clicked"<<std::endl;
	// }


	if (flagged_cells.size() + opened_cell.size() == num_of_cells){
		game_over = true;
		win_game = true;
		finish_time = (int)std::chrono::duration_cast<std::chrono::duration<float>>(t_now - t_start).count();
	}


}


void mouse_cursor_position_callback(GLFWwindow* window, double xpos, double ypos){
	float square_cell_length = 2.0/(1.0*square_cell_partition);

	glfwGetCursorPos(window, &xpos, &ypos);

	// Get the size of the window
	int width, height;
	glfwGetWindowSize(window, &width, &height);

	// Convert screen position to world coordinates

	Eigen::Vector4f p_screen(xpos,height-1-ypos,0,1);
	Eigen::Vector4f p_canonical((p_screen[0]/width)*2-1,(p_screen[1]/height)*2-1, -2,1);
	Eigen::Vector4f p_world = (view* projection_matrix* camera_matrix).inverse()*p_canonical;
	double xworld = p_world[0];
	double yworld = p_world[1];

	if (xworld>=-1 && xworld<=1 && yworld>=-1 && yworld<=1){
		int row_num = (int)floor((yworld+1)/square_cell_length);
		if (row_num >= square_cell_partition){
			row_num = square_cell_partition - 1;
		}

		int col_num = (int)floor((xworld+1)/square_cell_length);
		if (col_num >= square_cell_partition){
			col_num = square_cell_partition - 1;
		}

		int square_cell_index = square_cell_partition * row_num + col_num;


		std::set<int> square_cell_contained_cells= square_cell_maps[square_cell_index];

		for (auto it = square_cell_contained_cells.begin(); it!=square_cell_contained_cells.end(); ++it){
			int ind = (*it);
			if (PointInPolygon(cells[ind], cell_centers[ind], xworld, yworld)){
				cursor_position_cell = ind;
				break;
			}
		}
	}else{
		cursor_position_cell = -1;
	}


}


void sphere_mouse_button_callback(GLFWwindow* window, int button, int action, int mods){

	double xpos, ypos;
	glfwGetCursorPos(window, &xpos, &ypos);

	// Get the size of the window
	int width, height;
	glfwGetWindowSize(window, &width, &height);

	// Convert screen position to world coordinates

	Eigen::Vector4f p_screen(xpos,height-1-ypos,0,1);
	Eigen::Vector4f p_canonical((p_screen[0]/width)*2-1,(p_screen[1]/height)*2-1, -2,1);
	Eigen::Vector4f p_world = (view* projection_matrix* camera_matrix).inverse()*p_canonical;
	double xworld = p_world[0];
	double yworld = p_world[1];
	x_position_last = xworld;
	y_position_last = yworld;

	if (button == GLFW_MOUSE_BUTTON_LEFT && action == GLFW_PRESS){
		left_button_down = 1;

		if (xworld>=-1 && xworld<=1 && yworld>=-1 && yworld<=1){
			float t_min = 2;
			for (int i=0; i<spherical_cells.size(); i++){

				Eigen::Vector4f c = rotation * Eigen::Vector4f(spherical_cell_centers_unnormed[i](0), spherical_cell_centers_unnormed[i](1), spherical_cell_centers_unnormed[i](2), 1.0);
				// std::cout<<c<<std::endl;
				Eigen::Vector3f p = Eigen::Vector3f(xworld, yworld, sqrt(1 - pow(xworld, 2) - pow(yworld, 2)));
				Eigen::Vector3f ct = Eigen::Vector3f(c(0), c(1), c(2));
				// std::cout<<ct<<std::endl;
				float t = (p - ct).norm();
				if (t<t_min){
					t_min = t;
					left_clicked_cell = i;
				}
			}

		}else{
			left_clicked_cell = -1;
		}

	}
	if (button == GLFW_MOUSE_BUTTON_LEFT && action == GLFW_RELEASE){
		left_button_down = 0;
		if (xworld>=-1 && xworld<=1 && yworld>=-1 && yworld<=1){
			int ind = 2;
			float t_min = 2;
			for (int i=0; i<spherical_cells.size(); i++){

				Eigen::Vector4f c = rotation * Eigen::Vector4f(spherical_cell_centers_unnormed[i](0), spherical_cell_centers_unnormed[i](1), spherical_cell_centers_unnormed[i](2), 1.0);
				// std::cout<<c<<std::endl;
				Eigen::Vector3f p = Eigen::Vector3f(xworld, yworld, sqrt(1 - pow(xworld, 2) - pow(yworld, 2)));
				Eigen::Vector3f ct = Eigen::Vector3f(c(0), c(1), c(2));
				// std::cout<<ct<<std::endl;
				float t = (p - ct).norm();
				if (t<t_min){
					t_min = t;
					ind = i;
				}
			}


			if (ind != left_clicked_cell){
				left_clicked_cell = -1;
			}else if (flagged_cells.find(ind) == flagged_cells.end()){
				if (spherical_mine_cells.find(ind)!=spherical_mine_cells.end()){
					game_over = true;
					win_game = false;
					finish_time = (int)std::chrono::duration_cast<std::chrono::duration<float>>(t_now - t_start).count();

				}else{
					std::cout<<"Not flagged"<<std::endl;
					std::cout<<"Clicked on spherical cell: "<<ind<<std::endl;
					std::set<int> should_open = open_cells(ind, spherical_neighbor_mines, spherical_mine_cells, neighbors, flagged_cells);
					opened_cell.insert(should_open.begin(), should_open.end());
				}
			}

		}else{
			left_clicked_cell = -1;
		}

	}

	if (button == GLFW_MOUSE_BUTTON_RIGHT && action == GLFW_PRESS){
		double xpos, ypos;
		right_button_down = 1;


		if (xworld>=-1 && xworld<=1 && yworld>=-1 && yworld<=1){
			float t_min = 2;
			for (int i=0; i<spherical_cells.size(); i++){

				Eigen::Vector4f c = rotation * Eigen::Vector4f(spherical_cell_centers_unnormed[i](0), spherical_cell_centers_unnormed[i](1), spherical_cell_centers_unnormed[i](2), 1.0);
				// std::cout<<c<<std::endl;
				Eigen::Vector3f p = Eigen::Vector3f(xworld, yworld, sqrt(1 - pow(xworld, 2) - pow(yworld, 2)));
				Eigen::Vector3f ct = Eigen::Vector3f(c(0), c(1), c(2));
				// std::cout<<ct<<std::endl;
				float t = (p - ct).norm();
				if (t<t_min){
					t_min = t;
					right_clicked_cell = i;
				}
			}

		}else{
			right_clicked_cell = -1;
		}
	}

	if (button == GLFW_MOUSE_BUTTON_RIGHT && action == GLFW_RELEASE){
		right_button_down = 0;
		if (xworld>=-1 && xworld<=1 && yworld>=-1 && yworld<=1){
			int ind = -2;
			float t_min = 2;
			for (int i=0; i<spherical_cells.size(); i++){

				Eigen::Vector4f c = rotation * Eigen::Vector4f(spherical_cell_centers_unnormed[i](0), spherical_cell_centers_unnormed[i](1), spherical_cell_centers_unnormed[i](2), 1.0);
				// std::cout<<c<<std::endl;
				Eigen::Vector3f p = Eigen::Vector3f(xworld, yworld, sqrt(1 - pow(xworld, 2) - pow(yworld, 2)));
				Eigen::Vector3f ct = Eigen::Vector3f(c(0), c(1), c(2));
				// std::cout<<ct<<std::endl;
				float t = (p - ct).norm();
				if (t<t_min){
					t_min = t;
					ind = i;
				}
			}


			if (flagged_cells.find(ind)==flagged_cells.end()){
				if (flagged_cells.size()<number_of_mines && opened_cell.find(ind) == opened_cell.end()){
					flagged_cells.insert(ind);
				}
				

			}else{
				flagged_cells.erase(ind);
			}

		}
	}






	if (flagged_cells.size() + opened_cell.size() == num_of_cells){
		game_over = true;
		win_game = true;
		finish_time = (int)std::chrono::duration_cast<std::chrono::duration<float>>(t_now - t_start).count();
	}
}










void sphere_mouse_cursor_position_callback(GLFWwindow* window, double xpos, double ypos){

	glfwGetCursorPos(window, &xpos, &ypos);

	// Get the size of the window
	int width, height;
	glfwGetWindowSize(window, &width, &height);

	// Convert screen position to world coordinates

	Eigen::Vector4f p_screen(xpos,height-1-ypos,0,1);
	Eigen::Vector4f p_canonical((p_screen[0]/width)*2-1,(p_screen[1]/height)*2-1, -2,1);
	Eigen::Vector4f p_world = (view* projection_matrix* camera_matrix).inverse()*p_canonical;
	double xworld = p_world[0];
	double yworld = p_world[1];
	ry = -(xworld - x_position_last)*PI*10;
	rx = (yworld - y_position_last)*PI*10;


	if (left_button_down){
		rotation = rotation_scale_matrix(0.0, rx, ry, 1.0, Eigen::Vector3f(0.0, 0.0, 0.0))*rotation;
		y_position_last = yworld;
		x_position_last = xworld;
	}else{

		// here we will just use the property of delaunay triangle
		// to determine where we clicked
		if (xworld>=-1 && xworld<=1 && yworld>=-1 && yworld<=1){
			float t_min = 2;
			for (int i=0; i<spherical_cells.size(); i++){

				Eigen::Vector4f c = rotation * Eigen::Vector4f(spherical_cell_centers_unnormed[i](0), spherical_cell_centers_unnormed[i](1), spherical_cell_centers_unnormed[i](2), 1.0);
				// std::cout<<c<<std::endl;
				Eigen::Vector3f p = Eigen::Vector3f(xworld, yworld, sqrt(1 - pow(xworld, 2) - pow(yworld, 2)));
				Eigen::Vector3f ct = Eigen::Vector3f(c(0), c(1), c(2));
				// std::cout<<ct<<std::endl;
				float t = (p - ct).norm();
				if (t<t_min){
					t_min = t;
					cursor_position_cell = i;
				}
			}

		}else{
			cursor_position_cell = -1;
		}

	}


}

void key_callback(GLFWwindow* window, int key, int scancode, int action, int mods)
{
    switch (key)
    {	
    	case GLFW_KEY_R:

    	// reset the game
    		if (action==GLFW_RELEASE){

    			// clear normal grid cells
				cells.clear();
				cell_centers.clear();
				cell_center_squares.clear();
				neighbors.clear();
				mine_cells.clear();
				neighbor_mines.clear();
				cell_colors.clear();
				cell_complementary_colors.clear();
				number_centers_for_each_cell.clear();
				each_number_size_for_each_cell.clear();

				spherical_cells.clear(); 
				spherical_cell_centers.clear();
				spherical_cell_centers_unnormed.clear();
				spherical_cell_triangles.clear();
				sphere_center_squares.clear();
				spherical_mine_cells.clear();
				spherical_neighbor_mines.clear();
				spherical_layer_maps.clear();




				// clear game related parameters
				left_clicked_cell = -1;
				right_clicked_cell = -1;
				left_button_down = 0;
				right_button_down = 0;
				opened_cell.clear();
				game_over = false;
				flagged_cells.clear();
				cursor_position_cell = -1;
				win_game = false;

				t_start = std::chrono::high_resolution_clock::now();
				finish_time = 0;

				rx = 0;
				ry = 0;
				choice = 'g';
				x_position_last = 0;
				y_position_last = 0;

				rotation.resize(4, 4);
				rotation<<
				1.0, 0.0, 0.0, 0.0,
				0.0, 1.0, 0.0, 0.0,
				0.0, 0.0, 1.0, 0.0,
				0.0, 0.0, 0.0, 1.0;


				command = "python ../src/voronoi.py ";	
				sphere_command = "python ../src/spherical_voronoi.py ";

				std::cout<<"Please choose between sphere or square grid. Enter g for square grid and anything else for sphere."<<std::endl;

				std::cin>>choice;
				std::cout<<"Please enter number of cells"<<std::endl;
				
				std::cin>>num_of_cells;
				

				while (num_of_cells<=3){
					std::cout<<"The number of cells must be greater than 3, please enter another number"<<std::endl;
					std::cin>>num_of_cells;
				}
				
				if (choice !='g'){
					sphere_command+=std::to_string(num_of_cells);
					system(sphere_command.c_str());
				}else{
					command+= std::to_string(num_of_cells);
					system(command.c_str());
				}




				// ------------------------------------------------------------------------------------------------------------------------------------------------------------------
				// load voronoi
				
				
				std::cout<<"Please enter the number mines in this game"<<std::endl;
				std::cin>>number_of_mines;
				while(number_of_mines > 0.4 * num_of_cells || number_of_mines<=0){
					std::cout<<"The number of mines must be a positive number and cannot be greater than 40\% number of cells, please enter another number"<<std::endl;
					std::cin>>number_of_mines;
				}

				
				if (choice!='g'){
					// std::cout<<"Number of cells: "<<num_of_cells<<std::endl;
					load_spherical_voronoi(sphere_file,
			                            number_of_mines,
			                            spherical_cells, 
			                            spherical_cell_triangles,
			                            spherical_cell_centers,
			                            spherical_cell_centers_unnormed,
			                            sphere_center_squares,
			                            neighbors,
			                            spherical_mine_cells,
			                            spherical_neighbor_mines,
			                            spherical_layer_maps);
					// std::cout<<"1111111111111111111111111111111111111111"<<std::endl;
				}else{
					load_voronoi(file, 
								  number_of_mines, 
								  cells, 
								  cell_centers,
								  cell_center_squares,
								  neighbors, 
								  mine_cells, 
								  neighbor_mines,
								  square_cell_maps);

					square_cell_partition = (int)floor(cbrt(cell_centers.size()));
				}
				



				for (int i=0; i<num_of_cells; i++){
					cell_colors.push_back(Eigen::Vector3f(RandomFloat(0.5, 1.0), RandomFloat(0.5, 1.0), RandomFloat(0.5, 1.0)));
					cell_complementary_colors.push_back(complementary_color(cell_colors[i]));
				}



				if (choice=='g'){
					for (int i=0; i<num_of_cells; i++){
						std::vector<Eigen::Vector3f> number_centers;
						std::vector<float> each_number_size;
						int mine_number = neighbor_mines[i];

						getNumberCenterSize(cell_center_squares[i], 
												 mine_number, 
												 cell_centers[i], 
												 number_centers, 
												 each_number_size, 
												 all_numbers_width);
						number_centers_for_each_cell.push_back(number_centers);
						each_number_size_for_each_cell.push_back(each_number_size);

					}	
				}else{
					for (int i=0; i<num_of_cells; i++){
						std::vector<Eigen::Vector3f> number_centers;
						std::vector<float> each_number_size;
						int mine_number = spherical_neighbor_mines[i];

						getNumberCenterSize(sphere_center_squares[i], 
												 mine_number, 
												 Eigen::Vector3f(0, 0, 0), 
												 number_centers, 
												 each_number_size, 
												 all_numbers_width);
						number_centers_for_each_cell.push_back(number_centers);
						each_number_size_for_each_cell.push_back(each_number_size);

					}
				}
				// std::cout<<"222222222222222222222222222222222222222222222222"<<std::endl;


				break;	
    		}


    	default:
    		break;
    }



}







// The MAIN function, from here we start our application and run the Game loop
int main()
{

	// ------------------------------------------------------------------------------------------------------------------------------------------------------------------
	// initialization

	// Init GLFW
	// Initialize the library
	if (!glfwInit())
		return -1;

	// Activate supersampling
	glfwWindowHint(GLFW_SAMPLES, 8);
	glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
	glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
	glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
	glfwWindowHint(GLFW_RESIZABLE, GL_TRUE);

	// On apple we have to load a core profile with forward compatibility
	#ifdef __APPLE__
		glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
		glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE);
	#endif



	GLFWwindow* window = glfwCreateWindow(WIDTH, HEIGHT, "voronoi", nullptr, nullptr); // Windowed
	if (!window)
	{
		glfwTerminate();
		return -1;
	}

	glfwMakeContextCurrent(window);


	#ifndef __APPLE__
	  glewExperimental = true;
	  GLenum err = glewInit();
	  if(GLEW_OK != err)
	  {
		/* Problem: glewInit failed, something is seriously wrong. */
	   fprintf(stderr, "Error: %s\n", glewGetErrorString(err));
	  }
	  glGetError(); // pull and savely ignonre unhandled errors like GL_INVALID_ENUM
	  fprintf(stdout, "Status: Using GLEW %s\n", glewGetString(GLEW_VERSION));
	#endif



	// Initialize GLEW to setup the OpenGL Function pointers
	glewExperimental = GL_TRUE;



	// Set OpenGL options
	glEnable(GL_CULL_FACE);
	glEnable(GL_BLEND);
	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

	// Compile and setup the shader
	Shader shader("../src/shaders/text.vs", "../src/shaders/text.frag");



	// FreeType
	FT_Library ft;
	// All functions return a value different than 0 whenever an error occurred
	if (FT_Init_FreeType(&ft))
		std::cout << "ERROR::FREETYPE: Could not init FreeType Library" << std::endl;

	// Load font as face
	FT_Face face;
	if (FT_New_Face(ft, "../font/digital-7.ttf", 0, &face))
		std::cout << "ERROR::FREETYPE: Failed to load font" << std::endl;

	// Set size to load glyphs as
	FT_Set_Pixel_Sizes(face, 0, 48);

	// Disable byte-alignment restriction
	glPixelStorei(GL_UNPACK_ALIGNMENT, 1); 

	// Load first 128 characters of ASCII set
	for (GLubyte c = 0; c < 128; c++)
	{
		// Load character glyph 
		if (FT_Load_Char(face, c, FT_LOAD_RENDER))
		{
			std::cout << "ERROR::FREETYTPE: Failed to load Glyph" << std::endl;
			continue;
		}
		// Generate texture
		GLuint texture;
		glGenTextures(1, &texture);
		glBindTexture(GL_TEXTURE_2D, texture);
		glTexImage2D(
			GL_TEXTURE_2D,
			0,
			GL_RED,
			face->glyph->bitmap.width,
			face->glyph->bitmap.rows,
			0,
			GL_RED,
			GL_UNSIGNED_BYTE,
			face->glyph->bitmap.buffer
		);
		// Set texture options
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
		// Now store character for later use
		Character character = {
			texture,
			glm::ivec2(face->glyph->bitmap.width, face->glyph->bitmap.rows),
			glm::ivec2(face->glyph->bitmap_left, face->glyph->bitmap_top),
			face->glyph->advance.x
		};
		Characters.insert(std::pair<GLchar, Character>(c, character));
	}
	glBindTexture(GL_TEXTURE_2D, 0);
	// Destroy FreeType once we're finished
	FT_Done_Face(face);
	FT_Done_FreeType(ft);

	
	// Configure VAO/VBO for texture quads
	glGenVertexArrays(1, &VAO);
	glGenBuffers(1, &VBO);
	glBindVertexArray(VAO);
	glBindBuffer(GL_ARRAY_BUFFER, VBO);
	glBufferData(GL_ARRAY_BUFFER, sizeof(GLfloat) * 6 * 4, NULL, GL_DYNAMIC_DRAW);



	glEnableVertexAttribArray(0);
	glVertexAttribPointer(0, 4, GL_FLOAT, GL_FALSE, 4 * sizeof(GLfloat), 0);
	glBindBuffer(GL_ARRAY_BUFFER, 0);
	glBindVertexArray(0);
	

	VertexArrayObject vao;
	vao.init();
	vao.bind();

	// ------------------------------------------------------------------------------------------------------------------------------------------------------------------

	// ------------------------------------------------------------------------------------------------------------------------------------------------------------------
	// load numbers and draw numbers


	load_all_number(all_numbers_edges, all_numbers_width);

	Program program;
	const GLchar* number_vertex_shader =
				"#version 150 core\n"
					"in vec3 position;"
					"uniform vec3 center;"
					"uniform vec3 frag_color;"
					"out vec3 f_color;"
					"uniform mat4 view;"
					"uniform mat4 proj;"
					"uniform mat4 camera;"
					"uniform mat4 rot;"
					"void main()"
					"{"
					"    gl_Position = view*proj*camera*rot*(vec4((position + center), 1.0));"
					"    f_color = frag_color;"
					"}";
	const GLchar* number_fragment_shader =
			"#version 150 core\n"
					"in vec3 f_color;"
					"out vec4 outColor;"
					"void main()"
					"{"
					"    outColor = vec4(f_color, 1.0);"
					"}";

	// ------------------------------------------------------------------------------------------------------------------------------------------------------------------
	// load mines and flags
	Eigen::Matrix3Xf mine_edges;
	load_shape("../objects/mine.off", mine_edges);

	Eigen::Matrix3Xf flag_edges;
	load_shape("../objects/flag.off", flag_edges);




	// ------------------------------------------------------------------------------------------------------------------------------------------------------------------
	// ------------------------------------------------------------------------------------------------------------------------------------------------------------------
	// get the number of cells and store it in voronoi_data

	std::cout<<"Please choose between sphere or square grid. Enter g for square grid and anything else for sphere."<<std::endl;

	std::cin>>choice;



	std::cout<<"Please enter number of cells"<<std::endl;
	
	std::cin>>num_of_cells;

	while (num_of_cells<=3){
		std::cout<<"The number of cells must be greater than 3, please enter another number"<<std::endl;
		std::cin>>num_of_cells;
	}
	
	if (choice !='g'){
		sphere_command+=std::to_string(num_of_cells);
		system(sphere_command.c_str());
	}else{
		command+= std::to_string(num_of_cells);
		system(command.c_str());
	}




	// ------------------------------------------------------------------------------------------------------------------------------------------------------------------
	// load voronoi
	
	
	std::cout<<"Please enter the number mines in this game"<<std::endl;
	std::cin>>number_of_mines;
	while(number_of_mines > 0.4 * num_of_cells || number_of_mines<=0){
		std::cout<<"The number of mines must be a positive number and cannot be greater than 40\% number of cells, please enter another number"<<std::endl;
		std::cin>>number_of_mines;
	}
	

	if (choice!='g'){
		load_spherical_voronoi(sphere_file,
                            number_of_mines,
                            spherical_cells, 
                            spherical_cell_triangles,
                            spherical_cell_centers,
                            spherical_cell_centers_unnormed,
                            sphere_center_squares,
                            neighbors,
                            spherical_mine_cells,
                            spherical_neighbor_mines,
                            spherical_layer_maps);
	}else{
		load_voronoi(file, 
					  number_of_mines, 
					  cells, 
					  cell_centers,
					  cell_center_squares,
					  neighbors, 
					  mine_cells, 
					  neighbor_mines,
					  square_cell_maps);

		square_cell_partition = (int)floor(cbrt(cell_centers.size()));
	}
	



	for (int i=0; i<num_of_cells; i++){
		cell_colors.push_back(Eigen::Vector3f(RandomFloat(0.5, 1.0), RandomFloat(0.5, 1.0), RandomFloat(0.5, 1.0)));
		cell_complementary_colors.push_back(complementary_color(cell_colors[i]));
	}



	if (choice=='g'){
		for (int i=0; i<num_of_cells; i++){
			std::vector<Eigen::Vector3f> number_centers;
			std::vector<float> each_number_size;
			int mine_number = neighbor_mines[i];

			getNumberCenterSize(cell_center_squares[i], 
									 mine_number, 
									 cell_centers[i], 
									 number_centers, 
									 each_number_size, 
									 all_numbers_width);
			number_centers_for_each_cell.push_back(number_centers);
			each_number_size_for_each_cell.push_back(each_number_size);

		}	
	}else{
		for (int i=0; i<num_of_cells; i++){
			std::vector<Eigen::Vector3f> number_centers;
			std::vector<float> each_number_size;
			int mine_number = spherical_neighbor_mines[i];

			getNumberCenterSize(sphere_center_squares[i], 
									 mine_number, 
									 Eigen::Vector3f(0, 0, 0), 
									 number_centers, 
									 each_number_size, 
									 all_numbers_width);
			number_centers_for_each_cell.push_back(number_centers);
			each_number_size_for_each_cell.push_back(each_number_size);

		}
	}

	// ------------------------------------------------------------------------------------------------------------------------------------------------------------------

	// ------------------------------------------------------------------------------------------------------------------------------------------------------------------
	
	glfwSetFramebufferSizeCallback(window, framebuffer_size_callback);






	glfwSetKeyCallback(window, key_callback);
	t_start = std::chrono::high_resolution_clock::now();
	rotation = rotation_scale_matrix(0.0, 0.0, 0.0, 1.0, Eigen::Vector3f(0.0, 0.0, 0.0));



	// glfwSetWindowSizeLimits(window, WIDTH, HEIGHT, GLFW_DONT_CARE, GLFW_DONT_CARE);
	while (!glfwWindowShouldClose(window))
	{

		if (game_over){
			glfwSetMouseButtonCallback(window, empty_mouse_button_callback);
			glfwSetCursorPosCallback(window, empty_mouse_cursor_position_callback);
		}else{
			if (choice == 'g'){
				glfwSetMouseButtonCallback(window, mouse_button_callback);
				glfwSetCursorPosCallback(window, mouse_cursor_position_callback);
			}else{
				glfwSetMouseButtonCallback(window, sphere_mouse_button_callback);
				glfwSetCursorPosCallback(window, sphere_mouse_cursor_position_callback);
			}

		}
		

		// ------------------------------------------------------------------------------------------------------------------------------------------------------------------
		// initialization and text drawing
		// Get size of the window



		int width, height;
		glfwGetWindowSize(window, &width, &height);

		// if (width < WIDTH || height < HEIGHT){
		// 	glfwSetWindowSize(window, WIDTH, HEIGHT);
		// }

		glm::mat4 projection = glm::ortho(0.0f, static_cast<GLfloat>(width), 0.0f, static_cast<GLfloat>(height));
		shader.Use();
		glUniformMatrix4fv(glGetUniformLocation(shader.Program, "projection"), 1, GL_FALSE, glm::value_ptr(projection));

		// Check and call events
		glfwPollEvents();
		float aspect_ratio = float(height)/float(width); // corresponds to the necessary width scaling
		view <<
		aspect_ratio,0, 0, 0,
		0,           1, 0, 0,
		0,           0, 1, 0,
		0,           0, 0, 1;

		
		// Clear the colorbuffer
		
		glClearColor(background(0), background(1), background(2), 1.0f);
		glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

		Eigen::Vector3f text_color = complementary_color(background);


		t_now = std::chrono::high_resolution_clock::now();
		int time = (int)std::chrono::duration_cast<std::chrono::duration<float>>(t_now - t_start).count();
		glDisable(GL_DEPTH_TEST);
		RenderText(shader, "Time passed: "+std::to_string(time), 1, height-30, 0.5, glm::vec3(text_color(0), text_color(1), text_color(2)));
		RenderText(shader, "Flags left: "+std::to_string(number_of_mines - flagged_cells.size()), 1, height-60, 0.5, glm::vec3(text_color(0), text_color(1), text_color(2)));
		if (game_over){
			if (win_game){
				RenderText(shader, "You won using:", 1, height-90, 0.5, glm::vec3(text_color(0), text_color(1), text_color(2)));
				RenderText(shader, std::to_string(finish_time) + " seconds", 1, height-120, 0.5, glm::vec3(text_color(0), text_color(1), text_color(2)));
			}else{
				RenderText(shader, "Lost!", 1, height-90, 0.5, glm::vec3(text_color(0), text_color(1), text_color(2)));
			}
		}

		// ------------------------------------------------------------------------------------------------------------------------------------------------------------------
		// draw items
		glEnable(GL_DEPTH_TEST);

		camera_matrix = return_camera_matrix(return_rotation_matrix(0.0, 0.0, 0.0)*camera_position, return_rotation_matrix(0.0, 0.0, 0.0)*camera_look_up_endpoint, 0.0, 0.0);
		projection_matrix = return_orth_proj_matrix(1.0, camera_position);


		program.init(number_vertex_shader,number_fragment_shader,"outColor");
		vao.bind();
		program.bind();
		glUniformMatrix4fv(program.uniform("view"), 1, GL_FALSE, view.data());
		glUniformMatrix4fv(program.uniform("proj"), 1, GL_FALSE, projection_matrix.data());
		glUniformMatrix4fv(program.uniform("camera"), 1, GL_FALSE, camera_matrix.data());
		if (choice == 'g'){
			rotation.resize(4, 4);
			rotation<<
			1.0, 0.0, 0.0, 0.0,
			0.0, 1.0, 0.0, 0.0,
			0.0, 0.0, 1.0, 0.0,
			0.0, 0.0, 0.0, 1.0;
		}
		glUniformMatrix4fv(program.uniform("rot"), 1, GL_FALSE, rotation.data());



		glPolygonMode( GL_FRONT_AND_BACK, GL_FILL );
		VertexBufferObject vbo_cell;
		vbo_cell.init();
		vbo_cell.bind();


		if (choice!='g'){
			for (int i=0; i<spherical_cells.size(); i++){
				glUniformMatrix4fv(program.uniform("rot"), 1, GL_FALSE, rotation.data());
				glUniform3f(program.uniform("center"),0.0, 0.0, 0.0);


				if ((left_clicked_cell == i && left_button_down)){
					glUniform3f(program.uniform("frag_color"),cell_colors[i](0)/3 , cell_colors[i](1)/3 , cell_colors[i](2)/3);
				}else if (opened_cell.find(i)!=opened_cell.end()){
					// opened cell
					glUniform3f(program.uniform("frag_color"),cell_colors[i](0)/2 , cell_colors[i](1)/2 , cell_colors[i](2)/2);
				}else if (cursor_position_cell!=-1 && neighbors[cursor_position_cell].find(i) != neighbors[cursor_position_cell].end()){
					// neighbor cell
					glUniform3f(program.uniform("frag_color"),cell_colors[i](0) +0.2, cell_colors[i](1) +0.2, cell_colors[i](2) +0.2);
				}else if(cursor_position_cell == i){
					glUniform3f(program.uniform("frag_color"),0.5, 0.5, 0.5);
				}else{
					glUniform3f(program.uniform("frag_color"),cell_colors[i](0) , cell_colors[i](1) , cell_colors[i](2));
				}
				

				vbo_cell.update(spherical_cell_triangles[i]);
				program.bindVertexAttribArray("position", vbo_cell);
				glDrawArrays(GL_TRIANGLES, 0, spherical_cell_triangles[i].cols());


				

				if (spherical_mine_cells.find(i)!=spherical_mine_cells.end() && game_over){
					glUniform3f(program.uniform("frag_color"), 58.0/255.0, 50.0/255.0, 28.0/255.0);
					glUniform3f(program.uniform("center"),spherical_cell_centers[i](0), spherical_cell_centers[i](1), spherical_cell_centers[i](2));
					vbo_cell.update(mine_edges * each_number_size_for_each_cell[i][0] * 1.5);
					glDrawArrays(GL_TRIANGLES, 0, mine_edges.cols());
				}
				
				if (flagged_cells.find(i)!=flagged_cells.end()){
					glUniform3f(program.uniform("frag_color"), 165.0/255.0, 5.0/255.0, 29.0/255.0);
					glUniform3f(program.uniform("center"),1.1*spherical_cell_centers[i](0), 1.1*spherical_cell_centers[i](1), 1.1*spherical_cell_centers[i](2));
					vbo_cell.update(flag_edges * each_number_size_for_each_cell[i][0] * 2);
					glDrawArrays(GL_TRIANGLES, 0, flag_edges.cols());
				}
				Eigen::Matrix4Xf rot_tmp;
				rot_tmp.resize(4, 4);
				rot_tmp<<
				1.0, 0.0, 0.0, 0.0,
				0.0, 1.0, 0.0, 0.0,
				0.0, 0.0, 1.0, 0.0,
				0.0, 0.0, 0.0, 1.0;

				int mine_number = spherical_neighbor_mines[i];
				// glUniform3f(program.uniform("frag_color"), cell_complementary_colors[i](0) , cell_complementary_colors[i](1) , cell_complementary_colors[i](2));

				if (mine_number>0 && mine_cells.find(i) == mine_cells.end() && opened_cell.find(i)!=opened_cell.end()){
					if (mine_cells.find(i)!=mine_cells.end()){
						mine_number = 0;
					}
					glUniform3f(program.uniform("frag_color"), cell_complementary_colors[i](0) , cell_complementary_colors[i](1) , cell_complementary_colors[i](2));
					// draw numbers
					for (int j=0; j<std::to_string(mine_number).length(); j++){
						int b = (int)((std::to_string(mine_number))[j]) - 48;
						vbo_cell.update(all_numbers_edges[b] * each_number_size_for_each_cell[i][j] *0.9);
						program.bindVertexAttribArray("position", vbo_cell);


						glUniformMatrix4fv(program.uniform("rot"), 1, GL_FALSE, rot_tmp.data());
						Eigen::Vector4f mapped_location = (1+0.05)*rotation*Eigen::Vector4f(spherical_cell_centers[i](0), spherical_cell_centers[i](1), spherical_cell_centers[i](2), 1.0);

						glUniform3f(program.uniform("center"), number_centers_for_each_cell[i][j](0)+ mapped_location(0), number_centers_for_each_cell[i][j](1)+ mapped_location(1), 0+ mapped_location(2));
						// glUniform3f(program.uniform("center"), number_centers_for_each_cell[i][j](0), number_centers_for_each_cell[i][j](1), number_centers_for_each_cell[i][j](2));
						glDrawArrays(GL_TRIANGLES, 0, all_numbers_edges[b].cols());

					}
				}



				// std::cout<<"cell number :"<<i<<std::endl;
			}



	

		}else{
			for (int i=0; i<cells.size(); i++){
				vbo_cell.update(cells[i]);
				// std::cout<<i<<std::endl;
				// std::cout<<cells[i]<<std::endl<<std::endl;


				// draw the cell
				program.bindVertexAttribArray("position", vbo_cell);
				glUniform3f(program.uniform("center"),0.0, 0.0, 0.0);

				// mouse pressed cell
				if ((left_clicked_cell == i && left_button_down)){
					glUniform3f(program.uniform("frag_color"),cell_colors[i](0)/3 , cell_colors[i](1)/3 , cell_colors[i](2)/3);
				}else if (opened_cell.find(i)!=opened_cell.end()){
					// opened cell
					glUniform3f(program.uniform("frag_color"),cell_colors[i](0)/2 , cell_colors[i](1)/2 , cell_colors[i](2)/2);
				}else if (cursor_position_cell!=-1 && neighbors[cursor_position_cell].find(i) != neighbors[cursor_position_cell].end()){
					// neighbor cell
					glUniform3f(program.uniform("frag_color"),cell_colors[i](0) +0.2, cell_colors[i](1) +0.2, cell_colors[i](2) +0.2);
				}else if(cursor_position_cell == i){
					glUniform3f(program.uniform("frag_color"),0.5, 0.5, 0.5);
				}else{
					glUniform3f(program.uniform("frag_color"),cell_colors[i](0) , cell_colors[i](1) , cell_colors[i](2));
				}
				
				glDrawArrays(GL_TRIANGLE_FAN, 0, cells[i].cols());

				// draw the cell boundary
				glUniform3f(program.uniform("frag_color"),cell_complementary_colors[i](0) , cell_complementary_colors[i](1) , cell_complementary_colors[i](2));
				glDrawArrays(GL_LINE_LOOP, 0, cells[i].cols());


				// check mine counts for each cell
				int mine_number = neighbor_mines[i];
				// && mine_cells.find(i) == mine_cells.end()
				if (mine_number>0 && mine_cells.find(i) == mine_cells.end() && opened_cell.find(i)!=opened_cell.end()){

					if (mine_cells.find(i)!=mine_cells.end()){
						mine_number = 0;
					}

					glUniform3f(program.uniform("frag_color"), cell_complementary_colors[i](0) , cell_complementary_colors[i](1) , cell_complementary_colors[i](2));


					// draw numbers
					for (int j=0; j<std::to_string(mine_number).length(); j++){
						glUniform3f(program.uniform("center"),number_centers_for_each_cell[i][j](0), number_centers_for_each_cell[i][j](1), number_centers_for_each_cell[i][j](2));
						// std::cout<<std::to_string(number)[j];
						int b = (int)((std::to_string(mine_number))[j]) - 48;
						vbo_cell.update(all_numbers_edges[b] * each_number_size_for_each_cell[i][j]);
						glDrawArrays(GL_TRIANGLES, 0, all_numbers_edges[b].cols());

					}
				}

				if (mine_cells.find(i)!=mine_cells.end() && game_over ){
					glUniform3f(program.uniform("frag_color"), 58.0/255.0, 50.0/255.0, 28.0/255.0);
					glUniform3f(program.uniform("center"),cell_centers[i](0), cell_centers[i](1), cell_centers[i](2));
					vbo_cell.update(mine_edges * each_number_size_for_each_cell[i][0] * 2);
					glDrawArrays(GL_TRIANGLES, 0, mine_edges.cols());
				}

				if (flagged_cells.find(i)!=flagged_cells.end()){
					glUniform3f(program.uniform("frag_color"), 165.0/255.0, 5.0/255.0, 29.0/255.0);
					glUniform3f(program.uniform("center"),cell_centers[i](0), cell_centers[i](1), cell_centers[i](2)+1.0);
					vbo_cell.update(flag_edges * each_number_size_for_each_cell[i][0] * 2);
					glDrawArrays(GL_TRIANGLES, 0, flag_edges.cols());
				}

			}	
		}



	   
		// Swap the buffers
		glfwSwapBuffers(window);

	}

	glfwTerminate();
	return 0;
}

void RenderText(Shader &shader, std::string text, GLfloat x, GLfloat y, GLfloat scale, glm::vec3 color)
{
	// Activate corresponding render state	
	shader.Use();
	glUniform3f(glGetUniformLocation(shader.Program, "textColor"), color.x, color.y, color.z);
	glActiveTexture(GL_TEXTURE0);
	glBindVertexArray(VAO);

	// Iterate through all characters
	std::string::const_iterator c;
	for (c = text.begin(); c != text.end(); c++) 
	{
		Character ch = Characters[*c];

		GLfloat xpos = x + ch.Bearing.x * scale;
		GLfloat ypos = y - (ch.Size.y - ch.Bearing.y) * scale;

		GLfloat w = ch.Size.x * scale;
		GLfloat h = ch.Size.y * scale;
		// Update VBO for each character
		GLfloat vertices[6][4] = {
			{ xpos,     ypos + h,   0.0, 0.0 },            
			{ xpos,     ypos,       0.0, 1.0 },
			{ xpos + w, ypos,       1.0, 1.0 },

			{ xpos,     ypos + h,   0.0, 0.0 },
			{ xpos + w, ypos,       1.0, 1.0 },
			{ xpos + w, ypos + h,   1.0, 0.0 }           
		};
		// Render glyph texture over quad
		glBindTexture(GL_TEXTURE_2D, ch.TextureID);
		// Update content of VBO memory
		glBindBuffer(GL_ARRAY_BUFFER, VBO);
		glBufferSubData(GL_ARRAY_BUFFER, 0, sizeof(vertices), vertices); // Be sure to use glBufferSubData and not glBufferData

		glBindBuffer(GL_ARRAY_BUFFER, 0);
		// Render quad
		glDrawArrays(GL_TRIANGLES, 0, 6);
		// Now advance cursors for next glyph (note that advance is number of 1/64 pixels)
		x += (ch.Advance >> 6) * scale; // Bitshift by 6 to get value in pixels (2^6 = 64 (divide amount of 1/64th pixels by 64 to get amount of pixels))
	}
	glBindVertexArray(0);
	glBindTexture(GL_TEXTURE_2D, 0);
}