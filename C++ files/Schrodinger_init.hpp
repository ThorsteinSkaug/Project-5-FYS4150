#ifndef __Schrodinger_functions_hpp__
#define __Schrodinger_functions_hpp__
#include "Schrodinger_init.hpp"
#include <armadillo>

using namespace std;



/*
Function for calculate the total sum of u*u for a vector u
Input variables:
- arma::cx_vec u: A complex vector u
Output:
- The total sum of u*u
*/
double calculate_prob(arma::cx_vec u);


/*
Function for finding the correct index on a flatten 2D array given the 2D indices i and j.
Input variables:
- int i: The row in the 2D array
- int j: The column in the 2D array
- int M: The number of columns in the full matrix
Output:
- The index of the flatten array
*/
int find_index(int i, int j, int M);


/*
Function for making the A and B matrix.
Input variables:
- arma::sp_cx_mat A: A complex armadillo matrix where the real part is an identity matrix
- arma::sp_cx_mat B: A complex armadillo matrix where the real part is an identity matrix
- arma::cx_double r: The r value to be put on the off-diagonal entries
- arma::cx_vec a: A complex vector storing the diagonal entries of the A matrix
- arma::cx_vec b: A complex vector storing the diagonal entries of the B matrix
Output:
- Update of the A and B matrix, filled with correct values
*/
void make_A_B(arma::sp_cx_mat& A, arma::sp_cx_mat& B, arma::cx_double r, arma::cx_vec a, arma::cx_vec b);


/*
Function for filling A and B with the correct values
Input varaibles:
- arma::sp_cx_mat A: A complex armadillo matrix where the real part is an identity matrix
- arma::sp_cx_mat B: A complex armadillo matrix where the real part is an identity matrix
- int M: The number of columns in the full matrix
- double h: Step size in x and y direction
- double dt: The time step
- arma::sp_mat V: The wall V
Output:
- The correctlly filled A and B
*/
void fill_A_B(arma::sp_cx_mat& A, arma::sp_cx_mat& B, int M, double h, double dt, arma::sp_mat V);



/*
Funtion for normalizing a complex vector
Input variables:
- arma::cx_vec vek: A complex vector
Output:
- Update of the vector vek such that is now normalized
*/
void normalize_vec(arma::cx_vec& vek);


/*
Funtion to initialize the start values of u
Input variables:
- arma::cx_vec u0: The complex vector u0 to be initilized
- int M: The number of columns in the full matrix
- double h: Step size in x and y direction
- double x_c: The coordinates of the centre of the initial wave packet
- double y_c: The coordinates of the centre of the initial wave packet
- double sigma_x: The initial widths of the wave packet
- double sigma_y: The initial widths of the wave packet
- double p_x: The wave packet momenta
- double p_y: The wave packet momenta
Output:
- The initialized u0
*/
void init_u(arma::cx_vec& u0, int M, double h, double x_c, double y_c, double sigma_x, double sigma_y, double p_x, double p_y);



/*
Function for creating the wall (V)
Input variables:
- int nr_slits: The number of slits in the wall
- double v0: A big variable which creates the wall
- double h: Step size in x and y direction
- arma::sp_mat V: A matrix V to simulate the wall
Output:
An updated version of V
*/
void create_wall(int nr_slits, double v0, double h, arma::sp_mat& V);



#endif
