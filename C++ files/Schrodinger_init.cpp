#include "Schrodinger_init.hpp"
#include <math.h>



// For a clear explanation of all the functions and theirs input and outputs see the header file Schrodinger_init.hpp


double calculate_prob(arma::cx_vec u){
  double p=0.; //Variable for storing the sum
  for(int i=0; i<u.size(); i++){
    p += pow(real(u(i)),2) + pow(std::imag(u(i)),2); //Loop through the vector and add u*u to p
  }
  return p;
}



int find_index(int i, int j, int M){
  int m = M-2; //Find the number of internal points
  int begin = m*(j-1); //Get to the correct beginning of row
  return begin + (i -1); //The correct index
}



void make_A_B(arma::sp_cx_mat& A, arma::sp_cx_mat& B, arma::cx_double r, arma::cx_vec a, arma::cx_vec b){
  int m2 = a.size(); //The number of rows/columns in A and B
  int m = sqrt(m2); //The number of sub-matrices
  A.diag() = a; //Update the diagonal of A
  B.diag() = b; //Update the diagonal of B

  //Loop through the the m2-m upper/lower diagonal and add the r values
  for(int i=0; i<m2-m; i++){
    A(i, i+m) = -r;
    B(i, i+m) = r;
    A(i+m, i) = -r;
    B(i+m, i) = r;
  }

  //Loop through the upper and lower diagonal and add r values to the correct place
  for(int j=0; j<m2-1; j++){
    if((j+1)%m == 0){
      continue; //Drop adding r values
    }else{
      A(j,j+1) = -r;
      A(j+1,j) = -r;
      B(j,j+1) = r;
      B(j+1,j) = r;
    }
  }
}




void fill_A_B(arma::sp_cx_mat& A, arma::sp_cx_mat& B, int M, double h, double dt, arma::sp_mat V){
  arma::cx_double r (0,dt/(2*h*h)); //Tha r values
  arma::cx_vec a(pow(M-2,2)); //Initialize a
  arma::cx_vec b(pow(M-2,2)); //Initialize b

  //Loop through the inner matrix
  for(int i=0; i<M-2; i++){
    for(int j=0; j<M-2; j++){
      int k = find_index(i+1,j+1,M); //Find the correct index
      arma::cx_double a_k (1, 4*imag(r)+(dt/2)*V(i,j)); //Find the correct a values
      a(k) = a_k; //Add the correct a values to the a vector

      arma::cx_double b_k (1, -4*imag(r)-(dt/2)*V(i,j)); //Find the correct b values
      b(k) = b_k; //Add the correct b values to the b vector
    }
  }
  make_A_B(A, B, r, a, b); //Make the A and B matrices
}





void normalize_vec(arma::cx_vec& vek){
  double s = 0; //Variable storing the total sum
  for(int i=0; i<vek.size(); i++){
    s += pow(real(vek(i)),2) + pow(imag(vek(i)),2); //Loop through the vector and add to the total sum
  }
  for(int j=0; j<vek.size(); j++){
    vek(j) = vek(j)/sqrt(s); //Loop through the vector and divide by the square root of the total sum to normalize the values
  }
}



void init_u(arma::cx_vec& u0, int M, double h, double x_c, double y_c, double sigma_x, double sigma_y, double p_x, double p_y){

  //Loop through the u0 vector
  for(int i=0; i<sqrt(u0.size()); i++){
    double x = (i+1)*h; //Find the x coordinate
    for(int j=0; j<sqrt(u0.size()); j++){
      double y = (j+1)*h; //Find the y coordinate
      int k = find_index(i+1,j+1, M); //Find the corresponding flatten index
      double alpha = -pow(x-x_c,2)/(2*pow(sigma_x,2))-pow(y-y_c,2)/(2*pow(sigma_y,2));
      double beta = p_x*(x-x_c)+p_y*(y-y_c);
      double real = exp(alpha)*cos(beta);
      double ima = exp(alpha)*sin(beta);
      arma::cx_double uij (real, ima);
      u0(k) = uij; //Add the initial u values
    }
  }
}



void create_wall(int nr_slits, double v0, double h, arma::sp_mat& V){
  //Variables for storing the start and end of the slits
  double start_first_split = 0.;
  double end_first_split = 0.;
  //The second and thirs slit have been given default -1.5 values such that if the are not update they will not affect the next loop
  double start_second_split = -1.5;
  double end_second_split = -1.5;
  double start_third_split = -1.5;
  double end_third_split = -1.5;

  //Variables for storing where the wall in x-direction starts and ends
  double begin_x_wall = 0.5-(0.02/2);
  double end_x_wall = 0.5+(0.02/2);

  //Check how many slits and update the variables
  if(nr_slits == 1){
    start_first_split = 0.5-(0.05/2);
    end_first_split = 0.5+(0.05/2);
  }
  else if(nr_slits == 2){
    start_first_split = 0.5-(0.05/2)-0.05;
    end_first_split = 0.5-(0.05/2);

    start_second_split = 0.5+(0.05/2);
    end_second_split = 0.5+(0.05/2)+0.05;
  }else{
    start_second_split = 0.5-(0.05/2);
    end_second_split = 0.5+(0.05/2);

    start_first_split = start_second_split - 2*0.05;
    end_first_split = start_first_split + 0.05;

    start_third_split = end_second_split + 0.05;
    end_third_split = start_third_split + 0.05;
  }

  //Loop through the wall from the bottom to top in y direction
  for(double y=h; y<1; y+=h){
    //If we are between the start of the first split and end of the first split we jump the length of the split
    if(y >= start_first_split && end_first_split >= y){
      y += 0.05-h;

      //If 0.05 is not divisible by h we will jump to a point which we do not have, hence this if statement make sure we jump to the next possible point
      if(fmod(y,h) != 0){
        y += h-(fmod(y,h));
      }
    }
    //If we are between the start of the second split and end of the second split we jump the length of the split
    else if( y >= start_second_split && end_second_split >= y){
      y += 0.05-h;

      //If 0.05 is not divisible by h we will jump to a point which we do not have, hence this if statement make sure we jump to the next possible point
      if(fmod(y,h) != 0){
        y += h-(fmod(y,h));
      }
    }
    //If we are between the start of the third split and end of the third split we jump the length of the split
    else if( y >= start_third_split && end_third_split >= y){
      y += 0.05-h;

      //If 0.05 is not divisible by h we will jump to a point which we do not have, hence this if statement make sure we jump to the next possible point
      if(fmod(y,h) != 0){
        y += h-(fmod(y,h));
      }
    }
    //If we are at the wall we add v0 to the V matrix
    else{
      for(double x=begin_x_wall; x <= end_x_wall; x+=h){
        V(x/h - 1, y/h -1) = v0;
      }
    }
  }
}
