#include <armadillo>
#include <iostream>
#include <fstream>
#include "Schrodinger_init.cpp"


using namespace std;




int main(){

  //Defining the varaibles from the input file
  double h, dt, T, x_c, sigma_x, p_x, y_c, sigma_y, p_y, v0;
  string wall, file_name, input_file_name;
  int nr_slits;

  //Lets the user give a input file and name the output file
  cout << "Input file name: ";
  cin >> input_file_name;
  cout << "Output file name: ";
  cin >> file_name;

  //Open the input file
  ifstream infile;
  infile.open(input_file_name);


  if(infile.fail()){
    cout << "Invalid file name"; //Check if the input name is a valid file
  }
  else{
    //Here we read the file and give all the variables their correct value.
    //The input file need to be of the same shape as the examples given in the github repository. Else, this will either give an error or wrong calculations.
    string s;
    int c = 1;
    while(infile>>s){
      if(c == 3){
        cout << "This is h:" << s<< "\n";
        h = std::stod(s);
      }
      else if(c == 6){
        dt = stod(s);
      }
      else if(c == 9){
        T = stod(s);
      }
      else if(c == 12){
        x_c = stod(s);
      }
      else if(c == 15){
        sigma_x = stod(s);
      }
      else if(c == 18){
        p_x = stod(s);
      }
      else if(c == 21){
        y_c = stod(s);
      }
      else if(c == 24){
        sigma_y = stod(s);
      }
      else if(c == 27){
        p_y = stod(s);
      }
      else if(c == 30){
        wall = s;
      }
      else if(c == 33){
        nr_slits = stoi(s);
      }
      c++;
    }
  }

  //Check if the user wants the wall on or off. If wall is on we give v0 a huge value, else we set it to be 0.
  if(wall == "on"){
    v0 = 1.*pow(10,10);
  }else{
    v0 = 0;
  }



  //Calculate M and M-2 for the given h
  int M = 1/h + 1;
  int m = M-2;

  //Create wall
  arma::sp_mat V(m,m);
  create_wall(nr_slits, v0, h, V);


  //Initialize the A and B matrix before filling it with the correct values
  arma::sp_mat imag(pow(m,2), pow(m,2));
  arma::sp_cx_mat A(arma::speye (pow(m,2), pow(m,2)), imag);
  arma::sp_cx_mat B(arma::speye (pow(m,2), pow(m,2)), imag);

  //Fill A and B with the correct values
  fill_A_B(A, B, M, h, dt, V);

  //Initialize u, give it the correct values for t=0 and normalize it.
  arma::cx_vec u(pow(m,2));
  init_u(u, M, h, x_c, y_c, sigma_x, sigma_y, p_x, p_y);
  normalize_vec(u);

  //Write the u values to file when t=0
  std::ofstream myfile;
  myfile.open(file_name);
  for(int j=0; j<M-1; j++){
    myfile << std::scientific << 0 << " " << std::scientific<< 0 << " ";
  }
  for(int i=0; i<u.size(); i++){
    if(i%m==0){
      myfile << std::scientific << 0 << " " << std::scientific<< 0 << " ";
      myfile << std::scientific << 0 << " " << std::scientific<< 0 << " ";
    }
    myfile << std::scientific << real(u(i)) << " " << std::scientific<< std::imag(u(i)) << " ";
  }
  for(int j=0; j<M; j++){
    myfile << std::scientific << 0 << " " << std::scientific<< 0 << " ";
  }
  myfile << std::scientific << 0 << " " << std::scientific<< 0 << " ";
  double prob;

  //Calculate the total probability for time t=0 and write to file
  prob = calculate_prob(u);
  myfile << std::scientific << 0 << " ";
  myfile << std::scientific << 1.-prob << "\n";


  //Loop from time t=dt to t=T
  arma::cx_vec b;
  for(double t=dt; t<=T; t += dt){
    cout << t << "\n"; //Print t to show progress

    //Solve the matrix equations to find the next timesetp
    b = B*u;
    u = spsolve(A,b);

    //Write the next timestep to file
    for(int j=0; j<M-1; j++){
      myfile << std::scientific << 0 << " "<< std::scientific << 0 << " ";
    }
    for(int i=0; i<u.size(); i++){
      if(i%m==0){
        myfile << std::scientific << 0 << " " << std::scientific<< 0 << " ";
        myfile << std::scientific << 0 << " " << std::scientific<< 0 << " ";
      }
      myfile << std::scientific << real(u(i)) << " " << std::scientific << std::imag(u(i)) << " ";
    }
    for(int j=0; j<M; j++){
      myfile << std::scientific << 0 << " " << std::scientific<< 0 << " ";
    }
    myfile << std::scientific << 0 << " " << std::scientific<< 0 << " ";

    //Calculate the total probability for time t and write to file
    prob = calculate_prob(u);
    myfile << std::scientific << t << " ";
    myfile << std::scientific << 1.-prob << "\n";
  }

}
