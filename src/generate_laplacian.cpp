/*
 * ./src/generate_laplacian.cpp
 *
 * Author:  Fabian R. Lux
 * Date:    02/15/2023
 *
 * Read the basis file and construct the Laplacian.
 *
 */
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>

#include <vector>

#include <cstdlib>
#include <stdexcept>

#include "triangle_group.h"

int main()
{
  using namespace std;

  // Read the input files
  if (input_file.is_open())
  {
    while (input_file >> char_buffer >> int_buffer)
    {
      if (char_buffer == 'p')
      {
        p = int_buffer;
      }
      else if (char_buffer == 'q')
      {
        q = int_buffer;
      }
        else if (char_buffer == 'N')
      {
        N = int_buffer;
      }
    }
  }
  else
  {
   throw runtime_error("The file group_spec.inp was not found!");
  }
  input_file.close(); 
  
   cout << "#------------------------------------------------#" << endl;
  cout << "# C++ Spectra of Hyperbolic Cayley Crystals v0.1 #" << endl;
  cout << "#------------------------------------------------#" << endl;
  cout << "The construction of the proper triangle group Delta+(p,q,2) commences for:" << endl;
  cout << "p=" << p << endl;
  cout << "q=" << q << endl;

  bool periodic_boundary;
  if (N < 0)
  {
    cout << "Open boundary conditions are used." << endl;

    periodic_boundary = false;
    N = -N - 1;

    cout << "Max. number of generations = " << N+1 << endl;
  }
  else
  {
    cout << "Periodic boundary conditions are used." << endl;

    periodic_boundary = true;
    modulo = N;//pow(2,N);

    cout << "The matrix representation is treated modulo " << N << endl;
  }


  // Read the basis file 
  string basis_file_name;
  if(periodic_boundary)
  {
    basis_file_name = build_dir+"/bin/{"+to_string(p)+","+to_string(q)+"}_modulo_"+to_string(modulo)+".words";
  }
  else
  {
    basis_file_name = build_dir+"/bin/{"+to_string(p)+","+to_string(q)+"}_open_"+to_string(N+1)+".words";
  }
  
  ifstream basis_file(basis_file_name); 
  string line;
  while(getline(line,basis_file) )
  {
  }
  output_file.open(output_file_name, oddddddddfstream::out | ofstream::trunc);

  for(G elem: unordered_basis)
  {
    output_file << elem.word << endl;
  }

  output_file.close();
  
  // Set-up the Laplacians for A,B and AB
  
  // Store the Laplacians to file

  return 0;
}
