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
#include <algorithm>

#include <cstdlib>
#include <stdexcept>

#include "triangle_group.h"

int main()
{
  using namespace std;

  // -- reading group secifications -----------------------------------

  const char* env = std::getenv("HYPERBOLIC_BUILD");
  string build_dir = env;

  ifstream input_file(build_dir+"/bin/group_specs.inp");
  char char_buffer;
  int int_buffer;
  int p, q, N;
  int modulo;

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

  TriangleGroup T=TriangleGroup(p, q);
  
  ifstream basis_file(basis_file_name); 
  string line;

  G word=G(T.reduction);

  vector<G> basis;

  // Construct unit element
  G A = T.A;
  G B = T.B;
  G AB = A*B;

  if(periodic_boundary){
    A = A % modulo;
    B = B % modulo;
    AB = AB % modulo;
  } 
  
  if(basis_file.is_open())
  {
    while(getline(basis_file, line))
      {
        word=G(T.reduction);
        word.identity();
        for(char& c : line) {
          if(c=='A')
          {
            word = word * A;
          }
          else if(c=='B')
          {
            word = word * B;
          }
          else if(c=='E')
          {
          }
        }
        if(periodic_boundary) word = word % modulo;
        basis.push_back(word);
      }
  }
  else
  {
    throw runtime_error("Basis file not found! It needs to be generated first using bin/generate_group");
  }

  cout << "Basis file loaded succesfully! Basis dimension: " << basis.size() << endl;

  cout << "Generating the right regular representation." << endl;

  vector<G> operators;
  
  operators.push_back(A);
  operators.push_back(B);
  operators.push_back(AB);

  int i,j;
  G action=G(T.reduction);

  string output_file_name;
  if(periodic_boundary)
  {
    output_file_name = build_dir+"/bin/{"+to_string(p)+","+to_string(q)+"}_modulo_"+to_string(modulo)+"_";
  }
  else
  {
    output_file_name = build_dir+"/bin/{"+to_string(p)+","+to_string(q)+"}_open_"+to_string(N+1)+"_";
  }

  for(G op: operators)
  {
    j=0;
    ofstream output_file; 
    output_file.open(output_file_name+op.word+".reg", ofstream::out | ofstream::trunc); 

    for(G b: basis)
    {
      action = b * op;

      if(periodic_boundary) action = action % modulo;

      auto it = find(basis.begin(), basis.end(), action);
  
      if (it != basis.end()) 
      {
        i = it - basis.begin();
        output_file << i << " " << j << endl;
      } 
      else
      {
        if(periodic_boundary)
        {
          // -- must not happen for periodic boundary conditions
          throw runtime_error("Right regular representation failed.");
        }     
      }
      j += 1;
    }
    output_file.close();
  }
  
  // output_file.open(output_file_name, oddddddddfstream::out | ofstream::trunc);

  // for(G elem: unordered_basis)
  // {
  //   output_file << elem.word << endl;
  // }

  // output_file.close();
  
  // Set-up the Laplacians for A,B and AB
  
  // Store the Laplacians to file

  return 0;
}
