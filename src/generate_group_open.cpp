/*
 *  ./src/ring_extension.cpp
 *
 *  Author:  Fabian R. Lux
 *  Date:    12/01/2023
 *
 *  Constructs the ring extension Z[x] where x^d fullfills the relation
 *  x^d = r(0) x^0 + r(1) x^1 + ... r(d-1) x^(d-1).
 */
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>

#include <vector>

#include <stdexcept>

int main()
{

  using namespace std;

  // -- reading group secifications -----------------------------------
  ifstream input_file("group_specs.inp");
  char char_buffer;
  int int_buffer;
  int p, q, N;

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
    N = -N;

    cout << "Max. number of generations = " << N << endl;
  }
  else
  {
    cout << "Periodic boundary conditions are used." << endl;

    periodic_boundary = true;

    cout << "The matrix representation is treated modulo 2**N, with N=" << N << endl;
  }

  // -- reading ring secifications ------------------------------------

  string line;
  ifstream file("ring_reduction.inp");

  int n = 2 * p * q;

  cout << "Searching ring_reduction.inp for n=2pq=" << n << endl;

  vector<int> reduction;
  int buffer;
  int d = 0;

  bool listening = false;
  
  if (file.is_open())
  {
    while (getline(file, line))
    {
      //--start analysing after keyword is found
      if (line.find("BEGIN") != string::npos)
      {
        listening = true;
      }

      //--stop analysing after keyword is found
      if (line.find("END") != string::npos)
      {
        break;
      }

      if (listening)
      {
        istringstream iss(line);

        if (iss >> buffer && iss >> d && buffer == n)
        {
          while (iss >> buffer)
          {
            reduction.push_back(buffer);
          }
          //--no need to keep listening
          break;
        }
      }
    }
  }
  else
  {
    throw runtime_error("The file ring_reduction.inp was not found!");
  }
  file.close();

  cout << "Dimension of ring extension is d=" << d << endl;


  if(periodic_boundary)
  {
    // -- periodic boundary conditions
  }
  {
    // -- open boundary conditions
  }

  // cout << "The coefficients are:" << endl;

  // for (int c : reduction)
  // {
  //   cout << c << endl;
  // }

  return 0;
}
