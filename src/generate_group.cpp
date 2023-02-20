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

#include <unordered_set>
#include <vector>

#include <cmath>
#include <cstdlib>
#include <stdexcept>

#include "triangle_group.h"

void add_to_basis(G &candidate,  unordered_set<G, GHash> &unordered_basis,  vector<G> &next_generation){
  if (unordered_basis.find(candidate) == unordered_basis.end())
  {
    unordered_basis.insert(candidate);
    next_generation.push_back(candidate);
  }
}

int main()
{

  using namespace std;

  // -- reading group secifications -----------------------------------

  const char* env = std::getenv("HYPERBOLIC_DIR");
  string project_dir = env;

  ifstream input_file(project_dir+"/group_specs.inp");
  char char_buffer;
  int int_buffer;
  int p, q, N;
  int modulo;

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


  TriangleGroup T = TriangleGroup(p,q);

  G A = T.A;
  G B = T.B;
  G E = G(T.reduction);

  if(periodic_boundary)
  {
    A = A % modulo;
    B = B % modulo;
  }
  
  E.identity();

  // -- found elements are stored here
  std::unordered_set<G, GHash> unordered_basis;
  // vector<G> basis;

  // -- keep track of the iterative process
  vector<G> prev_generation;
  vector<G> next_generation;

  // -- initialize identity element
  // unordered_basis.insert(E);
  prev_generation.push_back(E);

  int generation=0;

  G generator = G(T.reduction);
  G candidate = G(T.reduction);
  int order;
  bool checked_twice = false;
  while(true)
  {
    if(generation%2==0)
    {
      generator = A;
      order = p;
    }
    else
    {
      generator = B;
      order = q;
    }

    if(periodic_boundary) generator = generator % modulo;

    for(G elem: prev_generation)
    {
      candidate = elem;

      for(int i=0; i<order; i++)
      {
        candidate = candidate * generator;

        if(periodic_boundary) candidate = candidate % modulo;

        add_to_basis(candidate, unordered_basis, next_generation);
      }
    }

    generation += 1;

    cout << generation << "\t | \t" << unordered_basis.size() << endl;
    
    if(periodic_boundary)
    {
      if (next_generation.size() == 0)
      {
	if(checked_twice) break;
	checked_twice= true;
      }
    }
    else
    {
      if(generation>N) break;
    }

    if(next_generation.size() != 0)
    {
      prev_generation = next_generation;
      next_generation.clear();
    }
  }

  string output_file_name;
  if(periodic_boundary)
  {
    output_file_name = project_dir+"/"+to_string(p)+","+to_string(q)+"_modulo_"+to_string(modulo)+".words";
  }
  else
  {
    output_file_name = project_dir+"/"+to_string(p)+","+to_string(q)+"_open_"+to_string(N+1)+".words";
  }
  
  ofstream output_file; 
  output_file.open(output_file_name, ofstream::out | ofstream::trunc);

  for(G elem: unordered_basis)
  {
    if(elem.word=="")
    {
      output_file << "E"<< endl;
    }
    else
    {
      output_file << elem.word << endl;
    }
  }

  output_file.close();

  return 0;
}
