/*
 *  ./src/tffg_generate_group.cpp
 *
 *  Author:  Fabian R. Lux
 *  Date:    8/1/2023
 *
 *  Generate the torsion free Fuchsian group
 */
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>

#include <unordered_map>
#include <vector>

#include <cmath>
#include <cstdlib>
#include <stdexcept>

#include "triangle_group.h"
#include "tffg_group.h"

void add_to_basis(int &count, G &candidate,   unordered_map<G,int,GHash> &unordered_basis,  vector<G> &next_generation){
  if (unordered_basis.find(candidate) == unordered_basis.end())
  {
    unordered_basis.insert( {candidate, count++});
    next_generation.push_back(candidate);
  }
}

int main()
{

  using namespace std;

  // -- reading group secifications -----------------------------------

  const char* env = std::getenv("HYPERBOLIC_DIR");
  string project_dir = env;

  ifstream input_file(project_dir+"/tffg_group_specs.inp");
  char char_buffer;
  int int_buffer;
  int g, N;
  int modulo;

  if (input_file.is_open())
  {
    while (input_file >> char_buffer >> int_buffer)
    {
      if (char_buffer == 'g')
      {
        g = int_buffer;
      }
      else if (char_buffer == 'N')
      {
        N = int_buffer;
      }
    }
  }
  else
  {
    throw runtime_error("The file tffg_group_spec.inp was not found!");
  }
  input_file.close();

  cout << "#------------------------------------------------#" << endl;
  cout << "# C++ Spectra of Hyperbolic Cayley Crystals v0.1 #" << endl;
  cout << "#------------------------------------------------#" << endl;
  cout << "The construction of the torsion-free Fuchsian group of genus " << endl;
  cout << "g=" << g << endl;

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

  TorsionFreeFuchsianGroup T = TorsionFreeFuchsianGroup(g);

  vector<G> generators;

  // G A = T.A;
  // G B = T.B;
  G E = G(T.reduction);
  E.identity();
  E.word = to_string(0) + " ";


  if(periodic_boundary)
  {
    for(G gamma: T.gamma)
    {
      generators.push_back(gamma % modulo);
    }
  }
  else
  {
    generators = T.gamma;
  }
  
  
  // -- found elements are stored here
  unordered_map<G,int,GHash> unordered_basis;
  vector<G> basis;

  // -- keep track of the iterative process
  vector<G> prev_generation;
  vector<G> next_generation;

  // -- initialize identity element
  // unordered_basis.insert(E);
  // prev_generation.push_back(E);


  int count = 0;
  int generation=0;

  // G generator = G(T.reduction);
  G candidate = G(T.reduction);

  add_to_basis(count, E, unordered_basis, prev_generation);

  cout << generation << "\t | \t" << unordered_basis.size() << endl;

  bool checked_twice = false;
  while(true)
  {
    
    for(G elem: prev_generation)
    {
      for(G gamma: generators)
      {
        candidate = elem;
        
        candidate = candidate * gamma;

        if(periodic_boundary) candidate = candidate % modulo;

        add_to_basis(count, candidate, unordered_basis, next_generation);
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
    output_file_name = project_dir+"/tffg_"+to_string(g)+"_modulo_"+to_string(modulo)+".words";
  }
  else
  {
    output_file_name = project_dir+"/tffg_"+to_string(g)+"_open_"+to_string(N+1)+".words";
  }
  
  ofstream output_file; 
  output_file.open(output_file_name, ofstream::out | ofstream::trunc);

  for(auto elem: unordered_basis)
  {
    output_file << to_string(elem.second) << " " ;
    if(elem.first.word=="")
    {
      output_file << "E"  << endl;
    }
    else
    {
      output_file << elem.first.word  << endl;
    }
  }

  output_file.close();

  return 0;
}
