/*
 * ./src/tffg_generate_reg.cpp
 *
 * Author:  Fabian R. Lux
 * Date:    9/7/2023
 *
 * Generate the right-regular representation of the generators of the Fuchsian group
 * of genus g.
 *
 */
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>

#include <vector>
#include <unordered_map>
#include <algorithm>

#include <cstdlib>
#include <stdexcept>

#include "triangle_group.h"
#include "tffg_group.h"

#pragma omp declare reduction(vec_int_plus : std::vector<int> : \
        std::transform(omp_out.begin(), omp_out.end(), omp_in.begin(), omp_out.begin(), std::plus<int>())) \
        initializer(omp_priv = decltype(omp_orig)(omp_orig.size()))

int main()
{
  using namespace std;

  // -- reading group secifications -----------------------------------

  const char *env = std::getenv("HYPERBOLIC_DIR");
  string project_dir = env;

  ifstream input_file(project_dir + "/tffg_group_specs.inp");
  char char_buffer;
  int int_buffer;
  int g = 2;
  int N=2;
  int modulo=2;

  // Read the input files
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
    throw runtime_error("The file group_spec.inp was not found!");
  }
  input_file.close();

  cout << "#------------------------------------------------#" << endl;
  cout << "# C++ Spectra of Hyperbolic Cayley Crystals v0.1 #" << endl;
  cout << "#------------------------------------------------#" << endl;
  cout << "Reconstruction of the torsion free Fuchsian group of genus" << endl;
  cout << "g=" << g << endl;

  bool periodic_boundary = false;
  if (N < 0)
  {
    cout << "Open boundary conditions are used." << endl;

    periodic_boundary = false;
    N = -N - 1;

    cout << "Max. number of generations = " << N + 1 << endl;
  }
  else
  {
    cout << "Periodic boundary conditions are used." << endl;

    periodic_boundary = true;
    modulo = N; // pow(2,N);

    cout << "The matrix representation is treated modulo " << N << endl;
  }

  // Read the basis file
  string basis_file_name;
  if (periodic_boundary)
  {
    basis_file_name = project_dir + "/tffg_" + to_string(g) + "_modulo_" + to_string(modulo) + ".words";
  }
  else
  {
    basis_file_name = project_dir + "/tffg_" + to_string(g) + "_open_" + to_string(N + 1) + ".words";
  }

  TorsionFreeFuchsianGroup T = TorsionFreeFuchsianGroup(g);

  ifstream basis_file(basis_file_name);
  string line;

  G word = G(T.reduction);

  unordered_map<G,int,GHash> basis;

  int word_count = 0;
  int index=0;
  vector<int> word_rep;

  vector<G> generators;

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

  if (basis_file.is_open())
  {
    while (getline(basis_file, line))
    {
      word = G(T.reduction);
      word.identity();
      word.word = to_string(0) + " ";

      stringstream ss(line);

      ss >> word_count;

  	  while(ss >> index)
      {
        if(index>0)
        {
          word = word * generators[index-1];
          if (periodic_boundary)
            word = word % modulo;
        }
      }

      // basis.push_back(word);
	    basis.insert( {word, word_count} );
    }
  }
  else
  {
    throw runtime_error("Basis file not found! It needs to be generated first the generate_group command");
  }

  cout << "Basis file loaded succesfully! Basis dimension: " << basis.size() << endl;

  cout << "Generating the right regular representation." << endl;

  
  // // int i;
  G action = G(T.reduction);

  string output_file_name;
  if (periodic_boundary)
  {
    output_file_name = project_dir + "/tffg_" + to_string(g) + "_modulo_" + to_string(modulo) + "_";
  }
  else
  {
    output_file_name = project_dir + "/tffg_" + to_string(g) + "_open_" + to_string(N + 1) + "_";
  }

  vector<int> j_map(basis.size(), 0);

  int op_number = 0;
  for (G op : generators)
  {
    {
      // j = 0;
      // i = 0;
      vector<int> zero(basis.size(), 0);
      j_map = zero;

      unordered_map<G,int,GHash>::iterator basis_iterator = basis.begin();
      while(basis_iterator!=basis.end())
      {
        action = (basis_iterator->first)*op;

        if (periodic_boundary) action = action % modulo;

        unordered_map<G,int, GHash>::const_iterator it = basis.find(action);

	if (it != basis.end())
        {
          j_map[basis_iterator->second] = it->second;
        }
        else
        {
          if (periodic_boundary)
          {
            // -- must not happen for periodic boundary conditions
            throw runtime_error("Right regular representation failed for the word " + op.word);
          }

          // -- intercept boundary
          j_map[basis_iterator->second] = -1;
        }
        basis_iterator++;
      }
    }

    ofstream output_file;
    output_file.open(output_file_name + to_string(op_number++) + ".reg", ofstream::out | ofstream::trunc);

    // -- write result to file
    for (vector<int>::size_type j = 0; j < basis.size(); j++)
    {
      output_file << j_map[j] << " " << j << endl;
    }
    output_file.close();
  }

  string spec_info_fname = "";
  if (periodic_boundary)
  {
    spec_info_fname = project_dir + "/tffg_" + to_string(g) + "_modulo_" + to_string(modulo) + ".info";
  }
  else
  {
    spec_info_fname = project_dir + "/tffg_" + to_string(g) + "_open_" + to_string(N + 1) + ".info";
  }

  ofstream spec_info_file;
  spec_info_file.open(spec_info_fname, ofstream::out | ofstream::trunc);
  spec_info_file << basis.size() << endl;
  spec_info_file.close();

  return 0;
}
