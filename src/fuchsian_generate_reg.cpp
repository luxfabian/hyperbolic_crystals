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
#include <unordered_map>
#include <algorithm>

#include <cmath>
#include <cstdlib>
#include <stdexcept>

#include "fuchsian_group.h"
#include "fuchsian_ring.h"

#pragma omp declare reduction(vec_int_plus : std::vector<int> : \
        std::transform(omp_out.begin(), omp_out.end(), omp_in.begin(), omp_out.begin(), std::plus<int>())) \
        initializer(omp_priv = decltype(omp_orig)(omp_orig.size()))

int main()
{
  using namespace std;

  // -- reading group secifications -----------------------------------

  const char *env = std::getenv("HYPERBOLIC_DIR");
  string project_dir = env;
 
  int p = 2;
  int N = 3;
  int m = pow(p, N);

  // -- read the input file -----------------------------------------

  string input_file_name("./fuchsian.inp");
  ifstream input_file(input_file_name);

  if (input_file.is_open())
  {
      // -- read
      input_file >> p;
      input_file >> N;
      m = pow(p, N);
  }
  else
  {
    throw runtime_error("The file fuchsian.inp was not found!");
  }
  input_file.close();

  cout << "#------------------------------------------------#" << endl;
  cout << "# C++ Spectra of Hyperbolic Cayley Crystals v0.1 #" << endl;
  cout << "#------------------------------------------------#" << endl;
  cout << "Fuchsian group of genus two modulo: " << m << endl;

  

  // Read the basis file
  cout << "Loading the basis file... " << endl;

  string basis_file_name(project_dir+"/fuchsian_"+to_string(m)+".words");
  ifstream basis_file(basis_file_name);

  string line;
  int gen = 0;
  Fuchsian F;

  vector<Fuchsian> basis;

  // -- initialize the group generators
	vector<Fuchsian> generators = get_generators(m);

  unordered_map<Fuchsian,int,FuchsianHash> unordered_basis;

  int word_count = 0;
  while (getline(basis_file, line))
  {
    istringstream input_stream(line);

    F = Fuchsian();

    while (input_stream)
    {
        // -- necessary to catch empty strings 
        gen = -1;

        // -- read letter
        input_stream >> gen;
        if(gen>0)
        {
            F = (F*generators[gen-1]) % m;
        }
    }
    basis.push_back(F);
    unordered_basis.insert( {F, word_count++} );
  }
  
  cout << "Done!" << endl;

  int d=basis.size();

  cout << "Dimension of the Hilbert space: " << d << endl;




  // TriangleGroup T = TriangleGroup(p, q);

  // ifstream basis_file(basis_file_name);
  // string line;

  // G word = G(T.reduction);

  // unordered_map<G,int,GHash> basis;

  // // Construct unit element
  // G A = T.A;
  // G B = T.B;
  // G AB = A * B;
  // G iA = G(T.reduction);
  
  // G iB = G(T.reduction);
  
  // iA.identity();
  // for (int k = 0; k < p - 1; k++)
  // {
  //   iA = iA * A;
  // }

  // iB.identity();
  // for (int k = 0; k < q - 1; k++)
  // {
  //   iB = iB * B;
  // }

  // if (periodic_boundary)
  // {
  //   A = A % modulo;
  //   B = B % modulo;
  //   AB = AB % modulo;
  //   iA = iA % modulo;
  //   iB = iB % modulo;
  // }

  // int word_count = 0;
  // string word_rep;
  // if (basis_file.is_open())
  // {
  //   while (getline(basis_file, line))
  //   {
  //     word = G(T.reduction);
  //     word.identity();
  //     word_rep = "";

  //     stringstream ss(line);

  //     ss >> word_count;
  //     ss >> word_rep;

  //     for (char &c : word_rep)
  //     {
  //       if (c == 'A')
  //       {
  //         word = word * A;
  //       }
  //       else if (c == 'B')
  //       {
  //         word = word * B;
  //       }
  //       else if (c == 'E')
  //       {
  //       }
  //       if (periodic_boundary)
  //         word = word % modulo;
  //     }
  //     // basis.push_back(word);
	//   basis.insert( {word, word_count} );
  //   }
  // }
  // else
  // {
  //   throw runtime_error("Basis file not found! It needs to be generated first the generate_group command");
  // }

  // cout << "Basis file loaded succesfully! Basis dimension: " << basis.size() << endl;

  // cout << "Generating the right regular representation." << endl;

  // vector<G> operators;

  // // Account for the inverse operation which appears in the right-regular representation
  // A.word = "iA";
  // B.word = "iB";
  // iA.word = "A";
  // iB.word = "B";

  // operators.push_back(A);
  // operators.push_back(B);
  // operators.push_back(AB);
  // operators.push_back(iA);
  // operators.push_back(iB);

  // // int i;
  // G action = G(T.reduction);

  // string output_file_name;
  // if (periodic_boundary)
  // {
  //   output_file_name = project_dir + "/" + to_string(p) + "_" + to_string(q) + "_modulo_" + to_string(modulo) + "_";
  // }
  // else
  // {
  //   output_file_name = project_dir + "/" + to_string(p) + "_" + to_string(q) + "_open_" + to_string(N + 1) + "_";
  // }

  vector<int> j_map(basis.size(), 0);

  Fuchsian action;

  int generator_index = 1;

  for (Fuchsian op : generators)
  {
    {
      // j = 0;
      // i = 0;
      vector<int> zero(basis.size(), 0);
      j_map = zero;

      unordered_map<Fuchsian,int,FuchsianHash>::iterator basis_iterator = unordered_basis.begin();
      while(basis_iterator!=unordered_basis.end())
      {
        action = (basis_iterator->first)*op % m;

        unordered_map<Fuchsian,int, FuchsianHash>::const_iterator it = unordered_basis.find(action);

	      if (it != unordered_basis.end())
        {
          j_map[basis_iterator->second] = it->second;
        }
        else
        {
          // -- must not happen for periodic boundary conditions
          throw runtime_error("Right regular representation failed. Group element mapped to nowhere");
        }
        basis_iterator++;
      }
    }
  
    ofstream output_file;
    output_file.open(project_dir+"/fuchsian_"+to_string(m) + "_" + to_string(generator_index++) + ".reg", ofstream::out | ofstream::trunc);

    // -- write result to file
    for (vector<int>::size_type j = 0; j < basis.size(); j++)
    {
      output_file << j_map[j] << " " << j << endl;
    }
    output_file.close();
  }

  string spec_info_fname = "";

  spec_info_fname = project_dir+"/fuchsian_"+to_string(m) + ".info";
  

  ofstream spec_info_file;
  spec_info_file.open(spec_info_fname, ofstream::out | ofstream::trunc);
  spec_info_file << basis.size() << endl;
  spec_info_file.close();

  return 0;
}
