#include <iostream>
#include <fstream>
#include <string>
#include <unordered_set>
#include "triangle_group.h"

void add_to_basis(string &candidate, TriangleGroup &T,  unordered_set<string> &unordered_basis,  vector<string> &basis,  vector<string> &next_generation){
  T.reduce(candidate);
  if (unordered_basis.find(candidate) == unordered_basis.end())
  {
    unordered_basis.insert(candidate);
    next_generation.push_back(candidate);
    basis.push_back(candidate);
  }
}

int main(){

  using namespace std;

  int p=5;
  int q=4;
  int max_generation=4;

  // -- read the input file -------------------------------------------

  string input_file_name("./input.inp");
  ifstream input_file(input_file_name);

  if (input_file.is_open())
  {
      // -- read
      input_file >> p;
      input_file >> q;
      input_file >> max_generation;
  }

  // -- set up the group ----------------------------------------------
  TriangleGroup T = TriangleGroup(p,q);

  // -- found elements are stored here
  unordered_set<string> unordered_basis;
  vector<string> basis;

  // -- keep track of the iterative process
	vector<string> prev_generation;
	vector<string> next_generation;

  // -- initialize identity element
	// unordered_basis.insert("");
  // basis.push_back("");
  prev_generation.push_back("");

  string generator;
  string candidate;
  int order;

  cout << "Generation | # of found words" << endl;

  for(int g=0; g<max_generation; g++)
  {
    cout << g << "\t | \t" << basis.size() << endl;

    if(g%2==0)
    {
      generator = "A";
      order = p;
    }
    else
    {
      generator = "B";
      order = q;
    }
   
    for(string word: prev_generation)
    {
      // -- complete cycles
      candidate = word;
      for(int k=0; k<order; k++)
      {
        candidate = candidate + generator;

        add_to_basis(candidate, T, unordered_basis,  basis, next_generation);
        
      }
    }

    prev_generation = next_generation;
  }
  
  string output_file_name("./{"+to_string(p)+","+to_string(q)+"}_"+to_string(max_generation)+".words");
	ofstream output_file;
	output_file.open(output_file_name, ofstream::out | ofstream::trunc);
	
	for(string word: basis)
	{
		output_file << word << endl;
	}

	output_file.close();

  return 0;
}
