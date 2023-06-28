/*
 * ./src/fuchsian_generate_group.cpp
 *
 * Author:  Fabian R. Lux
 * Date:    05/01/2023
 *
 * Find a finite approximation for the Fuchsian group of genus 2.
 * The approximation is obtained by first generation a faithful
 * representation of the group with the image in GL(2,F), where
 * F is the ring obtained by adjoining the square root of three
 * to the ring of integers.
 * 
 * This representation is subsequently modded w.r.t to m = p^N,
 * where p and N are integers. The default is p=2 and N=3. It
 * can be overwritten by placing the input file 'input.inp'
 * in the current directory with p in the first line and N in the
 * second.
 */
#include <iostream>
#include <vector>
#include <string>
#include <cmath>
#include <fstream>
#include <unordered_set>
#include "fuchsian_group.h"	
#include "fuchsian_ring.h"

int main(int argc, char **argv)
{
	using namespace std;

	int p = 2;
	int N = 3;
	int m = pow(p, N);

	// -- read the input file -----------------------------------------

	const char *env = std::getenv("HYPERBOLIC_DIR");
  	string project_dir = env;

	string input_file_name(project_dir+"/fuchsian.inp");
	ifstream input_file(input_file_name);

	if (input_file.is_open())
	{
		// -- read
		input_file >> p;
		input_file >> N;
		m = pow(p, N);
	}

	cout << "Fuchsian group of genus 2" << endl;

	cout << "Finite approximation w.r.t. mod " << m << endl;

	cout << "Initializing the basis set..." << endl;

	// -- found elements are stored here
	std::unordered_set<Fuchsian, FuchsianHash> basis;

	// -- keep track of the iterative process
	vector<Fuchsian> prev_generation;
	vector<Fuchsian> next_generation;

	// -- initialize identity element
	Fuchsian candidate;
	basis.insert(candidate);
	prev_generation.push_back(Fuchsian());

	// -- initialize the group generators
	vector<Fuchsian> generators = get_generators(m);

	// -- iterative group generation
	int g = 0;
	while (true)
	{
		cout << g << "\t |" << basis.size() << endl;

		next_generation.clear();
		for (auto it = prev_generation.begin(); it != prev_generation.end(); ++it)
		{
			for (int g = 0; g < 8; g++)
			{
				candidate = (*it) * generators[g];
				candidate = candidate % m;
				if (basis.find(candidate) == basis.end())
				{
					basis.insert(candidate);
					next_generation.push_back(candidate);
				}
			}
		}

		if (next_generation.size() == 0)
		{
			break;
		}

		prev_generation = next_generation;
		g += 1;
	}


	string output_file_name(project_dir+"/fuchsian_"+to_string(m)+".words");
	ofstream output_file;
	output_file.open(output_file_name, ofstream::out | ofstream::trunc);
	
	for(const auto& elem: basis)
	{
		for(int w: elem.word)
		{
			output_file << w << " ";
		}
		output_file << endl;
	}

	output_file.close();
	
	return 0;
}
