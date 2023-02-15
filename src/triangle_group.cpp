/*
 *    ./src/triangle_group.cpp
 *
 *    Author:  Fabian R. Lux
 *    Date:   04/02/2022
 *
 *    Implement the relation groups of the triangle group
 *    with signature {p,q,2}, which corresponds to tesselation
 *    of te hyperbolic plane.
 */
#include <stdexcept>
#include "ring_extension.h"
#include "triangle_group.h"
#include <boost/algorithm/string.hpp> // include Boost, a C++ library

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>

#include <vector>

#include <cstdlib>
 

using namespace std;

G::G(void)
{
}

G::G(const vector<int> reduction)
{
  this->reduction = reduction;
  this->word = "";

  // vector<int> unity(reduction.size(), 0);
  // unity[0] = 1;

  for (int i = 0; i < 3; i++)
  {
    for (int j = 0; j < 3; j++)
    {
      {
        this->mat[i][j] = Ring(reduction);
      }
    }
  }
}

G G::operator%(const int &m)
{
  G result = G(this->reduction);

  result.word = this->word;

  for (int i = 0; i < 3; i++)
  {
    for (int j = 0; j < 3; j++)
    {
      result.mat[i][j] = this->mat[i][j] % m;
    }
  }

  return result;
}

void G::operator=(const G &other)
{
  this->word = other.word;

  for (int i = 0; i < 3; i++)
  {
    for (int j = 0; j < 3; j++)
    {
      this->mat[i][j] = other.mat[i][j];
    }
  }
}

G G::operator*(const G &other)
{
  G result=G(this->reduction);

  result.word = this->word + other.word;

  for (int i = 0; i < 3; i++)
  {
    for (int j = 0; j < 3; j++)
    {
      for (int k = 0; k < 3; k++)
      {
	result.mat[i][j] = result.mat[i][j] + this->mat[i][k] * other.mat[k][j];
      }
    }
  }

  return result;
}

bool G::operator==(const G &other) const
{

  bool equal = true;

  for (int i = 0; i < 3; i++)
  {
    for (int j = 0; j < 3; j++)
    {
      equal &= (this->mat[i][j] == other.mat[i][j]);

      if (not equal)
      {
        break;
      }
    }
  }

  return equal;
}

string G::repr(void)
{
  string rep;

  for (int i=0; i<3; i++)
  {
    for (int j=0; j<3; j++)
    {
      rep += this->mat[i][j].repr() + "  ";
    }
    rep += "\n";
  }

  return rep;
}

/*
 * Generate a new hash from two known hashes lhs and rhs.
 *
 * The constants which are added represent
 * 
 *      1. 0x517cc1b727220a95 = inverse of pi 
 *      2. 0x9e3779b9 = inverse golden ratio
 * 
 * Their bit representations consist of a random sequence of
 * 0's and 1's which is ideal for scrambling up the hashes.
 * 
 * Together, the bit-shifts and the addition of the constant above
 * break the symmetry of the XOR operator ^ and lead to a reasonable
 * way of combining to hashes.
*/
size_t hash_combine(size_t lhs, size_t rhs)
{
    // https://stackoverflow.com/a/27952689/7236657
    if constexpr (sizeof(size_t) >= 8)
    {
        // 64-bit architecture
        lhs ^= rhs + 0x517cc1b727220a95 + (lhs << 6) + (lhs >> 2);
    }
    else
    {
        // 32-bit architecture
        lhs ^= rhs + 0x9e3779b9 + (lhs << 6) + (lhs >> 2);
    }
    return lhs;
}

size_t GHash::operator()(const G &obj) const
{
    // https://stackoverflow.com/a/2595226/7236657

    std::hash<int> hasher;
    size_t seed = 0;

    vector<int> vals;

    for (int i = 0; i < 3; i++)
    {
      for(int j = 0; j < 3; j++)
      {
	for(int c: obj.mat[i][j].representation)
	{
          seed = hash_combine(seed, hasher(c));
	}
      }
    }

    return seed;
}



// concatenate n copies of a word
string word_power(const string &word, const int &n)
{
  string result("");

  for (int i = 0; i < n; i++)
  {
    result += word;
  }

  return result;
}

TriangleGroup::TriangleGroup(const int &p, const int &q)
{
  // Signature
  this->p = p;
  this->q = q;  

  this->reduction = read_numberfield_reduction(p,q);

  bool triangle_is_admissble = (1.0 / p + 1.0 / q < 0.5);

  if (not triangle_is_admissble)
  {
    throw std::invalid_argument("Triangle is not admissible since 1/p +1/q < 1/2 is not fulfilled.");
  }

  // -- generate relations
  this->relations.push_back(word_power("A", this->p));
  this->relations.push_back(word_power("B", this->q));
  this->relations.push_back(word_power("AB", 2));

  this->read_generators();
}

// apply relations of the group
void TriangleGroup::reduce(string &word)
{
  for (string relation : this->relations)
  {
    boost::replace_all(word, relation, "");
  }
}

void TriangleGroup::read_generators(void)
{


  const char* env = std::getenv("HYPERBOLIC_BUILD");
  string build_dir = env;

  string file_name_A = build_dir +"/generators/"+to_string(this->p)+"_"+to_string(this->q)+".A";
  string file_name_B = build_dir +"/generators/"+to_string(this->p)+"_"+to_string(this->q)+".B";
  ifstream file_A(file_name_A);
  ifstream file_B(file_name_B);

  string line;
  int buffer;
  vector<int> c;

  this->A = G(this->reduction);
  this->B = G(this->reduction);

  // 1D representation of the matrix; needs to be re-shaped at the end
  vector<vector<int>> A_flat;
  vector<vector<int>> B_flat;

  // Read in A
  if(file_A.is_open())
  {
    while(getline(file_A, line))
    {
      istringstream iss(line);
      c.clear();
      while(iss>>buffer)
      {
        c.push_back(buffer);
      }
      A_flat.push_back(c);
    }
  }
  else
  {
    throw runtime_error(file_name_A+" not found!");
  }
  file_A.close();

  // Read in B
  if(file_B.is_open())
  {
    while(getline(file_B, line))
    {
      istringstream iss(line);
      c.clear();
      while(iss>>buffer)
      {
        c.push_back(buffer);
      }
      B_flat.push_back(c);
    }
  }
  else
  {
    throw runtime_error(file_name_B+" not found!");
  }
  file_B.close();

  int k=0;
  for(int i=0; i<3; i++)
  {
    for(int j=0; j<3; j++)
    {
      this->A.mat[i][j] = Ring(A_flat[k], this->reduction);
      // cout << (this->A.mat[i][j].repr()) << endl;
      this->B.mat[i][j] = Ring(B_flat[k], this->reduction);
      k++;
    }
  }

  this->A.word = "A";
  this->B.word = "B";
}
