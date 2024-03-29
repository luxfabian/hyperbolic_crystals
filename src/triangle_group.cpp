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
// #include <boost/algorithm/string/replace.hpp> // include Boost, a C++ library

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

G::G(const vector<long> reduction)
{
  this->reduction = reduction;
  this->word = "";

  for (long i = 0; i < 3; i++)
  {
    for (long j = 0; j < 3; j++)
    {
      {
        this->mat[i][j] = Ring(reduction);
      }
    }
  }
}

void G::identity(void)
{
  this->word="";

  vector<long> c(this->reduction.size(), 0);

  Ring zero = Ring(c, this->reduction);
  
  c[0]=1;

  Ring id = Ring(c,this->reduction);
  

  for (long i = 0; i < 3; i++)
  {
    for (long j = 0; j < 3; j++)
    {
      if(i==j)
      {
        this->mat[i][i] = id;
      }
      else
      {
        this->mat[i][j] = zero;
      }
    }
  }
}

G G::operator%(const long &m) const
{
  G result = G(this->reduction);

  result.word = this->word;

  for (long i = 0; i < 3; i++)
  {
    for (long j = 0; j < 3; j++)
    {
      result.mat[i][j] = this->mat[i][j] % m;
    }
  }

  return result;
}

void G::operator=(const G &other)
{
  this->word = other.word;
  this->reduction = other.reduction;

  for (long i = 0; i < 3; i++)
  {
    for (long j = 0; j < 3; j++)
    {
      this->mat[i][j] = other.mat[i][j];
    }
  }
}

G G::operator*(const G &other) const
{
  G result=G(this->reduction);

  result.word = this->word + other.word;
  // result.reduction = other.reduction;

  for (long i = 0; i < 3; i++)
  {
    for (long j = 0; j < 3; j++)
    {
      for (long k = 0; k < 3; k++)
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

  for (long i = 0; i < 3; i++)
  {
    for (long j = 0; j < 3; j++)
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

string G::repr(void) const
{
  string rep;

  for (long i=0; i<3; i++)
  {
    for (long j=0; j<3; j++)
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

    std::hash<long> hasher;
    size_t seed = 0;

    vector<long> vals;

    for (long i = 0; i < 3; i++)
    {
      for(long j = 0; j < 3; j++)
      {
	for(long c: obj.mat[i][j].representation)
	{
          seed = hash_combine(seed, hasher(c));
	}
      }
    }

    return seed;
}


G power(const G &groupelement, const long &n)
{
  G result = G(groupelement.reduction);
  result.identity();
  
  for(long i=0; i<n; i++)
  {
    result = result * groupelement;
  }

  return result;
}

G power(const G &groupelement, const long &n, const long &mod)
{
  G result = G(groupelement.reduction);
  result.identity();
  
  for(long i=0; i<n; i++)
  {
    result = (result * groupelement) % mod;
  }

  return result;
}


// concatenate n copies of a word
string word_power(const string &word, const long &n)
{
  string result("");

  for (long i = 0; i < n; i++)
  {
    result += word;
  }

  return result;
}

TriangleGroup::TriangleGroup(const long &p, const long &q)
{
  // Signature
  this->p = p;
  this->q = q;  

  this->reduction = read_numberfield_reduction(p,q);

  // bool triangle_is_admissble = (1.0 / p + 1.0 / q < 0.5);

  // if (not triangle_is_admissble)
  // {
  //   throw std::invalid_argument("Triangle is not admissible since 1/p +1/q < 1/2 is not fulfilled.");
  // }

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
    // boost::algorithm::replace_all(word, relation, "");
  }
}

void TriangleGroup::read_generators(void)
{

  const char* env = std::getenv("HYPERBOLIC_BUILD");
  string build_dir = env;

  string file_name_X = build_dir +"/generators/"+to_string(this->p)+"_"+to_string(this->q)+".X";
  string file_name_Y = build_dir +"/generators/"+to_string(this->p)+"_"+to_string(this->q)+".Y";
  string file_name_Z = build_dir +"/generators/"+to_string(this->p)+"_"+to_string(this->q)+".Z";

  ifstream file_X(file_name_X);
  ifstream file_Y(file_name_Y);
  ifstream file_Z(file_name_Z);

  string line;
  long buffer;
  vector<long> c;

  this->X = G(this->reduction);
  this->Y = G(this->reduction);
  this->Z = G(this->reduction);

  // 1D representation of the matrix; needs to be re-shaped at the end
  vector<vector<long>> X_flat;
  vector<vector<long>> Y_flat;
  vector<vector<long>> Z_flat;

  // -- Read in X -----------------------------------------------------
  if(file_X.is_open())
  {
    while(getline(file_X, line))
    {
      istringstream iss(line);
      c.clear();
      while(iss>>buffer)
      {
        c.push_back(buffer);
      }
      X_flat.push_back(c);
    }
  }
  else
  {
    throw runtime_error(file_name_X+" not found!");
  }
  file_X.close();

  // -- Read in Y -----------------------------------------------------
  if(file_Y.is_open())
  {
    while(getline(file_Y, line))
    {
      istringstream iss(line);
      c.clear();
      while(iss>>buffer)
      {
        c.push_back(buffer);
      }
      Y_flat.push_back(c);
    }
  }
  else
  {
    throw runtime_error(file_name_Y+" not found!");
  }
  file_Y.close();

  // -- Read in Z -----------------------------------------------------
  if(file_Z.is_open())
  {
    while(getline(file_Z, line))
    {
      istringstream iss(line);
      c.clear();
      while(iss>>buffer)
      {
        c.push_back(buffer);
      }
      Z_flat.push_back(c);
    }
  }
  else
  {
    throw runtime_error(file_name_Z+" not found!");
  }
  file_Z.close();

  long k=0;
  for(long i=0; i<3; i++)
  {
    for(long j=0; j<3; j++)
    {
      this->X.mat[i][j] = Ring(X_flat[k], this->reduction);
      this->Y.mat[i][j] = Ring(Y_flat[k], this->reduction);
      this->Z.mat[i][j] = Ring(Z_flat[k], this->reduction);
      k++;
    }
  }

  this->X.word = "X";
  this->Y.word = "Y";
  this->Z.word = "Z";

  // Construct generator of proper triangle group
  this->A = X*Y;
  this->B = Y*Z;
  this->C = X*Z;

  // Overwrite the word
  this->A.word = "A";
  this->B.word = "B";
  this->C.word = "C";
}
