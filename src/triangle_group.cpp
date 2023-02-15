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

  vector<int> unity(reduction.size(), 0);
  unity[0] = 1;

  for (int i = 0; i < 3; i++)
  {
    for (int j = 0; j < 3; j++)
    {
      if (i == j)
      {
        this->mat[i][i] = unity;
      }
      else
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
      this->B.mat[i][j] = Ring(B_flat[k], this->reduction);
      k++;
    }
  }

  this->A.word +"A";
  this->B.word +"B";
}
