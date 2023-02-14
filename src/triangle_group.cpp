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
#include "triangle_group.h"
#include <boost/algorithm/string.hpp> // include Boost, a C++ library

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
  // -- set signature
  this->p = p;
  this->q = q;

  bool triangle_is_admissble = (1.0 / p + 1.0 / q < 0.5);

  if (not triangle_is_admissble)
  {
    throw std::invalid_argument("Triangle is not admissible since 1/p +1/q < 1/2 is not fulfilled.");
  }

  // generate relations
  this->relations.push_back(word_power("A", this->p));
  this->relations.push_back(word_power("B", this->q));
  this->relations.push_back(word_power("AB", 2));
}

// apply relations of the group
void TriangleGroup::reduce(string &word)
{
  for(string relation: this->relations)
  {
    boost::replace_all(word, relation, "");
  }
}
