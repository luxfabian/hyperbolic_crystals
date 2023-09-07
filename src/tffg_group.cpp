/*
 *    ./src/tffg_group.cpp
 *
 *    Author:  Fabian R. Lux
 *    Date:   08/01/2023
 *
 */
#include <stdexcept>
#include "tffg_group.h"
#include "triangle_group.h"

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>

#include <vector>

#include <cstdlib>
 

using namespace std;

TorsionFreeFuchsianGroup::TorsionFreeFuchsianGroup(const long &g)
{
  // initialize parent group
  TriangleGroup Delta = TriangleGroup(4*g, 4*g);
  this->reduction = Delta.reduction;
  
  // genus
  this->g = g;
  
  // first translation operator
  G gamma = G(this->reduction);
  
  // rotations
  G A = Delta.A;
  G A_inverse =  power(Delta.A, 4*g-1);

  // initializing first gamma
  gamma = Delta.C * power(Delta.A, 2*g);
  gamma.word = to_string(1) + " ";

  // add first element to list of generators
  (this -> gamma).push_back(gamma);

  G buffer;
  for(int n=1; n<4*g; n++)
  {
    buffer = ( power(A_inverse, (2*g-1) * n) * gamma ) * power(A, (2*g-1) * n);
    buffer.word = to_string(n+1) + " ";
    (this -> gamma).push_back(buffer);
  }
}