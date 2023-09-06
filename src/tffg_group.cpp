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
  TriangleGroup Delta = TriangleGroup(4*g, 4*g);
  // initialize parent group

  this->g = g;

  // first translation operator
  G gamma;
  G buffer;

  // rotations
  G A = Delta.A;
  G A_inverse =  power(Delta.A, 4*g-1);

  // initializing first gamma
  gamma = Delta.C * power(Delta.A, 2*g);

  // add first element to list of generators
  (this -> gamma).push_back(gamma);

  for(int n=1; n<4*g; n++)
  {
    buffer = ( power(A_inverse, (2*g-1) * n) * gamma ) * power(A, (2*g-1) * n);
    (this -> gamma).push_back(buffer);
  }
}