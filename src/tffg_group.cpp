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

TorsionFreeFuchsianGroup(const long &g)
{
  // initialize parent group
  this->Delta = TriangleGroup(4*g,4*g);

  // first translation operator
  G gamma;

  // rotations
  G A = Delta.A;
  G A_inverse =  Delta.A;

  // inverse can be constructed algebraically
  for(int i=0; i<4*g-1; i++)
  {
    A_inverse = A_inverse * Delta.A;
  }

  // initializing gamma
  gamma = Delta.C;

  for(int i=0; i<2*g; i++)
  {
    gamma = gamma * Delta.A;
  }

  // add first element to list of generators
  (this -> gamma).push_back(gamma);

  for(int i=0; i<4*g-1; i++)
  {
    gamma = A_inverse * (this -> gamma).back() * A;
    (this -> gamma).push_back(gamma);
  }
}