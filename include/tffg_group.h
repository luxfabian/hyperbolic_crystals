/*
 *    ./include/tffg_group.h
 *
 *    Author: Fabian R. Lux
 *    Date:   8/1/2023
 *
 */
#ifndef TFFG_GROUP_H
#define TFFG_GROUP_H

#include <vector>
#include <string>
#include "triangle_group.h"

using namespace std;


class TorsionFreeFuchsianGroup
{
public:
    // Signature of the triangle group
    long g;

    // Polynomimal reduction for ring extension
    vector<long> reduction;

    // Generators of the group
    vector<G> gamma;
    
    // Constructor
    TorsionFreeFuchsianGroup(const long &g);
};

#endif
