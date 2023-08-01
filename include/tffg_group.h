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
    long p, q;

    //Parent group
    TriangleGroup Delta;

    // Polynomimal reduction for ring extension
    vector<long> reduction;

    // Generators of the group
    vector<G> gamma;
    
    // Constructor
    TorsionFreeFuchsianGroup(const long &p, const long &q);

private:
    // Reads A and B for given p and q. Is called by constructor
    void read_generators(void);
};

#endif
