/*
 *    ./include/triangle.group.h
 *
 *    Author: Fabian R. Lux
 *    Date:   04/02/2022
 *
 *    Implement the relation groups of the triangle group
 *    with signature {p,q,2}, which corresponds to tesselation
 *    of te hyperbolic plane.
 */
#ifndef TRIANGLE_GROUP_H
#define TRIANGLE_GROUP_H

#include <vector>
#include <string>

using namespace std;


// concatenate n copies of a word
string word_power(const string &word, const int &n);

class TriangleGroup
{
public:
    // -- signature
    int p, q;
    const int r = 2;

    vector<string> relations;

    // -- constructor
    TriangleGroup(const int &p, const int &q);

    // -- apply relations of the group
    void reduce(string &word);
};


#endif
