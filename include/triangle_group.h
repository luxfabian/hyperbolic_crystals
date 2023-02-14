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
#include "ring_extension.h"

using namespace std;

// Represents an element of the triangle group
class G
{
public:
    string word;
    vector<int> reduction;

    // -- matrix representation of the word
    Ring mat[3][3];

    G(const vector<int> reduction);

    // Modulo operator
    G operator%(const int &m);

    // Multiplication operator
    G operator*(const G &other);

    // Assignment operator
    void operator=(const G &other);

    // Equality operator
    bool operator==(const G &other) const;
};

// extern size_t hash_combine(size_t lhs, size_t rhs);

// class GHash
// {
//     public:

//     size_t operator()(constG &fuchs) const;
// };


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
