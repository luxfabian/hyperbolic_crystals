/*
 *    ./include/triangle.group.h
 *
 *    Author: Fabian R. Lux
 *    Date:   04/02/2022
 *
 *    Implements a representation of the proper triangle group Delta+(p,q,2)
 *    which can be used to describe each group element g as a matrix in the 
 *    space SL(3,Z[2 cos(pi/pq)]). The group is given by the presentation
 *
 *    Delta+(p,q,2) = < A,B | A**p = B**q = (AB)**2 >. 
 *
 *    The representation is only given for the generators. This means that
 *    other group elements still need to be generated by an external routine.
 *
 *    The classes presented here only enable the basic algebraic manipulations
 *    which are necessary for this goal. 
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
    string word = "";
    vector<int> reduction;

    // -- matrix representation of the word
    Ring mat[3][3];

    // Default constructor
    G(void);

    // Initialize with known field specs
    G(const vector<int> reduction);

    // Modulo operator
    G operator%(const int &m) const;

    // Multiplication operator
    G operator*(const G &other) const;

    // Assignment operator
    void operator=(const G &other);

    // Equality operator
    bool operator==(const G &other) const;

    // Pretty string representation
    string repr(void) const;

    // Set G to be equal to the identity element
    void identity(void);
};

extern size_t hash_combine(size_t lhs, size_t rhs);

class GHash
{
public:
  size_t operator()(const G &obj) const;
};


// concatenate n copies of a word
string word_power(const string &word, const int &n);

class TriangleGroup
{
public:
    // Signature of the triangle group
    int p, q;
    const int r = 2;

    // Polynomimal reduction for ring extension
    vector<int> reduction;

    // Generators of the triangle group
    G X;
    G Y;
    G Z;

    // Generators of the proper triangle group
    G A;
    G B;

    // Encoding of the group relations
    vector<string> relations;

    // Constructor
    TriangleGroup(const int &p, const int &q);

    // Uses relations of the group to reduce the given word
    void reduce(string &word);
private:
    // Reads A and B for given p and q. Is called by constructor
    void read_generators(void);
};

#endif
