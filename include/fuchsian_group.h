/*
 *    ./include/fuchsian_group.h
 *
 *    Author: Fabian R. Lux
 *    Date:   6/28/2023
 *
 *    Definition of the Fuchsian group of genus two.
 */
#ifndef FUCHSIAN_H
#define FUCHSIAN_H

#include <vector>
#include "fuchsian_ring.h"

using namespace std;

extern Field generator[8][2][2];

class Fuchsian
{
public:
    // An unreduced word in the Fuchsian lexikon
    vector<int> word = {0};

    // Matrix representation of the word
    Field mat[4];

    // Constructor
    Fuchsian(void);

    // Right group action w.r.t. generator of given index
    void right_action(int generator_id);

    Fuchsian operator%(const int &m) const;

    Fuchsian operator*(const Fuchsian &other) const;

    void operator=(const Fuchsian &other);

    bool operator==(const Fuchsian &other) const;

    // Update mat
    void assign(Field F[4]);
};

extern size_t hash_combine(size_t lhs, size_t rhs);

class FuchsianHash
{
    public:

    size_t operator()(const Fuchsian &fuchs) const;
};

extern vector<Fuchsian> get_generators(const int &m);

#endif
