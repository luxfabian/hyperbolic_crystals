/*
 *    ./include/fuchsian_ring.h
 *
 *    Author: Fabian R. Lux
 *    Date:   6/28/2023
 *
 *    Definition of the ring extension required for the Fuchsian group of genus two.
 */
#ifndef FUCHSIAN_RING_H
#define FUCHSIAN_RING_H

class Field
{

public:
    int a;
    int b;

    Field(void);

    Field(int a_in, int b_in);

    Field operator+(const Field &other) const;

    Field operator*(const Field &other) const;

    Field operator%(const int &m) const;

    void operator=(const Field other);

    bool operator==(const Field &other) const;

    void print(void);
};

#endif
