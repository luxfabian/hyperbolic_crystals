/*
 *    ./include/ring_extension.h
 *
 *    Author: Fabian R. Lux
 *    Date:   12/02/2022
 *
 *    Constructs the ring extension Z[x] where x^d fullfills the relation
 *    x^d = r(0) x^0 + r(1) x^1 + ... r(d-1) x^(d-1).
 */
#ifndef RING_EXTENSION_H
#define RING_EXTENSION_H

#include <vector>

using namespace std;

class Ring
{
public:
  vector<int> representation;
  vector<int> reduction;
  int dim;

  //-- constructor
  Ring(void);

  //-- constructor with initialization
  Ring(const vector<int> &coeffs);

  //-- ring operations
  Ring operator+(const Ring &other);
  Ring operator*(const Ring &other);

  //-- modulo operator
  Ring operator%(const int &m);

  //-- assignment
  void operator=(const Field other);

  //-- equality
  bool operator==(const Field &other) const;
};

#endif
