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

#include <string>
#include <vector>

using namespace std;

extern vector<int> read_numberfield_reduction(const int &p, const int &q);

class Ring
{
public:
  vector<int> representation;
  vector<int> reduction;
  int dim;

  // Default constructor
  Ring(void);

  // Initializes zero
  Ring(const vector<int> &reduction);

  // Constructor with initialization
  Ring(const vector<int> &coeffs, const vector<int> &reduction);

  // Ring operations
  Ring operator+(const Ring &other);
  Ring operator*(const Ring &other);
  Ring operator%(const int &m);
  void operator=(const Ring other);
  bool operator==(const Ring &other) const;

  // Pretty string representation
  string repr(void);
};

#endif
