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

extern vector<long> read_numberfield_reduction(const long &p, const long &q);

class Ring
{
public:
  vector<long> representation;
  vector<long> reduction;
  long dim;

  // Default constructor
  Ring(void);

  // Initializes zero
  Ring(const vector<long> &reduction);

  // Constructor with initialization
  Ring(const vector<long> &coeffs, const vector<long> &reduction);

  // Ring addition
  Ring operator+(const Ring &other) const;

  // Ring multiplication
  Ring operator*(const Ring &other) const;

  //Take modulus w.r.t to longeger m
  Ring operator%(const long &m) const;

  //Copy one Ring object to another
  void operator=(const Ring other);

  //Check if two ring objects are the same
  bool operator==(const Ring &other) const;

  // Pretty string representation
  string repr(void) const;
};

#endif
