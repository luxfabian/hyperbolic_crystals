/*
 *    ./src/ring_extension.cpp
 *
 *    Author:  Fabian R. Lux
 *    Date:    12/01/2023
 *
 *    Constructs the ring extension Z[x] where x^d fullfills the relation
 *    x^d = r(0) x^0 + r(1) x^1 + ... r(d-1) x^(d-1).
 */

#include "ring_extension.h"
#include <stdexcept>

using namespace std;

// Initialize as zero
Ring::Ring(void)
{
  vector<int> empty(this->dim,0);
  this->representation = empty;
}

// Initialize with a given representation
Ring::Ring(const vector<int> &coeffs)
{
  if(coeffs.size()==this->dim)
  {
    this->representation = coeffs;
  }
  else
  {
    throw invalid_argument("Representation vector has wrong size");
  }
}

// Ring addition
Ring Ring::operator+(const Ring &other)
{
  return Ring(this->representation + other.representation);
}

// Ring multiplication
Ring Ring::operator*(const Ring &other)
{
  const vector<int> &a = this->representation;
  const vector<int> &b = other.representation;

  const int &d = this->dim;

  vector<int> buffer(2*d-1,0);

  // -- unreduced product of polynomials
  for(int i=0; i<d; i++)
  {
    for(int j=0; j<d; j++)
    {
      buffer[i+j] += a[i]*b[j];
    }
  }

  // -- folding back
  int b=0;
  for(int k=0; k<d-1; k++)
  {
    // -- store buffer element temporarily
    b=buffer[2*d-2-k]
    buffer[2*d-2-k]=0;

    // -- fold this deleted index back onto the buffer
    for(int l=0; l<d-1; l++)
    {
      // -- index runs from id - (d-1) until id-1
      buffer[d-2-k+l] += b * this->reduction[l]
    }
  }

  vector<int> c;
  for(int l=0; l<d-1; l++)
  {
    c.push_back(buffer[l])  
  }

  return Ring(c);
}

//Take modulus of the ring w.r.t to m
Ring Ring::operator%(const int &m)
{
  //https://stackoverflow.com/a/44197900/7236657

  vector <int> mod;

  for(int i=0; i<(this->dim); i++)
  {
    mod.push_back( m + ((this->representation)[i] % m) % m);
  }

  return mod;
}

//Copy one field value to another
void Ring::operator=(const Ring other)
{
  this->representation = other.representation;
}

//Check if two field values are the same
bool Ring::operator==(const Ring &other) const
{
  return this->representation == other.representation
}
