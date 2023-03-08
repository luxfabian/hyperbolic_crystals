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

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

#include <cstdlib>

using namespace std;

vector<int> read_numberfield_reduction(const int &p, const int &q)
{
  string line;

  const char *env = std::getenv("HYPERBOLIC_BUILD");
  string build_dir = env;

  ifstream file(build_dir + "/ring_reduction.inp");

  int n = 2 * p * q;

  // cout << "Searching ring_reduction.inp for n=2pq=" << n << endl;

  vector<int> reduction;
  int buffer;
  int d = 0;

  bool listening = false;

  if (file.is_open())
  {
    while (getline(file, line))
    {
      //--start analysing after keyword is found
      if (line.find("BEGIN") != string::npos)
      {
        listening = true;
      }

      //--stop analysing after keyword is found
      if (line.find("END") != string::npos)
      {
        break;
      }

      if (listening)
      {
        istringstream iss(line);

        if (iss >> buffer && iss >> d && buffer == n)
        {
          while (iss >> buffer)
          {
            reduction.push_back(buffer);
          }
          //--no need to keep listening
          break;
        }
      }
    }
  }
  else
  {
    throw runtime_error("The file ring_reduction.inp was not found!");
  }
  file.close();

  // cout << "Dimension of ring extension is d=" << d << endl;

  return reduction;
}

Ring::Ring(void)
{
  this->dim = 0;
}

Ring::Ring(const vector<int> &reduction)
{
  this->dim = reduction.size();
  vector<int> empty(this->dim, 0);
  this->representation = empty;
  this->reduction = reduction;
}

Ring::Ring(const vector<int> &coeffs, const vector<int> &reduction)
{
  if (coeffs.size() == reduction.size())
  {
    this->dim = reduction.size();
    this->representation = coeffs;
    this->reduction = reduction;
  }
  else
  {
    throw invalid_argument("Representation vector has wrong size");
  }
}

Ring Ring::operator+(const Ring &other)
{
  vector<int> add;

  for (vector<int>::size_type i = 0; i < (this->representation.size()); i++)
  {
    add.push_back(this->representation[i] + other.representation[i]);
  }

  return Ring(add, this->reduction);
}

Ring Ring::operator*(const Ring &other) const
{
  const vector<int> &a = this->representation;
  const vector<int> &b = other.representation;

  const int &d = this->dim;

  vector<int> buffer(2 * d - 1, 0);

  // -- unreduced product of polynomials
  for (int i = 0; i < d; i++)
  {
    for (int j = 0; j < d; j++)
    {
      buffer[i + j] += a[i] * b[j];
    }
  }

  // -- folding back
  int b_id = 0;
  int id;
  for (int k = 0; k < d - 1; k++)
  {
    // -- store buffer element temporarily
    id = 2 * d - 2 - k;
    b_id = buffer[id];
    buffer[id] = 0;

    // -- fold this deleted index back onto the buffer
    for (int l = 0; l < d; l++)
    {
      // -- index runs from id - (d-1) until id-1
      buffer.at(d - 2 + l - k) += b_id * this->reduction[l];
      // cout << b_id << " | " << this->reduction[l] << endl;
    }
  }

  vector<int> c;
  c.clear();
  for (int l = 0; l < d; l++)
  {
    c.push_back(buffer[l]);
  }

  return Ring(c, this->reduction);
}

Ring Ring::operator%(const int &m)
{
  // https://stackoverflow.com/a/44197900/7236657

  vector<int> mod;

  for (int i = 0; i < (this->dim); i++)
  {
    mod.push_back((m + ((this->representation)[i] % m)) % m);
  }

  return Ring(mod, this->reduction);
}

void Ring::operator=(const Ring other)
{
  this->representation = other.representation;
  this->reduction = other.reduction;
  this->dim = other.dim;
}

bool Ring::operator==(const Ring &other) const
{
  bool equality = true;

  equality &= (this->representation == other.representation);
  equality &= (this->reduction == other.reduction);
  equality &= (this->dim == other.dim);

  return equality;
}

string Ring::repr(void)
{
  string rep = "";
  int c;
  string x;
  for (vector<int>::size_type i = 0; i < this->representation.size(); i++)
  {
    c = this->representation[i];

    if (i == 0)
    {
      x = "";
    }
    else
    {
      x = "x^" + to_string(i);
    }

    if (c < 0)
    {
      rep += to_string(c) + x;
    }
    else if (c > 0)
    {
      if (i == 0)
      {
        rep += to_string(c) + x;
      }
      else
      {
        rep += "+" + to_string(c) + x;
      }
    }
  }

  if (rep == "")
  {
    rep = "0";
  }

  return rep;
}
