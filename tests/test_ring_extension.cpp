
#include "ring_extension.h"
#include <iostream>
#include <vector>

int main(){

  using namespace std;

  const int p=5;
  const int q=4;

  vector<int> reduction;

  reduction = read_numberfield_reduction(p,q);

  int d = reduction.size();

  Ring a(reduction);
  Ring b(reduction);

  a = Ring(reduction,reduction);
  b = a*a;

  cout << "a = " << a.repr() << endl;

  return 0;
}
