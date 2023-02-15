#include <iostream>
#include "triangle_group.h"

int main(){

  using namespace std;
  
  TriangleGroup T = TriangleGroup(5,4);

  string word("ABABAAAABABBBBBB");

  cout << "Some generic word: " << word << endl;
  T.reduce(word);
  cout << "Its reduction in the {5,4,2} van Dyck triangle group: " << word << endl;

  cout << "A is given by: " << endl;

  G C=G(T.reduction);
  cout << T.A.repr() << endl;

  G A=T.A;

  C = A*A;

  cout << "A^p" << endl;

  cout << C.repr() << endl;
  return 0;
}
