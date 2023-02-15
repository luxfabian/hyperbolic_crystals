#include <iostream>
#include "triangle_group.h"

int main(){

  using namespace std;

  TriangleGroup T = TriangleGroup(5,4);

  string word("ABABAAAABABBBBBB");

  cout << "Some generic word: " << word << endl;
  T.reduce(word);
  cout << "Its reduction in the {5,4,2} van Dyck triangle group: " << word << endl;

  

  G X=T.X;
  G Y=T.Y;
  G Z=T.Z;

  cout << "X is given by: " << endl;
  cout << X.repr() << endl;

  cout << "X^2 is given by: " << endl;
  cout << (X*X).repr() << endl;
  
  cout << "Y is given by: " << endl;
  cout << Y.repr() << endl;

  cout << "Y^2 is given by: " << endl;
  cout << (Y*Y).repr() << endl;

  cout << "Z is given by: " << endl;
  cout << Z.repr() << endl;

  cout << "Z^2 is given by: " << endl;
  cout << (Z*Z).repr() << endl;

  cout << "A is given by: " << endl;
  G A = X*Y;
  cout << A.repr() << endl;

  cout << "A^2 is given by: " << endl;
  cout << (A*A).repr() << endl;

  cout << "A^3 is given by: " << endl;
  cout << (A*A*A).repr() << endl;

  cout << "A^4 is given by: " << endl;
  cout << (A*A*A*A).repr() << endl;

  cout << "A^5 is given by: " << endl;
  cout << (A*A*A*A*A).repr() << endl;

  cout << "B is given by: " << endl;
  G B = Y*Z;
  cout << B.repr() << endl;

  cout << "B^2 is given by: " << endl;
  cout << (B*B).repr() << endl;

  cout << "B^3 is given by: " << endl;
  cout << (B*B*B).repr() << endl;

  cout << "B^4 is given by: " << endl;
  cout << (B*B*B*B).repr() << endl;

  cout << "C is given by: " << endl;
  G C = A*B;
  cout << C.repr() << endl;

  cout << "C^2 is given by: " << endl;
  cout << (C*C).repr() << endl;

  return 0;
}
