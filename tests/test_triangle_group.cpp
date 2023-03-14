#include <iostream>
#include "triangle_group.h"

int main(){

  using namespace std;

  TriangleGroup T = TriangleGroup(7,3);

  G X = T.X;
  G Y = T.Y;
  G Z = T.Z;

  cout << "X" << endl;
  cout << X.repr() << endl;

  cout << "X^2" << endl;
  cout << (X*X).repr() << endl;

  cout << "Y" << endl;
  cout << Y.repr() << endl;

  cout << "Y^2" << endl;
  cout << (Y*Y).repr() << endl;

  cout << "Z" << endl;
  cout << Z.repr() << endl;

  cout << "Z^2" << endl;
  cout << (Z*Z).repr() << endl;

  G A = X*Y;

  G B = Y*Z;

  G C = X*Z;
  
  // G A = T.A;
  // G A7 = T.A;

  cout << "A" << endl;
  cout << A.repr() << endl;

  cout << "A^7" << endl;
  cout << (A*A*A*A*A*A*A).repr() << endl;

  cout << "C" << endl;
  cout << C.repr() << endl;

  cout << "C^2" << endl;
  cout << (C*C).repr() << endl;

  cout << "B^3" << endl;
  cout << (B*B*B).repr() << endl;


  return 0;
}
