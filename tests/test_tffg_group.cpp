#include <iostream>
#include "tffg_group.h"
#include "triangle_group.h"

int main(){

  using namespace std;

  long g=2; //genus
  int count = 0;

  TorsionFreeFuchsianGroup TFFG = TorsionFreeFuchsianGroup(g);

  for(G gamma: TFFG.gamma)
  {
    cout << "gamma_" << count++ + 1 << "=" << endl;
    cout << gamma.repr() << endl;
  }

  return 0;
}
