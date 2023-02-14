#include <boost/algorithm/string.hpp> // include Boost, a C++ library
#include <iostream>

int main(){
 
  using namespace std;
  string target("Would you like a foo of chocolate. Two foos of chocolate?");
 
  cout << target << endl;
  boost::replace_all(target, "foo", "bar");
  cout << target << endl;

  return 0;
}
