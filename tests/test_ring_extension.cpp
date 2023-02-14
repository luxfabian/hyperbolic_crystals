#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

int main(){

  using namespace std;

  string line;
  ifstream file("ring_reduction.inp");

  const int p = 5;
  const int q = 4;
  const int n = 2*p*q;

  vector<int> reduction; 
  int buffer;
  int d=0;

  bool listening = false;

  if(file.is_open())
  {
    while(getline(file,line))
    {
      //--start analysing after keyword is found
      if(line.find("BEGIN") != string::npos)
      {
	listening = true;
      }

      //--stop analysing after keyword is found
      if(line.find("END") != string::npos)
      {
	break;
      }

      if(listening)
      {
        istringstream iss(line);
      
        if(iss>>buffer && iss>>d && buffer==n)
        {
          while(iss>>buffer)
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
    cout << "ohoh" << endl;
  }

  cout << "Succesfully read-in at index n=" << n << endl;
  cout << "Dimension of ring extension is d=" << d; 

  cout << "The coefficients are:" << endl;

  for(int c: reduction)
  {
    cout << c << endl;
  }

  return 0;
}
