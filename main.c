#include "ss.h"



int main(int argc, char *argv[])
{
  double gamma = 1.4;
  double M = 1.0;
  int choice = 0;

  choice = getM(argc, argv, &gamma, &M);

  if(choice == 0)
    calc_isentropics(gamma, M);
  else
    calc_normal_shock(gamma, M);

  return 0;
}

