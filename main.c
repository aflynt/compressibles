#include "ss.h"


int main(int argc, char *argv[])
{
  double gamma = 1.4;
  double M = 1.0;
  double popt;
  double qopt;
  double tott;
  double rort;
  double aoas;

  getM(argc, argv, &gamma, &M);


  popt = ppt(gamma, M);
  qopt = qpt(gamma, M);
  tott = ttt(gamma, M);
  rort = rrt(gamma, M);
  aoas = aoastar(gamma, M);

  printf("gama %6.4f ", gamma);
  printf("Mach %6.4f ", M);
  printf("popt %6.4f ", popt);
  printf("qopt %6.4f ", qopt);
  printf("tott %6.4f ", tott);
  printf("rort %6.4f ", rort);
  printf("aoas %6.4f\n", aoas);

  return 0;
}

