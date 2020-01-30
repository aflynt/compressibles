#include "ss.h"

int calc_isentropics(double gamma, double M)
{

  double popt = ppt(gamma, M);
  double qopt = qpt(gamma, M);
  double tott = ttt(gamma, M);
  double rort = rrt(gamma, M);
  double aoas = aoastar(gamma, M);

  printf("gama %6.4f ", gamma);
  printf("Mach %6.4f ", M);
  printf("popt %6.4f ", popt);
  printf("qopt %6.4f ", qopt);
  printf("tott %6.4f ", tott);
  printf("rort %6.4f ", rort);
  printf("aoas %6.4f\n", aoas);

  return 0;
}

int calc_normal_shock(double gamma, double M)
{
  double p2p1   = ns_p2p1(gamma, M);
  double r2r1   = ns_r2r1(gamma, M);
  double t2t1   = ns_t2t1(gamma, M);
  double pt2pt1 = ns_pt2pt1(gamma, M);
  double p1pt2  = ns_p1pt2(gamma, M);
  //double m2     = ns_m2(gamma, M);

  printf("gama %6.4f ", gamma);
  printf("Mach %6.4f ", M);
  printf("p2/p1 %6.4f ", p2p1);
  printf("r2/r1 %6.4f ", r2r1);
  printf("t2/t1 %6.4f ", t2t1);
  printf("pt2/pt1 %6.4f ", pt2pt1);
  printf("p1/pt2 %6.4f ", p1pt2);
  //printf("aoas %6.4f\n", aoas);
  putchar('\n');

  return 0;
}


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

