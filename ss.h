#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#define PI 3.14159265359

#define CHECKPT {printf("Checkpoint: %s, line %d\n", __FILE__, __LINE__);\
  fflush(stdout);}

// get command line arguments for mach and gamma
int getM(int argc, char *argv[], double *gamma, double *M);

// isentropic flows
double rrt(double gamma, double M); //rho / rho_t
double ppt(double gamma, double M); //P / P_t
double qpt(double gamma, double M); //q / P_t
double ttt(double gamma, double M); //T / T_t 
double aoastar(double gamma, double M); //A / A* 

// normal shocks
double ns_p2p1(double y, double M1); // P2 / P1
double ns_r2r1(double y, double M1); // rho2 / rho1
double ns_t2t1(double y, double M1); // T2 / T1
double ns_pt2pt1(double y, double M1); // pt2 / pt1
double ns_p1pt2(double y, double M1); // p1 / pt2




