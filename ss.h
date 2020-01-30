#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#define PI 3.14159265359

#define CHECKPT {printf("Checkpoint: %s, line %d\n", __FILE__, __LINE__);\
  fflush(stdout);}

double rrt(double gamma, double M); //rho / rho_t
double ppt(double gamma, double M); //P / P_t
double qpt(double gamma, double M); //q / P_t
double ttt(double gamma, double M); //T / T_t 
double aoastar(double gamma, double M); //A / A* 

int getM(int argc, char *argv[], double *gamma, double *M);

