#include "ss.h"

void print_usage(void)
{
  printf("Alt-usage: ss -g 1.4 -m 2.0 [takes gamma, mach]\n");
  printf("Alt-usage: ss 2.0 [takes mach]\n");
  printf("Alt-usage: ss -g 1.4 -m 2.0 -s [takes gamma, mach, normal shock]\n");
}

int getM(int argc, char *argv[], double *gamma, double *M)
{
  // get command line arguments for mach and gamma
  int i;
  int gotg=0;
  int gotm=0;
  int TYPE = 0;
  // TYPES:
  // 0 == isentropic flows
  // 1 == normal shocks
  char strg[] = "-g";  // gamma
  char strm[] = "-m";  // mach
  char strh[] = "-h";  // help
  char strs[] = "-s";  // shock

  // if only given 1 arg
  if(argc == 2){
    // given Mach number ?
    if(isdigit(argv[1][0])){
      *M = atof(argv[1]);
      gotm = 1;
    }
    // else print usage
    else {
      print_usage();
    }
  }
  // loop through all cmd line args
  else{
    for(i = 0; i < argc; i++){

      // found gamma
      if(strcmp(strg,argv[i]) == 0){
        *gamma = atof(argv[i+1]);
        gotg = 1;
      }

      // found M
      if(strcmp(strm,argv[i]) == 0){
        *M = atof(argv[i+1]);
        gotm = 1;
      }
      // usage
      if(strcmp(strh,argv[i]) == 0){
        print_usage();
      }
      // normal shock check
      if(strcmp(strs,argv[i]) == 0){
        TYPE = 1;
      }
    }
  }
  if(!gotg){
    //printf("Input Gamma:");
    //scanf("%lf",gamma);
    gotg = 1;
    *gamma = 1.4;
  }

  if(!gotm){
    printf("Input Mach:");
    scanf("%lf",M);
  }

  return TYPE;
}

// isentropic flows
double aoastar(double gam, double M)
{
  // compute A / A*
  double man1 = (gam + 1.0) / 2.0;
  double exp1 = -(gam + 1.0) / (2.0 * (gam - 1.0));
  double R1   = pow(man1, exp1);

  //double man2 = (1.0 + (gam - 1.0) / 2.0 * M * M);
  double man2 = 1.0 + M*M * (gam - 1.0) / 2.0 ;
  double exp2 = (gam + 1.0) / 2.0 / (gam - 1.0);
  double R2   = pow(man2, exp2);

  //{[1 + M^2 * (gam-1)/2]  EXP  [(gam+1)/(gam-1)/2]} 
  //MULT 
  //
  //{[(gam+1)/2]  EXP  -[(gam+1)/(gam-1)/2]} / M

  return (R1 * R2 / M);
}

double ppt(double gamma, double M)
{
  // compute p / P_t isentropic
  double exp = -gamma / (gamma - 1.0);
  double man = 1.0 + (gamma - 1.0) / 2.0 * M * M;
  return pow(man, exp);
}

double rrt(double gamma, double M)
{
  // compute rho / rho_t isentropic
  double exp = -1.0 / (gamma - 1.0);
  double man = 1.0 + (gamma - 1.0) / 2.0 * M * M;
  return pow(man, exp);
}

double qpt(double gamma, double M)
{
  // compute q / P_t isentropic
  double A = gamma * 0.5* M * M;
  double B = (1 + (gamma - 1) * 0.5 * M * M);
  double E = -gamma / (gamma - 1);
  return A * pow(B,E);
}

double ttt(double gamma, double M)
{
  // compute T / T_t isentropic
  return 1.0 / (1 + (gamma - 1.0) / 2.0 * M * M);
}

// normal shocks
double ns_p2p1(double y, double M1)
{
  // computes P2 / P1 across normal shock
  return (2 * y * M1 * M1 - (y - 1)) / (y + 1);
}

double ns_r2r1(double y, double M1)
{
  // computes rho2 / rho1 across normal shock
  return ((y+1)*M1*M1) / (2+(y-1)*M1*M1);
}

double ns_t2t1(double y, double M1)
{
  // computes T2 / T1 across normal
  // shock
  double A = 1.0 + 0.5*(y-1.0)*M1*M1;
  double B = 2.0*y / (y-1.0)*M1*M1 - 1.0;
  double C = (y+1)*(y+1)/(2.0*(y-1.0))*M1*M1;
  return A*B/C;
}

double ns_pt2pt1(double y, double M1)
{
  double A = (y+1) / (2*y*M1*M1-(y-1));
  double B = (y+1)*M1*M1 / ((y-1)*M1*M1+2.0);
  return pow(A,1/(y-1)) * pow(B, y/(y-1));
}

double ns_p1pt2(double y, double M1)
{
  double A = ppt(y, M1);
  double B = ns_pt2pt1(y,M1);
  return A / B;
}



