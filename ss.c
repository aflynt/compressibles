#include "ss.h"

void print_usage(void){
  printf("Alt-usage: ss -g 1.4 -m 2.0 \n");
  printf("Alt-usage: ss 2.0 \n");
}

double aoastar(double gam, double M)
{
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
  double exp = -gamma / (gamma - 1.0);
  double man = 1.0 + (gamma - 1.0) / 2.0 * M * M;
  return pow(man, exp);
}

double rrt(double gamma, double M)
{
  double exp = -1.0 / (gamma - 1.0);
  double man = 1.0 + (gamma - 1.0) / 2.0 * M * M;
  return pow(man, exp);
}

double qpt(double gamma, double M)
{
  double A = gamma * 0.5* M * M;
  double B = (1 + (gamma - 1) * 0.5 * M * M);
  double E = -gamma / (gamma - 1);
  return A * pow(B,E);
}

double ttt(double gamma, double M)
{
  return 1.0 / (1 + (gamma - 1.0) / 2.0 * M * M);
}

int getM(int argc, char *argv[], double *gamma, double *M)
{
  int i;
  int gotg=0;
  int gotm=0;
  char strg[] = "-g";
  char strm[] = "-m";
  char strh[] = "-h";

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

  return 0;
}
