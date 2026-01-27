#ifndef MEMSYS_H
#define MEMSYS_H

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <math.h>

#define TOP -1
#define BOT 1
#define EMPTY -100
#define ALOT 100000000.0

#define max(a,b) \
  ({ __typeof__ (a) _a = (a); \
      __typeof__ (b) _b = (b); \
    _a > _b ? _a : _b; })

#define min(a,b) \
  ({ __typeof__ (a) _a = (a); \
      __typeof__ (b) _b = (b); \
    _a < _b ? _a : _b; })

#define sgn(a) \
  ({ __typeof__ (a) _a = (a); \
    _a < 0 ? -1 : 1; })
  
typedef struct vec{
  double x, y, z;
} vec;

typedef struct Meshpoint{
  vec Pos;
  vec nhat;
  double nNorm;
  double lc_x, lc_y;
  double CurvatureEnergy;
  double veff, veffdummy;
  double **ParticleProbs, *BondNums;
} MeshPoint;

typedef struct hist{
  double *CompHhist[2], *Hhist[2];
  double zmin, zmax, dz, histnorm;
  int binnums;
} hist;

typedef struct WangLandau{
  double **FreeEnergy, **Hitsz;
  double zmax, zmin, zmaxcut, zmincut, dz, eps;
  int    binnums;
  char   PotentialFile[1024];
} WangLandau;

typedef struct PopBias{
  int bnums, *pInds;
  double *Locs, *Strs;
  double *PBVals, *PBDummyVals;
} PopBias;


typedef struct meshbias{
  char   BiasType[1024];
  double BiasStr, VeffStr, theta0;

  //Current state
  double thetaTop, thetaBot;
  int imin, jmin;
  
  //Suggested state
  int idummy, jdummy;
  double thetadummy, thetadummy2;
} meshbias;  

typedef struct sys{
//MC Ingredients
  long MCSteps, CurrStep;
  long DoF, SweepIndex, SweepNum;
  double MDisp, GlobShiftSize, meshglobalrate ;
  gsl_rng *rng;
  
//Mesh properties
  int    Nx, Ny, MTot, twoMTot;
  double Lx, Ly, Lxon2, Lyon2, A;
  double kc, MeshEnergyPrefac;
  double dx, dy, dA;
  double dxsqr, dysqr;
  
  MeshPoint **TopMesh;
  MeshPoint **BotMesh;
  
  //Mesh MC move vars
    int i0, j0, choice, MeshShift;
    double suggestedz;
    //Local derivatives dummy vars
      vec lcCent, lcRight, lcLeft, lcIn, lcOut;
      vec nhatRight, nhatLeft, nhatIn, nhatOut;
      double nNormRight, nNormLeft, nNormIn, nNormOut;
      double veffRight, veffLeft, veffIn, veffOut, veffC;
      
      
  //MeshBiasParams
    meshbias   MeshBias;
    WangLandau WL;
    PopBias    PB;
    
//Particle properties
  int NumPTypes, *PTypePop;
  double *LPType[2], *BindPType, *RanPType;
  double *Fugacities[2], VeffZero, pDiam, VeffCGFact;
    
//Observables
  double AcceptanceRate, NodeAccept, ShiftAccept, ShiftAttempt;
  double TotEnergy, MeshEnergy, BiasEnergy, TotVeff;
  double TotCurPop[2], *CurParPop[2], *BondNums, TotAreaTop, TotAreaBot;
  
//Input/OutPut options
  int datarate, framerate;
  char Trajfile[1024], Datafile[1024], Restartfile[1024], Histfile[1024], SimInfo[1024];
  FILE *fileptr_traj, *fileptr_data, *fileptr_restart, *fileptr_input;
  
} sys;

#endif
