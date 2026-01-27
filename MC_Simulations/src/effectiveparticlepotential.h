#ifndef EFFPOT_H
#define EFFPOT_H

#include "membranesys.h"
#include "mesh.h"
#include <math.h>


double (*ReturnParticlePotentialLatticeSite)();
double ReturnConvParticlePotentialLatticeSite(sys *Sys, double dz, double sqrtnormbot, double sqrtnormtop);
double ReturnParticleLambPotentialLatticeSite(sys *Sys, double dz, double sqrtnormbot, double sqrtnormtop);

double ReturnTotalParticlePotential(sys *Sys);
double VeffChange(sys *Sys);
double VeffGlobalShiftChange(sys *Sys);
double ComputeVeffZero(sys *Sys);
double ReturnBondNumbers(double z, int ind, sys *Sys);

void SetVeffLatticeSite(sys *Sys, int i0, int j0);
void InitializeVeff(sys *Sys);
void StoreVeffChange(sys *Sys, double DeltaEVeff);
void StoreVeffGlobalShiftChange(sys *Sys, double DeltaEVeff);
void ComputeParticlePops(sys *Sys);

int ErrorCheck_Veff(sys *Sys);

#endif
