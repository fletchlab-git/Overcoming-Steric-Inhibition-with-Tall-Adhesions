#ifndef MESHBIAS_H
#define MESHBIAS_H

#include "parseparameters.h"
#include "membranesys.h"
#include "effectiveparticlepotential.h"
#include "utilities.h"

void   (*ComputeBiasPotential)();
double (*MeshBiasChange)();
double (*MeshBiasGlobalChange)();
void   (*StoreBiasChange)();
double (*PrintBiasPotential)();
int    (*ErrorCheck_Bias)();
void   (*UpdateBiases)();
void   E();

void   ComputeBiasPotentialMinDist(sys *Sys);
double MeshBiasChangeMinDist(sys *Sys);
double MeshBiasGlobalMinDist(sys *Sys);
void   StoreBiasChangeMinDist(sys *Sys, double DeltaEBias);
int    ErrorCheckMinDistBias(sys *Sys);

void   ComputeBiasPotentialAverage(sys *Sys);
double MeshBiasChangeAverage(sys *Sys);
void   StoreBiasChangeAverage(sys *Sys, double DeltaEBias);

void   ComputeBiasPotentialHarmonic(sys *Sys);
double MeshBiasChangeHarmonic(sys *Sys);
double MeshBiasGlobalHarmonic(sys *Sys);
void   StoreBiasChangeHarmonic(sys *Sys, double DeltaEBias);
int    ErrorCheckHarmonicBias(sys *Sys);

void   InitializePopSampling(sys *Sys, struct arguments *arguments);
void   ComputeBiasPotentialPop(sys *Sys);
double MeshBiasChangePop(sys *Sys);
double MeshBiasGlobalPop(sys *Sys);
void   StoreBiasChangePop(sys *Sys, double DeltaEBias);
int    ErrorCheckPopBias(sys *Sys);

void   InitializeWLSampling(sys *Sys, struct arguments *arguments);
void   ComputeBiasPotentialWL(sys *Sys);
double MeshBiasChangeWL(sys *Sys);
double MeshBiasGlobalWL(sys *Sys);
void   StoreBiasChangeWL(sys *Sys, double DeltaEBias);
int    ErrorCheckWLBias(sys *Sys);
void   UpdateWL(sys *Sys);

void   ComputeBiasPotentialHarmonicPlusVeff(sys *Sys);
double MeshBiasChangeHarmonicPlusVeff(sys *Sys);
double MeshBiasGlobalHarmonicPlusVeff(sys *Sys);
void   StoreBiasChangeHarmonicPlusVeff(sys *Sys, double DeltaEBias);
int    ErrorCheckHarmonicPlusVeffBias(sys *Sys);

double PrintBiasPotentialMinDist(sys *Sys);
double PrintBiasPotentialAverage(sys *Sys);
double PrintBiasPotentialWL(sys *Sys);
double PrintBiasPotentialHarmonic(sys *Sys);
double PrintBiasPotentialPop(sys *Sys);
double PrintBiasPotentialHarmonicPlusVeff(sys *Sys);

#endif
