#ifndef UTILITIES_H
#define UTILITIES_H

#include "membranesys.h"
#include <string.h>

void PrintConfiguration(FILE *fp, sys *Sys);
void PrintData(FILE *fp, sys *Sys, long currentMCstep);
void PrintDataToScreen(FILE *fp, sys *Sys, long currentMCstep);
void PrintVec(vec target, char *note);

void ConfigurationRestart(FILE *fp, sys *Sys);
void ReadConfiguration(FILE *fp, sys *Sys);
int  computeindexdistsqr(int di, int dj, int Nx, int Ny);
void CopyVec(vec *target, vec source);
void ExtractValuesFromStringList( void *array, char *opt, char *str, int flavornums);

double ReturnHeightStd(sys *Sys);
void ReturnMeshHeights(double *zbot, double *ztop, sys *Sys);

#endif
