#ifndef SYSDEFS_H
#define SYSDEFS_H

#include <stdlib.h>
#include <argp.h>
#include "string.h"

struct arguments{
/* Struct which stores all input arguments. */

// MC properties
  int      MCSweeps;
  double   MDisp, z0;
  long int RandSeed;
  
// BiasProperties
  double BiasStr, theta0, WLzmin, WLzmax;
  int WLbins;
  char   BiasType[1024], PopStr[1024], PopLoc[1024], BC[1024];
    
// MeshProperties
  double kc, Lx, Ly;
  int    Nx, Ny;
  
// ParticleProperties
  int  NumPType;
  char PTypePop[1024], LPType[1024], BindPType[1024], RanPType[1024];
  double pDiam;
  
// OutPutProperties
  char OutputDir[1024], Tag[1024]; 
  int  framenums, datnums, Histbnums;
  
// GCMCProperties
  char Densities[2024];
  
// Restart/Equilibration option
  char Inputfile[1024], Equil[1024];
  
// ErrorCheckRun
  char ErrorCheck[1024];
};

static struct argp_option options[] = {
/* The options we understand.
First is the parameter name, then is _a single letter with single quotes_
for the code that is to be used in parse_opt, the variable type description,
and the help message. End it with an empty item {0}.
*/

  {"Lx"        , 'a', "DOUBLE", 0,  "Mesh length in x direction [nm]" },
  {"Ly"        , 'b', "DOUBLE", 0,  "Mesh length in y direction [nm]" },
  {"Nx"        , 'c', "INT"   , 0,  "Meshpoint nums in x direction" },
  {"Ny"        , 'd', "INT"   , 0,  "Meshpoint nums in y direction" },
  {"z0"        , 'e', "DOUBLE", 0,  "Initial membrane seperation [nm]" },
  {"MeshDisp"  , 'f', "DOUBLE", 0,  "Size of Mesh height suggestions [nm]" },
  {"kc"        , 'g', "DOUBLE", 0,  "Bending modulus [kBT]" },
  
  {"NumPType"  , 'h', "INT"   , 0, "Number of protein flavors" },
  {"LPType"    , 'i', "LIST"  , 0, "Length of each protein type" },
  {"BindPType" , 'j', "LIST"  , 0, "Binding strength of each protein type" },
  {"RanPType"  , 'k', "LIST"  , 0, "Range of binding of each protein type" },
  
  {"Densities" , 'l', "LIST"  , 0, "Densities of the particle resarvoir, in molecules per micron squared. The first [NumPType] entries are assigned to the bottom mesh, the second set to the top" },
  
  {"OutDir"    , 'm', "DIR"   , 0, "Output directory" },
  {"Tag"       , 'n', "STR"   , 0, "Cosmetic tag to add to filenames" },
  
  {"framenums" , 'o', "INT"   , 0, "Number of frames in video" },
  {"datnums"   , 'p', "INT"   , 0, "Number of datapoints in datafile"},
  {"histbnums" , 'q', "INT"   , 0, "Number of bins in dz histogram"},
  
  {"Input"     , 'r', "FILE"  , 0, "Input filename" },
  {"Equil"     , 's', "STR"   , 0, "Equilibration option [YES/NO]"},
  
  {"MCSweeps"  , 't', "INT"   , 0, "Number of MC sweeps to perform" },
  
  {"BiasStr"   , 'u', "DOUBLE", 0, "Strength of bias applied on mesh distance [kBT]"},
  {"BiasType"  , 'v', "STR"   , 0, "Type of Bias to apply on membranes: Harmonic, Minimum, Pop, WangLandau, or Thermodynamic integration"},
  {"BiasCenter", 'w', "DOUBLE", 0, "Shift parameter of bias potential"},
  {"WLbinnums" , 'x', "DOUBLE", 0, "Number of WL bins"},
  {"WLzmin"    , 'y', "DOUBLE", 0, "Beginning of WL sampling range"},
  {"WLzmax"    , 'z', "DOUBLE", 0, "End of WL sampling range"},
  
  {"PopStr", 'A', "STR"   , 0, "Spring constants for bond biases"},
  {"PopLoc", 'B', "STR"   , 0, "Locations of bond biases"},
  
  {"ErrorCheck", 'C', "STR"   , 0, "Check all configurations"},
  {"SEED"      , 'D', "LONG INT", 0, "Random seed for rng"},
  {"BC"        , 'E', "STR", 0, "Either PBC or Frame"},
  {"pDiam"     , 'F', "DOUBLE", 0, "Particle diameters"},
  {0}
};

void InterpretInputCommands( struct arguments *arguments, int argc, char **argv);
void PrintInputCommands(FILE *fp, struct arguments *arguments);
void ReadInputCommands(FILE *fp, struct arguments *arguments);

#endif
