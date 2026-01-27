#include <stdio.h>	
#include <stdlib.h>
#include <argp.h>
#include <math.h>

#include "parseparameters.h"
#include "initialize.h"
#include "membranesys.h"
#include "utilities.h"
#include "mesh.h"
#include "meshbiasfuncs.h"
#include "errorcheck.h"
#include "effectiveparticlepotential.h"

void ConductMC(sys *Sys, struct arguments *args){
//   Extract MC vars for shorter code
  long MCStepNums      = Sys->MCSteps; 
  long *ptr_SweepIndex = &(Sys->SweepIndex);
  long datarate        = Sys->datarate;
  long vidrate         = Sys->framerate;
  long printfrate      = max(1, Sys->SweepNum/500);
  
  Sys->AcceptanceRate = Sys->NodeAccept = Sys->ShiftAccept = Sys->ShiftAttempt = 0.0;
  
//  Helper vars
  long step;
  double acc;
//  Initialize all energies and lists
  InitializeMeshEnergies(Sys);
  Sys->TotEnergy=Sys->MeshEnergy + Sys->BiasEnergy + Sys->TotVeff;
  Sys->CurrStep = 0;
//  Perform MC simulation
  for (step=0;step<MCStepNums;step++){
    acc = PerformMeshMove(Sys);
    
    ErrorCheckEvery(Sys, step, acc);
    Sys->AcceptanceRate += acc;
    if(Sys->MeshShift == 1) Sys->ShiftAccept += acc;
    else Sys->NodeAccept += acc;
    Sys->ShiftAttempt += Sys->MeshShift;
      
    if(step % Sys->DoF == 0){
      if((*ptr_SweepIndex%printfrate==0 || *ptr_SweepIndex%datarate==0) && Sys->NumPTypes>0 ) ComputeParticlePops(Sys);
      ErrorCheckIntermittent(Sys, step, acc);
      if(*ptr_SweepIndex%printfrate==0) {PrintDataToScreen(stdout, Sys, step);WriteRestart(Sys, args);}
      if(*ptr_SweepIndex%datarate==0)    PrintData(Sys->fileptr_data, Sys, step);
      if(*ptr_SweepIndex%vidrate==0)     PrintConfiguration(Sys->fileptr_traj, Sys);
      (*ptr_SweepIndex)++;
    }
    Sys->CurrStep++;
  }
  
  if( Sys->NumPTypes>0 ) ComputeParticlePops(Sys);
  PrintDataToScreen(stdout, Sys, step);
  PrintData(Sys->fileptr_data, Sys, step);
  PrintConfiguration(Sys->fileptr_traj, Sys);
}

void Equilibrate(sys *Sys, struct arguments *arguments){
  int EquiTime  = 5000;
  double mrate  = Sys->meshglobalrate;
  int SweepNum  = Sys->SweepNum;
  long MCSteps  = Sys->MCSteps;
  
  Sys->meshglobalrate = 0.0;
  Sys->SweepNum       = EquiTime;
  Sys->MCSteps        = (long)Sys->DoF*EquiTime;
  
  ModifyInitialConfig(Sys, arguments->z0);
  
  SuggestMeshMove = SuggestMeshMoveEQ;
  
  printf("Equilibration run\n");
  ConductMC(Sys, arguments);  
  printf("Equilibration done!\n");
  SetOutPutFiles(Sys, arguments);
  
  Sys->meshglobalrate = mrate;
  Sys->SweepNum       = SweepNum;
  Sys->MCSteps        = MCSteps;
  
  Sys->SweepIndex     = 0;
  SuggestMeshMove     = SuggestMeshMoveR;
}

int main(int argc, char **argv){
  struct arguments arguments;
  sys Sys;
  InterpretInputCommands(&arguments, argc, argv);
  
  InitializeStructs(&Sys, &arguments);
  SuggestMeshMove = SuggestMeshMoveR;
    
  if (strcmp(arguments.Equil, "Yes")==0  ) Equilibrate(&Sys, &arguments);
  PrintInputCommands(Sys.fileptr_data, &arguments);
  ConductMC(&Sys, &arguments);
  
  PrintInputCommands(stdout, &arguments);
}
    
   


