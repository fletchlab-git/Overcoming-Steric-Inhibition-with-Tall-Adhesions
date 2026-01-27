#include "initialize.h"

int CountNumberOfCompletedRuns(char *OutputDir, char *tmp){
  int index=0;
  DIR * dirp;
  struct dirent * entry;
  char RestartDir[1024];
  
  strcpy(RestartDir, OutputDir);
  strcat(RestartDir, "/Restarts/");

//   printf("%s %s\n", RestartDir, OutputDir);
  dirp = opendir(RestartDir);
  if(dirp){
    while ((entry = readdir(dirp)) != NULL) {
      if (entry->d_type == DT_REG && strstr(entry->d_name, tmp) ){ 
        index++;
      }
    }
  }
  else{printf("\nERROR! Could not find %s on this machine\n", RestartDir);exit(1);}
  closedir(dirp);
  return index;
}

void SetOutPutFiles(sys *Sys, struct arguments *arguments ){
//Set filepaths
  char tmp[1024], tmp2[1024];
  if(0 == strcmp(arguments->Tag,"")){tmp[0] = '\0';}
  else{snprintf(tmp, 1024, "%s_", arguments->Tag);}
   
  if(strcmp(arguments->BiasType, "Pop")==0){
    int i=0;char tmp3[1024];strncpy(tmp3, arguments->PopLoc, 1024);
    while(tmp3[i]){
      if (isspace(tmp3[i]))
        tmp3[i]='_';
      i++;
    }
    
    snprintf(tmp2, 1024, "%sNodes_%d_NumPTypes_%d_kc_%.2lf_A_%.1lf_PopLoc_%s", 
             tmp,
             arguments->Nx*arguments->Ny,
             arguments->NumPType,
             arguments->kc,
             arguments->Lx*arguments->Ly,
             tmp3
            );
  }
  
  else{
    snprintf(tmp2, 1024, "%sNodes_%d_NumPTypes_%d_kc_%.2lf_A_%.1lf_theta0_%lf", 
             tmp,
             arguments->Nx*arguments->Ny,
             arguments->NumPType,
             arguments->kc,
             arguments->Lx*arguments->Ly,
             arguments->theta0
            );
  }
  
  strcpy(Sys->SimInfo, tmp2);
  
  int index = CountNumberOfCompletedRuns(arguments->OutputDir, tmp2);
  
  snprintf(Sys->Trajfile,    1024, "%s/Trajectories/Traj_%s_Run_%d_.xyz",
           arguments->OutputDir, tmp2, index);
  snprintf(Sys->Datafile,    1024, "%s/DataFiles/Data_%s_Run_%d_.dat",
           arguments->OutputDir, tmp2, index);
  snprintf(Sys->Histfile,    1024, "%s/DataFiles/Hists/Zhist_%s_Run_%d_.hist",
           arguments->OutputDir, tmp2, index);
  snprintf(Sys->Restartfile, 1024, "%s/Restarts/Restart_%s_Run_%d_.rest",
           arguments->OutputDir, tmp2, index);
  
  if (Sys->fileptr_traj!=NULL){fclose(Sys->fileptr_traj);}
  if (Sys->fileptr_data!=NULL){fclose(Sys->fileptr_data);}
  
  Sys->fileptr_traj = fopen(Sys->Trajfile, "w");
  Sys->fileptr_data = fopen(Sys->Datafile, "w");
  
}

void GetLengths(double *Llong, double *Lshort, sys *Sys){
  double L_curr;
  *Lshort = ALOT;*Llong = -ALOT;
  for(int i=0; i<Sys->NumPTypes; i++){
    if((Sys->Fugacities[0][i] + Sys->Fugacities[1][i])> 1e-6){
      L_curr = Sys->LPType[0][i] + Sys->LPType[1][i] + 0.5*Sys->RanPType[i];
      if(L_curr > *Llong  ) *Llong  = L_curr;
      if(L_curr < *Lshort ) *Lshort = L_curr;
    }
  }
}

void ShortInLong(double z0, double Llong, double Lshort, sys *Sys){
  double    f = (1.0 - z0/Llong)/(1.0 - Lshort/Llong);
  if( f>0 ) f = sqrt(f);
  else      f = 0.0;
  int Nmin = (int)round(0.5*Sys->Nx*(1. - f));
  int Nmax = (int)round(0.5*Sys->Nx*(1. + f));

  printf("Initializing geometry to: ShortInLong\n");
  
  for(int i=0; i<Sys->Nx; i++){
    for(int j=0; j<Sys->Ny; j++){
      if(i>Nmin && i<Nmax && j>Nmin && j<Nmax){
	Sys->TopMesh[i][j].Pos.z =  Lshort/2;
	Sys->BotMesh[i][j].Pos.z = -Lshort/2;
      }
      else{
	Sys->TopMesh[i][j].Pos.z =  Llong/2;
	Sys->BotMesh[i][j].Pos.z = -Llong/2;
      }      
    }
  }
}

void Slab(double z0, double Llong, double Lshort, sys *Sys){
  double f = (1.0 - z0/Llong)/(1.0 - Lshort/Llong);
  int Nmin = (int)round(0.5*Sys->Nx*(1.0 - f));
  int Nmax = (int)round(0.5*Sys->Nx*(1.0 + f));

  printf("Initializing geometry to: Slab\n");

  for(int i=0; i<Sys->Nx; i++){
    for(int j=0; j<Sys->Ny; j++){
      if(i>Nmin && i<Nmax){
	Sys->TopMesh[i][j].Pos.z =  Lshort/2;
	Sys->BotMesh[i][j].Pos.z = -Lshort/2;
      }
      else{
	Sys->TopMesh[i][j].Pos.z =  Llong/2;
	Sys->BotMesh[i][j].Pos.z = -Llong/2;
      }      
    }
  }
}

void LongInShort(double z0, double Llong, double Lshort, sys *Sys){
  double f = (z0 - Lshort)/(Llong - Lshort);f = sqrt(f);
  int Nmin = (int)round(0.5*Sys->Nx*(1. - f));
  int Nmax = (int)round(0.5*Sys->Nx*(1. + f));

  printf("Initializing geometry to: LongInShort\n");

  for(int i=0; i<Sys->Nx; i++){
   for(int j=0; j<Sys->Ny; j++){
      if(i>Nmin && i<Nmax && j>Nmin && j<Nmax){
	Sys->TopMesh[i][j].Pos.z =  Llong/2;
	Sys->BotMesh[i][j].Pos.z = -Llong/2;
      }
      else{
	Sys->TopMesh[i][j].Pos.z =  Lshort/2;
	Sys->BotMesh[i][j].Pos.z = -Lshort/2;
      }      
    }
  }
}

void MakeIntermediateConfig(double z0, double Lshort, sys *Sys){
  int Nmin = (int)round(0.25*Sys->Nx);
  int Nmax = (int)round(0.75*Sys->Nx);
  double Llong = (4.0*z0 - Lshort)/3.;
  
  for(int i=0; i<Sys->Nx; i++){
   for(int j=0; j<Sys->Ny; j++){
      if(i>Nmin && i<Nmax && j>Nmin && j<Nmax){
	Sys->TopMesh[i][j].Pos.z =  Lshort/2;
	Sys->BotMesh[i][j].Pos.z = -Lshort/2;
      }
      else{
	Sys->TopMesh[i][j].Pos.z =  Llong/2;
	Sys->BotMesh[i][j].Pos.z = -Llong/2;
      }      
    }
  }
} 

void MakeLBConfig(double z0, double Llong, sys *Sys){
  int Nmin = (int)round(0.25*Sys->Nx);
  int Nmax = (int)round(0.75*Sys->Nx);
  double LVlong = (4.0*z0 - Llong)/3.;
  
  for(int i=0; i<Sys->Nx; i++){
   for(int j=0; j<Sys->Ny; j++){
      if(i>Nmin && i<Nmax && j>Nmin && j<Nmax){
	Sys->TopMesh[i][j].Pos.z =  Llong/2;
	Sys->BotMesh[i][j].Pos.z = -Llong/2;
      }
      else{
	Sys->TopMesh[i][j].Pos.z =  LVlong/2;
	Sys->BotMesh[i][j].Pos.z = -LVlong/2;
      }      
    }
  }
}

void MakeShortInLong(double Llong, double Lshort, sys *Sys){
  int Nmin = (int)round(0.25*Sys->Nx);
  int Nmax = (int)round(0.75*Sys->Nx);
  for(int i=0; i<Sys->Nx; i++){
   for(int j=0; j<Sys->Ny; j++){
      if(i>Nmin && i<Nmax && j>Nmin && j<Nmax){
	Sys->TopMesh[i][j].Pos.z =  Lshort/2;
	Sys->BotMesh[i][j].Pos.z = -Lshort/2;
      }
      else{
	Sys->TopMesh[i][j].Pos.z =  Llong/2;
	Sys->BotMesh[i][j].Pos.z = -Llong/2;
      }      
    }
  }
}

void ModifyInitialConfig(sys *Sys, double z0){
  double Llong, Lshort, high, low;
  GetLengths(&Llong, &Lshort, Sys);
  high = Llong  - (Llong - Lshort)/M_PI;
  low  = Lshort + (Llong - Lshort)/M_PI;
  
//  printf("long: %lf short: %lf high: %lf low: %lf z0: %lf\n",
// 	  Llong, Lshort, high, low, z0
//        );exit(1);
  
//   if (z0 < Llong) MakeIntermediateConfig(z0, Lshort, Sys);
//   else            MakeLBConfig(z0, Llong, Sys);
  MakeShortInLong(Llong, Lshort, Sys);
  low = low;high = high;
  
//   if (z0 < Llong){
//     if      (z0 > high)   ShortInLong(z0, Llong, Lshort, Sys);
//     else if (z0 > low )   Slab(z0, Llong, Lshort, Sys);
//     else if (z0 > Lshort) LongInShort(z0, Llong, Lshort, Sys);
//   }
}

void CreateMeshes(sys *Sys, struct arguments *arguments){
  int i, j, k;
  double x;
  
  Sys->Nx      = arguments->Nx;
  Sys->Ny      = arguments->Ny;
  Sys->MTot    = arguments->Nx*arguments->Ny;
  Sys->twoMTot = 2*Sys->MTot;
    
  Sys->Lx    = arguments->Lx;
  Sys->Ly    = arguments->Ly;
  Sys->Lxon2 = arguments->Lx/2;
  Sys->Lyon2 = arguments->Ly/2;
  Sys->A     = arguments->Lx*arguments->Ly;
  
  Sys->dx    = Sys->Lx/Sys->Nx;
  Sys->dy    = Sys->Ly/Sys->Ny;
  Sys->dxsqr = Sys->dx*Sys->dx;
  Sys->dysqr = Sys->dy*Sys->dy;
  
  Sys->dA = Sys->dx*Sys->dy;
  Sys->kc = arguments->kc;
  Sys->MeshEnergyPrefac = (Sys->kc/2)*Sys->dA;
  
  Sys->TopMesh = (MeshPoint **)malloc(Sys->Nx*sizeof(MeshPoint *));
  Sys->BotMesh = (MeshPoint **)malloc(Sys->Nx*sizeof(MeshPoint *));
  
  for (i=0;i<Sys->Nx;i++){
    Sys->TopMesh[i]= (MeshPoint *)malloc(Sys->Ny*sizeof(MeshPoint));
    x=i*Sys->dx;
    for (j=0; j<Sys->Ny; j++){
      Sys->TopMesh[i][j].Pos.x            = x;
      Sys->TopMesh[i][j].Pos.y            = Sys->dy*j;
      Sys->TopMesh[i][j].Pos.z            = arguments->z0/2;
      Sys->TopMesh[i][j].ParticleProbs    = (double **)malloc((arguments->NumPType + 1)*sizeof(double *));
      for(k=0;k<(arguments->NumPType + 1);k++) Sys->TopMesh[i][j].ParticleProbs[k] = (double *)malloc((arguments->NumPType + 1)*sizeof(double));
      Sys->TopMesh[i][j].BondNums         = (double *)malloc((arguments->NumPType)*sizeof(double));
    }
  }
  
  for (i=0;i<Sys->Nx;i++){
    Sys->BotMesh[i]= (MeshPoint *)malloc(Sys->Ny*sizeof(MeshPoint));
    x=i*Sys->dx;
    for (j=0; j<Sys->Ny; j++){
      Sys->BotMesh[i][j].Pos.x            = x;
      Sys->BotMesh[i][j].Pos.y            = Sys->dy*j;
      Sys->BotMesh[i][j].Pos.z            = -arguments->z0/2;
      Sys->BotMesh[i][j].ParticleProbs    = (double **)malloc((arguments->NumPType + 1)*sizeof(double *));
      for(k=0;k<(arguments->NumPType + 1);k++) Sys->BotMesh[i][j].ParticleProbs[k] = (double *)malloc((arguments->NumPType + 1)*sizeof(double));
      Sys->BotMesh[i][j].BondNums         = (double *)malloc((arguments->NumPType)*sizeof(double));
    }
  }
  
  Sys->MDisp  = arguments->MDisp;
}

void CreateBiasPotential(sys *Sys, struct arguments *arguments){
  if (arguments->BiasStr<1e-6)   Sys->meshglobalrate = 0.0;
  else if (arguments->kc>1000.0) Sys->meshglobalrate = 1.0;
  else Sys->meshglobalrate = 1.0/Sys->twoMTot;

  UpdateBiases = E;
  strcpy(Sys->MeshBias.BiasType, arguments->BiasType);
  
  if(strcmp(arguments->BiasType, "Harmonic")==0 /*&& arguments->NumPType==0*/){
    Sys->MeshBias.theta0  = arguments->theta0/2;
    Sys->MeshBias.BiasStr = arguments->BiasStr;
    ComputeBiasPotential  = ComputeBiasPotentialHarmonic;
    MeshBiasChange        = MeshBiasChangeHarmonic;
    MeshBiasGlobalChange  = MeshBiasGlobalHarmonic;
    StoreBiasChange       = StoreBiasChangeHarmonic;
    PrintBiasPotential    = PrintBiasPotentialHarmonic;
    ErrorCheck_Bias       = ErrorCheckHarmonicBias;
    Sys->GlobShiftSize    = min(Sys->MDisp, sqrt(2.0/Sys->MeshBias.BiasStr));
  }
  
  else if(strcmp(arguments->BiasType, "Minimum")==0){
    Sys->MeshBias.theta0  = arguments->theta0;
    Sys->MeshBias.BiasStr = arguments->BiasStr;
    ComputeBiasPotential  = ComputeBiasPotentialMinDist;
    MeshBiasChange        = MeshBiasChangeMinDist;
    MeshBiasGlobalChange  = MeshBiasGlobalMinDist;
    StoreBiasChange       = StoreBiasChangeMinDist;
    PrintBiasPotential    = PrintBiasPotentialMinDist;
    ErrorCheck_Bias       = ErrorCheckMinDistBias;
    Sys->GlobShiftSize    = min(Sys->MDisp, sqrt(2.0/Sys->MeshBias.BiasStr));
  }
  
  else if(strcmp(arguments->BiasType, "Pop")==0 /*&& arguments->NumPType==0*/){
    printf("Population Bias not yet implemented!\nExiting...\n");exit(1);
    Sys->MeshBias.theta0    = 0.0;
    Sys->MeshBias.BiasStr   = 0.0001;
    ComputeBiasPotential    = ComputeBiasPotentialPop;
    MeshBiasChange          = MeshBiasChangePop;
    MeshBiasGlobalChange    = MeshBiasGlobalPop;
    StoreBiasChange         = StoreBiasChangePop;
    PrintBiasPotential      = PrintBiasPotentialPop;
    ErrorCheck_Bias         = ErrorCheckPopBias;
    InitializePopSampling(Sys, arguments);
    Sys->GlobShiftSize      = Sys->MDisp;    
  }
  
  else if(strcmp(arguments->BiasType, "WangLandau")==0){
    printf("Warning, screen output and restart file format not implemented for this bias type...\n");
    Sys->MeshBias.theta0    = arguments->theta0/2;
    Sys->MeshBias.BiasStr   = arguments->BiasStr;
    ComputeBiasPotential    = ComputeBiasPotentialWL;
    MeshBiasChange          = MeshBiasChangeWL;
    MeshBiasGlobalChange    = MeshBiasGlobalWL;
    StoreBiasChange         = StoreBiasChangeWL;
    PrintBiasPotential      = PrintBiasPotentialWL;
    ErrorCheck_Bias         = ErrorCheckWLBias;
    UpdateBiases            = UpdateWL;
    InitializeWLSampling(Sys, arguments);
    Sys->GlobShiftSize      = (Sys->WL.zmax - Sys->WL.zmin)/4.;
  }
  
  else if(strcmp(arguments->BiasType, "HarmonicVeff")){
    Sys->MeshBias.theta0  = arguments->theta0/2;
    Sys->MeshBias.BiasStr = arguments->BiasStr;
    Sys->MeshBias.VeffStr = 0.0*arguments->Nx*arguments->Ny;
    ComputeBiasPotential  = ComputeBiasPotentialHarmonicPlusVeff;
    MeshBiasChange        = MeshBiasChangeHarmonicPlusVeff;
    MeshBiasGlobalChange  = MeshBiasGlobalHarmonicPlusVeff;
    StoreBiasChange       = StoreBiasChangeHarmonicPlusVeff;
    PrintBiasPotential    = PrintBiasPotentialHarmonicPlusVeff;
    ErrorCheck_Bias       = ErrorCheckHarmonicPlusVeffBias;
    snprintf(arguments->BiasType, 1024, "HarmonicVeff");
    Sys->GlobShiftSize    = min(Sys->MDisp, sqrt(2.0/Sys->MeshBias.BiasStr));
  }
  
  else if (strcmp(arguments->BiasType, "Point")==0){
    Sys->MeshBias.theta0  = arguments->z0/2;
    Sys->MeshBias.BiasStr = arguments->BiasStr;
    ComputeBiasPotential  = ComputeBiasPotentialMinDist;
    MeshBiasChange        = MeshBiasChangeMinDist;
    StoreBiasChange       = StoreBiasChangeMinDist; 
    Sys->GlobShiftSize    = 0.4;
  }
  
  else if (strcmp(arguments->BiasType, "Linear")==0){
    Sys->MeshBias.theta0  = Sys->MTot*arguments->z0/2;  
    Sys->MeshBias.BiasStr = arguments->BiasStr/Sys->MTot;
    ComputeBiasPotential  = ComputeBiasPotentialAverage;
    MeshBiasChange        = MeshBiasChangeAverage;
    StoreBiasChange       = StoreBiasChangeAverage;
    PrintBiasPotential    = PrintBiasPotentialAverage;
    Sys->GlobShiftSize    = 0.4;
  }
  
  else{
    printf("Error: Could not parse bias type option!\nImplemented are: Point, Linear, and Harmonic\nYou Supplied: %s\n", arguments->BiasType);
    exit(1);
  }
}

void SetChemPotParams(sys *Sys, struct arguments *arguments){
  char Densities[1024], LPType[1024], BindPType[1024], RanPType[1024];
  int pflavors = Sys->NumPTypes = arguments->NumPType, i;
//   double alpha = 1e6*arguments->Nx*arguments->Ny/(arguments->Lx*arguments->Ly);
  double alpha = 1.e6/(arguments->pDiam*arguments->pDiam);
  Sys->VeffCGFact = Sys->dA/(arguments->pDiam*arguments->pDiam);
  
  ReturnParticlePotentialLatticeSite = ReturnConvParticlePotentialLatticeSite;
  
  if(pflavors>0){
    strcpy(Densities, arguments->Densities);
    strcpy(BindPType, arguments->BindPType);
    strcpy(RanPType,  arguments->RanPType);
    strcpy(LPType,    arguments->LPType);
    
    Sys->CurParPop[0]     = (double *)malloc((pflavors+1)*sizeof(double));
    Sys->CurParPop[1]     = (double *)malloc((pflavors+1)*sizeof(double));
    
    Sys->LPType[0]     = (double *)malloc(pflavors*sizeof(double));
    Sys->LPType[1]     = (double *)malloc(pflavors*sizeof(double));
    Sys->Fugacities[0] = (double *)malloc(pflavors*sizeof(double));
    Sys->Fugacities[1] = (double *)malloc(pflavors*sizeof(double)); 
    Sys->BindPType     = (double *)malloc(pflavors*sizeof(double));
    Sys->BondNums      = (double *)malloc(pflavors*sizeof(double));
    Sys->RanPType      = (double *)malloc(pflavors*sizeof(double));
  
    ExtractValuesFromStringList( (void *) Sys->BindPType,"DOUBLE", BindPType, pflavors);
    ExtractValuesFromStringList( (void *) Sys->RanPType, "DOUBLE", RanPType,  pflavors);
  
    double *tmp = (double *)malloc(2*pflavors*sizeof(double));
    ExtractValuesFromStringList( (void *) tmp, "DOUBLE", Densities,  2*pflavors);
    
    double zbar;
    zbar=0.0;
    for(i=0;i<pflavors;i++) zbar+=tmp[i];
    zbar = zbar/(alpha-zbar);
    for(i=0;i<pflavors;i++) Sys->Fugacities[0][i] = (1.0+zbar)*tmp[i]/alpha;
    
    zbar=0.0;
    for(i=pflavors;i<2*pflavors;i++) zbar+=tmp[i];
    zbar = zbar/(alpha-zbar);
    for(i=pflavors;i<2*pflavors;i++) Sys->Fugacities[1][i-pflavors] = (1.0+zbar)*tmp[i]/alpha;
    
    ExtractValuesFromStringList( (void *) tmp, "DOUBLE", LPType,  2*pflavors);
    for(i=0;i<pflavors;i++){Sys->LPType[0][i] = tmp[i];}
    for(i=pflavors;i<2*pflavors;i++){Sys->LPType[1][i-pflavors] = tmp[i];} 
  }
}

void SetMonteCarloParams(sys *Sys, struct arguments *arguments){
  Sys->SweepNum       = arguments->MCSweeps;
  Sys->DoF            = Sys->twoMTot;
  Sys->MCSteps        = (long)arguments->MCSweeps*Sys->DoF;
  
  if(arguments->datnums == 0)  Sys->datarate = Sys->SweepNum;
  else{Sys->datarate  = max(1, Sys->SweepNum/arguments->datnums);}
  
  if(arguments->framenums == 0) Sys->framerate = Sys->SweepNum;
  else{Sys->framerate = max(1, Sys->SweepNum/arguments->framenums);}
}

void InitializeStructs(sys *Sys, struct arguments *arguments){
  Sys->rng=gsl_rng_alloc(gsl_rng_taus2);
  gsl_rng_set(Sys->rng, arguments->RandSeed);

  if ( strcmp(arguments->Inputfile, "None")!=0 ){
    if ((Sys->fileptr_input = fopen(arguments->Inputfile, "r"))){}
    else{printf("ERROR:\n Could not find file: check if,\n%s\nis correct\n", arguments->Inputfile);exit(1);}
    char tmp[1024];strcpy(tmp, arguments->Inputfile);
    ReadInputCommands(Sys->fileptr_input, arguments);strcpy(arguments->Inputfile, tmp);
  }
  
  if      ( strcmp(arguments->BC, "PBC")==0 )   ChooseNode = ChooseNodePBC;
  else if ( strcmp(arguments->BC, "Frame")==0 ) ChooseNode = ChooseNodeFrame;
  
  CreateMeshes(Sys, arguments);
  Sys->fileptr_traj = Sys->fileptr_data = NULL;SetOutPutFiles(Sys, arguments);
  SetChemPotParams(Sys, arguments);
  CreateBiasPotential(Sys, arguments);
  SetMonteCarloParams(Sys, arguments);
  if( strcmp(arguments->Inputfile, "None")!=0 ){
    printf("Reading input file: %s\n\n", arguments->Inputfile);
    ReadConfiguration(Sys->fileptr_input, Sys);
    fclose(Sys->fileptr_input);
  }
  if(strcmp(arguments->ErrorCheck, "Yes")==0){
    ErrorCheckIntermittent = ErrorCheckConfiguration;
    ErrorCheckEvery = DontErrorCheck;
    printf("ErrorCheck Intermittently\n");
  }
  else if(strcmp(arguments->ErrorCheck, "YES")==0){
    ErrorCheckIntermittent = DontErrorCheck;
    ErrorCheckEvery = ErrorCheckConfiguration;
    printf("ErrorCheck Every step\n");
  }
  else {
    ErrorCheckIntermittent = DontErrorCheck;
    ErrorCheckEvery = DontErrorCheck;
  }
}
