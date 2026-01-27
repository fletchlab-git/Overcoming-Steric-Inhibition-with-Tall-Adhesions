#include "meshbiasfuncs.h"

void E(){}

//MinDist Sampling

int  ErrorCheckMinDistBias(sys *Sys){
  int Err = 0, i, j, idummy = -1, jdummy = -1;
  double zdiff = ALOT, tmpzdiff, zsumAve = 0.0;
  
  for(i=0; i<Sys->Nx;i++){
    for(j=0;j<Sys->Ny;j++){
      zsumAve += (Sys->TopMesh[i][j].Pos.z + Sys->BotMesh[i][j].Pos.z)/Sys->twoMTot;
      tmpzdiff = Sys->TopMesh[i][j].Pos.z - Sys->BotMesh[i][j].Pos.z;
      if(tmpzdiff < zdiff){
        zdiff = tmpzdiff; idummy = i; jdummy = j;
      }
    }
  }
  
  if(idummy != Sys->MeshBias.imin || jdummy != Sys->MeshBias.jmin){
    printf("\nStored: imin %d jmin %d zdiff %f\nfound: imin %d jmin %d zdiff %f\n", Sys->MeshBias.imin, Sys->MeshBias.jmin, Sys->MeshBias.thetaTop, idummy, jdummy, zdiff);
    Err=6;
  }
  
  if(fabs(zdiff - Sys->MeshBias.thetaTop)>1e-6){
    printf("\nzdiffMin is %lf stored %lf\ntheta0: %lf\n", zdiff, Sys->MeshBias.thetaTop, Sys->MeshBias.theta0);
    Err=6; 
  }
  
  if(fabs(zsumAve - Sys->MeshBias.thetaBot)>1e-6){
    printf("\nCenter of mass is %lf stored %lf\n", zsumAve, Sys->MeshBias.thetaBot);
    Err=6; 
  }
  
  double En = PrintBiasPotentialMinDist(Sys);
  
  if(fabs(En-Sys->BiasEnergy)>1e-6){
    printf("Bias energy off computed %lf stored %lf delta: %e\n", En, Sys->BiasEnergy, fabs(En-Sys->BiasEnergy));
    Err=6;
  }
  
  if(Err!=0){
    printf("\nSummary\nthetaTop: %lf zdiff: %lf\nthatBot: %lf zsumAve: %lf\nimin %d jmin %d i %d j %d\nE stored %lf Ecomp %lf\n\n", Sys->MeshBias.thetaTop, zdiff, Sys->MeshBias.thetaBot, zsumAve, 
           Sys->MeshBias.imin, Sys->MeshBias.jmin, idummy, jdummy,
           Sys->BiasEnergy, En 
    );
  }
  return Err;
}

double PrintBiasPotentialMinDist(sys *Sys){
  double zdiff = ALOT, tmpzdiff, zsumAve = 0.0;
  for(int i=0; i<Sys->Nx;i++){
    for(int j=0;j<Sys->Ny;j++){
      zsumAve += (Sys->TopMesh[i][j].Pos.z + Sys->BotMesh[i][j].Pos.z)/Sys->twoMTot;
      tmpzdiff = Sys->TopMesh[i][j].Pos.z - Sys->BotMesh[i][j].Pos.z;
      if(tmpzdiff < zdiff) zdiff = tmpzdiff;
    }
  }
  return zsumAve*zsumAve + Sys->MeshBias.BiasStr*(zdiff - Sys->MeshBias.theta0)*(zdiff - Sys->MeshBias.theta0);
}

void ComputeBiasPotentialMinDist(sys *Sys){
  double zsumAve = 0.0, minz0 = ALOT, dummy;
  for(int i=0; i<Sys->Nx;i++){
    for(int j=0; j<Sys->Ny;j++){
      zsumAve += (Sys->TopMesh[i][j].Pos.z + Sys->BotMesh[i][j].Pos.z)/Sys->twoMTot;
      dummy = Sys->TopMesh[i][j].Pos.z - Sys->BotMesh[i][j].Pos.z;
      if (dummy<minz0){
        Sys->MeshBias.imin=i;Sys->MeshBias.jmin=j;
        minz0=dummy;
        Sys->MeshBias.thetaTop=minz0;
      }
    }
  }
  Sys->MeshBias.thetaBot = zsumAve;
  Sys->BiasEnergy= Sys->MeshBias.thetaBot*Sys->MeshBias.thetaBot + Sys->MeshBias.BiasStr*(Sys->MeshBias.thetaTop - Sys->MeshBias.theta0)*(Sys->MeshBias.thetaTop - Sys->MeshBias.theta0);
}

double MeshBiasChangeMinDist(sys *Sys){
  double dummy, deltaEBias;
  int *idummy = &(Sys->MeshBias.idummy), i0=Sys->i0;
  int *jdummy = &(Sys->MeshBias.jdummy), j0=Sys->j0;
  double *thetadummy = &(Sys->MeshBias.thetadummy);
  
  double dzSuggested = (Sys->choice == TOP) ? 
                       Sys->suggestedz - Sys->BotMesh[i0][j0].Pos.z: 
                       Sys->TopMesh[i0][j0].Pos.z - Sys->suggestedz;
  
  if(i0 == Sys->MeshBias.imin && j0 == Sys->MeshBias.jmin){
    *thetadummy=dzSuggested;*idummy=i0;*jdummy=j0;
    for(int i=0; i<Sys->Nx;i++){
      for(int j=0; j<Sys->Ny;j++){
        dummy = Sys->TopMesh[i][j].Pos.z - Sys->BotMesh[i][j].Pos.z;
        if (dummy<*thetadummy && (i!=i0 || j!=j0) ){
          *idummy=i;*jdummy=j;
          *thetadummy=dummy;
        }
      }
    }
    deltaEBias = Sys->MeshBias.BiasStr*( (*thetadummy - Sys->MeshBias.theta0)*(*thetadummy - Sys->MeshBias.theta0) - (Sys->MeshBias.thetaTop- Sys->MeshBias.theta0)*(Sys->MeshBias.thetaTop - Sys->MeshBias.theta0) );
  }
  
  else{
    dummy=dzSuggested;
    if( dummy < Sys->MeshBias.thetaTop ){
      *idummy = i0;*jdummy = j0;*thetadummy = dummy;
      deltaEBias = Sys->MeshBias.BiasStr*( (*thetadummy - Sys->MeshBias.theta0)*(*thetadummy - Sys->MeshBias.theta0) - (Sys->MeshBias.thetaTop - Sys->MeshBias.theta0)*(Sys->MeshBias.thetaTop - Sys->MeshBias.theta0));
    }
    else{
      deltaEBias=0.0;
      *idummy=Sys->MeshBias.imin;*jdummy=Sys->MeshBias.jmin;
      *thetadummy=Sys->MeshBias.thetaTop;
    }
  }
  
  dzSuggested = (Sys->choice==TOP) ? (Sys->suggestedz - Sys->TopMesh[i0][j0].Pos.z) :
(Sys->suggestedz - Sys->BotMesh[i0][j0].Pos.z); 
  dzSuggested /= Sys->twoMTot;
  Sys->MeshBias.thetadummy2 = Sys->MeshBias.thetaBot + dzSuggested;

  return deltaEBias + dzSuggested*(dzSuggested + 2.0*Sys->MeshBias.thetaBot);
}

double MeshBiasGlobalMinDist(sys *Sys){
  double prefac = ((Sys->choice==BOT) ? -1. : 1.), dz = prefac*Sys->suggestedz;
  Sys->MeshBias.thetadummy  = Sys->MeshBias.thetaTop + dz;
  Sys->MeshBias.thetadummy2 = Sys->MeshBias.thetaBot + Sys->suggestedz/2;
  Sys->MeshBias.idummy      = Sys->MeshBias.imin;
  Sys->MeshBias.jdummy      = Sys->MeshBias.jmin;
  return   Sys->MeshBias.BiasStr*dz*( 2.*(Sys->MeshBias.thetaTop - Sys->MeshBias.theta0) + dz )
         + 0.5*Sys->suggestedz*(2.*Sys->MeshBias.thetaBot + 0.5*Sys->suggestedz);
}

void StoreBiasChangeMinDist(sys *Sys, double DeltaEBias){
  Sys->MeshBias.imin      = Sys->MeshBias.idummy;
  Sys->MeshBias.jmin      = Sys->MeshBias.jdummy;
  Sys->MeshBias.thetaTop  = Sys->MeshBias.thetadummy;
  Sys->MeshBias.thetaBot  = Sys->MeshBias.thetadummy2;
  Sys->BiasEnergy        += DeltaEBias;
}

//Constant Force Sampling

void ComputeBiasPotentialAverage(sys *Sys){
  double zave=0.0;
  for(int i = 0;i<Sys->Nx;i++){
    for(int j = 0;j<Sys->Ny;j++){
      zave+=(Sys->TopMesh[i][j].Pos.z - Sys->BotMesh[i][j].Pos.z);
    }
  }
  Sys->MeshBias.thetaTop = zave;
  Sys->BiasEnergy=Sys->MeshBias.BiasStr*(zave - Sys->MeshBias.theta0);
}

double PrintBiasPotentialAverage(sys *Sys){
  double zave=0.0;
  for(int i = 0;i<Sys->Nx;i++){
    for(int j = 0;j<Sys->Ny;j++){
      zave+=(Sys->TopMesh[i][j].Pos.z - Sys->BotMesh[i][j].Pos.z);
    }
  }
  return Sys->MeshBias.BiasStr*(zave - Sys->MeshBias.theta0);
}

double MeshBiasChangeAverage(sys *Sys){
  double dzSuggested = (Sys->choice == TOP) ? 
                       Sys->suggestedz - Sys->TopMesh[Sys->i0][Sys->j0].Pos.z: 
                       Sys->BotMesh[Sys->i0][Sys->j0].Pos.z - Sys->suggestedz;

  Sys->MeshBias.thetadummy = Sys->MeshBias.thetaTop + dzSuggested;
  return Sys->MeshBias.BiasStr*dzSuggested;      
}

void StoreBiasChangeAverage(sys *Sys, double DeltaEBias){
  Sys->MeshBias.thetaTop= Sys->MeshBias.thetadummy;
  Sys->BiasEnergy+=DeltaEBias;
}

//Harmonic Bias Sampling routines

int  ErrorCheckHarmonicBias(sys *Sys){
  int Err=0;
  double ztop=0.0, zbot=0.0;
  
  ReturnMeshHeights(&zbot, &ztop, Sys);
  ztop-=Sys->MeshBias.theta0;zbot+=Sys->MeshBias.theta0;
  
  if(fabs(ztop-Sys->MeshBias.thetaTop)>1e-6){
    printf("displacement of top is %lf stored %lf\ntheta0: %lf\n", ztop, Sys->MeshBias.thetaTop, Sys->MeshBias.theta0);
    Err=6; 
  }
  if(fabs(zbot-Sys->MeshBias.thetaBot)>1e-6){
    printf("displacement of bot is %lf stored %lf\ntheta0: %lf\n", zbot, Sys->MeshBias.thetaBot,-Sys->MeshBias.theta0);
    Err=6; 
  }
  double En = Sys->MeshBias.BiasStr*(zbot*zbot + ztop*ztop);
  if(fabs(En-Sys->BiasEnergy)>1e-6){
    printf("Bias energy off computed %lf stored %lf\n", En, Sys->BiasEnergy);
    Err=6;
  }
  return Err;
}

double PrintBiasPotentialHarmonic(sys *Sys){
  double ztop=0.0, zbot=0.0;
  ReturnMeshHeights(&zbot, &ztop, Sys);
  zbot+=Sys->MeshBias.theta0;ztop-=Sys->MeshBias.theta0;
  return Sys->MeshBias.BiasStr*(zbot*zbot + ztop*ztop);
}

void ComputeBiasPotentialHarmonic(sys *Sys){
  double ztop=0.0, zbot=0.0;
  ReturnMeshHeights(&zbot, &ztop, Sys);
  zbot+=Sys->MeshBias.theta0;ztop-=Sys->MeshBias.theta0;
  Sys->MeshBias.thetaBot = zbot;Sys->MeshBias.thetaTop = ztop;
  Sys->BiasEnergy        = Sys->MeshBias.BiasStr*(zbot*zbot + ztop*ztop);
}

double MeshBiasChangeHarmonic(sys *Sys){
  double dzSuggested;
  if(Sys->choice == TOP){
    dzSuggested = (Sys->suggestedz - Sys->TopMesh[Sys->i0][Sys->j0].Pos.z)/Sys->MTot;
    Sys->MeshBias.thetadummy = Sys->MeshBias.thetaTop + dzSuggested;
    return Sys->MeshBias.BiasStr*dzSuggested*(dzSuggested + 2.0*Sys->MeshBias.thetaTop);
  }
  else{
    dzSuggested = (Sys->suggestedz - Sys->BotMesh[Sys->i0][Sys->j0].Pos.z)/Sys->MTot;
    Sys->MeshBias.thetadummy = Sys->MeshBias.thetaBot + dzSuggested;
    return Sys->MeshBias.BiasStr*dzSuggested*(dzSuggested + 2.0*Sys->MeshBias.thetaBot);
  }
}

double MeshBiasGlobalHarmonic(sys *Sys){
  if(Sys->choice == BOT){
    Sys->MeshBias.thetadummy = Sys->MeshBias.thetaBot + Sys->suggestedz;
    return Sys->MeshBias.BiasStr*Sys->suggestedz*(Sys->suggestedz + 2.0*Sys->MeshBias.thetaBot);
  }
  else{
    Sys->MeshBias.thetadummy = Sys->MeshBias.thetaTop + Sys->suggestedz;
    return Sys->MeshBias.BiasStr*Sys->suggestedz*(Sys->suggestedz + 2.0*Sys->MeshBias.thetaTop);
  }  
}

void StoreBiasChangeHarmonic(sys *Sys, double DeltaEBias){
  if (Sys->choice == BOT) Sys->MeshBias.thetaBot = Sys->MeshBias.thetadummy;
  else                   {Sys->MeshBias.thetaTop = Sys->MeshBias.thetadummy;}
                          Sys->BiasEnergy       += DeltaEBias;
}

//Population Bias Sampling routines

void InitializePopSampling(sys *Sys, struct arguments *arguments){
  int i, bnums = 0;
  char Locs[1024], Strs[1024];
  
  for (i=0;i<Sys->NumPTypes;i++){
    if(Sys->BindPType[i] > 1e-6) bnums++;
  }
  
  Sys->PB.bnums       = bnums;
  Sys->PB.Locs        = (double *)malloc(bnums*sizeof(double));
  Sys->PB.Strs        = (double *)malloc(bnums*sizeof(double));
  Sys->PB.PBVals      = (double *)malloc(bnums*sizeof(double));
  Sys->PB.PBDummyVals = (double *)malloc(bnums*sizeof(double));
  Sys->PB.pInds       = (int    *)malloc(bnums*sizeof(int));
  
  int g=0;
  for (i=0;i<Sys->NumPTypes;i++){
    if(Sys->BindPType[i] > 1e-6){
      Sys->PB.pInds[g] = i;g++;
    }
  }
  
  strcpy(Locs, arguments->PopLoc);strcpy(Strs, arguments->PopStr);
  ExtractValuesFromStringList( (void *) Sys->PB.Locs,"DOUBLE", Locs, bnums);
  ExtractValuesFromStringList( (void *) Sys->PB.Strs,"DOUBLE", Strs, bnums);
}

int ErrorCheckPopBias(sys *Sys){
  int Err=0;
  double ztop=0.0, zbot=0.0;
  for(int i = 0;i<Sys->Nx;i++){
    for(int j = 0;j<Sys->Ny;j++){
      ztop += Sys->TopMesh[i][j].Pos.z/Sys->MTot;
      zbot += Sys->BotMesh[i][j].Pos.z/Sys->MTot;
    }
  }
  ztop-=Sys->MeshBias.theta0;
  zbot+=Sys->MeshBias.theta0;
  if(fabs(ztop-Sys->MeshBias.thetaTop)>1e-6){
    printf("displacement of top is %lf stored %lf\ntheta0: %lf\n", ztop, Sys->MeshBias.thetaTop, Sys->MeshBias.theta0);
    Err=6; 
  }
  if(fabs(zbot-Sys->MeshBias.thetaBot)>1e-6){
    printf("displacement of bot is %lf stored %lf\ntheta0: %lf\n", zbot, Sys->MeshBias.thetaBot,-Sys->MeshBias.theta0);
    Err=6; 
  }
  double En = Sys->MeshBias.BiasStr*(zbot*zbot + ztop*ztop);
  if(fabs(En-Sys->BiasEnergy)>1e-6){
    printf("Bias energy off computed %lf stored %lf\n", En, Sys->BiasEnergy);
    Err=6;
  }
  return Err;
}

double PrintBiasPotentialPop(sys *Sys){
  double ztop=0.0, zbot=0.0;
  for(int i = 0;i<Sys->Nx;i++){
    for(int j = 0;j<Sys->Ny;j++){
      ztop += Sys->TopMesh[i][j].Pos.z/Sys->MTot;
      zbot += Sys->BotMesh[i][j].Pos.z/Sys->MTot;
    }
  }
  ztop-=Sys->MeshBias.theta0;
  zbot+=Sys->MeshBias.theta0;
  return Sys->MeshBias.BiasStr*(zbot*zbot + ztop*ztop);
}

void ComputeBiasPotentialPop(sys *Sys){
 double ztop=0.0, zbot=0.0;
  for(int i = 0;i<Sys->Nx;i++){
    for(int j = 0;j<Sys->Ny;j++){
      ztop += Sys->TopMesh[i][j].Pos.z/Sys->MTot;
      zbot += Sys->BotMesh[i][j].Pos.z/Sys->MTot;
    }
  }
  ztop-=Sys->MeshBias.theta0;
  zbot+=Sys->MeshBias.theta0;
  
  Sys->MeshBias.thetaTop = ztop;
  Sys->MeshBias.thetaBot = zbot;
  Sys->BiasEnergy        = Sys->MeshBias.BiasStr*(zbot*zbot + ztop*ztop);
}

double MeshBiasChangePop(sys *Sys){
  double dzSuggested;
  if(Sys->choice == TOP){
    dzSuggested = (Sys->suggestedz - Sys->TopMesh[Sys->i0][Sys->j0].Pos.z)/Sys->MTot;
    Sys->MeshBias.thetadummy = Sys->MeshBias.thetaTop + dzSuggested;
    return Sys->MeshBias.BiasStr*dzSuggested*(dzSuggested + 2.0*Sys->MeshBias.thetaTop);
  }
  else{
    dzSuggested = (Sys->suggestedz - Sys->BotMesh[Sys->i0][Sys->j0].Pos.z)/Sys->MTot;
    Sys->MeshBias.thetadummy = Sys->MeshBias.thetaBot + dzSuggested;
    return Sys->MeshBias.BiasStr*dzSuggested*(dzSuggested + 2.0*Sys->MeshBias.thetaBot);
  }
}

double MeshBiasGlobalPop(sys *Sys){
  if(Sys->choice == BOT){
    Sys->MeshBias.thetadummy = Sys->MeshBias.thetaBot + Sys->suggestedz;
    return Sys->MeshBias.BiasStr*Sys->suggestedz*(Sys->suggestedz + 2.0*Sys->MeshBias.thetaBot);
  }
  else{
    Sys->MeshBias.thetadummy = Sys->MeshBias.thetaTop + Sys->suggestedz;
    return Sys->MeshBias.BiasStr*Sys->suggestedz*(Sys->suggestedz + 2.0*Sys->MeshBias.thetaTop);
  }
}

void StoreBiasChangePop(sys *Sys, double DeltaEBias){
  if (Sys->choice == BOT) Sys->MeshBias.thetaBot = Sys->MeshBias.thetadummy;
  else                   {Sys->MeshBias.thetaTop = Sys->MeshBias.thetadummy;}
                          Sys->BiasEnergy       += DeltaEBias;
}

//Wang Landau sampling

int Return_WL_Index(double Theta, sys *Sys){ 
  return  max(0, min(Sys->WL.binnums-1, (int)round( (Theta - Sys->WL.zmin)/Sys->WL.dz)));
}

double WLEnergy(double Theta, sys *Sys){
  return Sys->WL.FreeEnergy[1][Return_WL_Index(Theta, Sys)];
}

double BiasPot(double ztop, double zbot, double z0, double Amp){
  return Amp*(ztop + zbot)*(ztop + zbot);
}

void AllocateWLArrays(sys *Sys){
  Sys->WL.FreeEnergy = (double **)malloc(2*sizeof(double *));
  Sys->WL.Hitsz      = (double **)malloc(2*sizeof(double *));
    
  for(int i=0;i<2;i++){
    Sys->WL.FreeEnergy[i] = (double *)malloc(Sys->WL.binnums*sizeof(double));
    Sys->WL.Hitsz[i]      = (double *)malloc(Sys->WL.binnums*sizeof(double));
  }

  Sys->WL.dz      = (Sys->WL.zmax - Sys->WL.zmin)/(Sys->WL.binnums-1);
  Sys->WL.zmincut =  Sys->WL.zmin - Sys->WL.dz/2;
  Sys->WL.zmaxcut =  Sys->WL.zmax + Sys->WL.dz/2;

  for(int i=0;i<Sys->WL.binnums;i++){
    Sys->WL.FreeEnergy[0][i] =  Sys->WL.zmin + i*Sys->WL.dz;
    Sys->WL.Hitsz[0][i]      =  Sys->WL.zmin + i*Sys->WL.dz;
    Sys->WL.Hitsz[1][i]      = 0.0;
  }
}

void SetWLParamsFromFile(sys *Sys, FILE *fp){
  char line[1024], *dummy;
  int bnum = 1;
  double zmax, d1, d2;

  dummy = fgets(line, 1024, fp);
  sscanf(line, "#z(%*d) FreeEnergy(%*d) Probs(%*d) current eps: %lf min: %*s\n", &(Sys->WL.eps));
  dummy=dummy;
 
  dummy = fgets(line, 1024, fp);
  sscanf(line, "%lf %lf %lf\n", &(Sys->WL.zmin), &d1, &d2);
 
  while( fgets(line, 1024, fp) !=NULL){
    sscanf(line, "%lf %lf %lf\n", &zmax, &d1, &d2);bnum++;
  }

  Sys->WL.binnums = bnum;
  Sys->WL.zmax    = zmax;
}

int ReadWLFile(sys *Sys, char *fname){
  int ind, res = 0;
  FILE *fp = fopen(fname, "r");
  char line[1024], *dummy;
  double zval, fval, d;
  if (fp != NULL){
    printf("Found WL restart file\n");res=1;
    SetWLParamsFromFile(Sys, fp);
    AllocateWLArrays(Sys);    

    fclose(fp);fp = fopen(fname, "r");
    dummy = fgets(line, 1024, fp);dummy=dummy;
    while( fgets(line, 1024, fp) !=NULL){
      sscanf(line, "%lf %lf %lf\n", &zval, &fval, &d);
      ind = Return_WL_Index(zval, Sys);Sys->WL.FreeEnergy[1][ind] = fval;
    }
    fclose(fp);
  }
  printf("checking for file: %s\n result: %d\n", fname, res);

  return res;
}

void InitializeWLSampling(sys *Sys, struct arguments *arguments){
  int res = 0;
  snprintf(Sys->WL.PotentialFile, 1024, "./Restarts/WLFE_%s_.dat", Sys->SimInfo);
  res += ReadWLFile(Sys, Sys->WL.PotentialFile);
    
  snprintf(Sys->WL.PotentialFile, 1024, "./Restarts/WLFE_.dat");
  res += ReadWLFile(Sys, Sys->WL.PotentialFile);
  
  if(res==0){
    Sys->WL.eps     = 1.0;
    Sys->WL.binnums = (int)arguments->WLbins;
    Sys->WL.zmin    = arguments->WLzmin;
    Sys->WL.zmax    = arguments->WLzmax;

    AllocateWLArrays(Sys);    
    for(int i=0;i<Sys->WL.binnums;i++){
      Sys->WL.FreeEnergy[1][i] = -( Sys->MTot*
				    ReturnParticlePotentialLatticeSite(
				      Sys, Sys->WL.FreeEnergy[0][i], 1.0, 1.0
				    )
				  );  
    }
  }
  snprintf( Sys->WL.PotentialFile, 1024, "%sDataFiles/WLFE_%s_.dat", 
  	    arguments->OutputDir, Sys->SimInfo
    	  );
}

int ErrorCheckWLBias(sys *Sys){
  int Err=0;
  double ztop=0.0, zbot=0.0;
  for(int i = 0;i<Sys->Nx;i++){
    for(int j = 0;j<Sys->Ny;j++){
      ztop += Sys->TopMesh[i][j].Pos.z/Sys->MTot;
      zbot += Sys->BotMesh[i][j].Pos.z/Sys->MTot;
    }
  }
//   ztop-=Sys->MeshBias.theta0;zbot+=Sys->MeshBias.theta0;

  if(fabs(ztop-Sys->MeshBias.thetaTop)>1e-6){
    printf("Position of top is %lf stored %lf\n", ztop, Sys->MeshBias.thetaTop);
    Err=6; 
  }
  if(fabs(zbot-Sys->MeshBias.thetaBot)>1e-6){
    printf("Position of bot is %lf stored %lf\n", zbot, Sys->MeshBias.thetaBot);
    Err=6; 
  }
  double En = WLEnergy(ztop - zbot, Sys) + BiasPot(ztop, zbot, Sys->MeshBias.theta0, Sys->MeshBias.BiasStr);
  if(fabs(En-Sys->BiasEnergy)>1e-6){
    printf("Bias energy off computed %lf stored %lf\n", En, Sys->BiasEnergy);
    Err=6;
  }
  return Err;
}

double PrintBiasPotentialWL(sys *Sys){
  double ztop=0.0, zbot=0.0;
  for(int i = 0;i<Sys->Nx;i++){
    for(int j = 0;j<Sys->Ny;j++){
      ztop += Sys->TopMesh[i][j].Pos.z/Sys->MTot;
      zbot += Sys->BotMesh[i][j].Pos.z/Sys->MTot;
    }
  }
//   ztop-=Sys->MeshBias.theta0;zbot+=Sys->MeshBias.theta0;
  return BiasPot(ztop, zbot, Sys->MeshBias.theta0, Sys->MeshBias.BiasStr) + WLEnergy(ztop - zbot, Sys);
}

void ComputeBiasPotentialWL(sys *Sys){
 double ztop=0.0, zbot=0.0;
  for(int i = 0;i<Sys->Nx;i++){
    for(int j = 0;j<Sys->Ny;j++){
      ztop += Sys->TopMesh[i][j].Pos.z/Sys->MTot;
      zbot += Sys->BotMesh[i][j].Pos.z/Sys->MTot;
    }
  }
//   ztop -= Sys->MeshBias.theta0;zbot += Sys->MeshBias.theta0;

  Sys->MeshBias.thetaTop = ztop;Sys->MeshBias.thetaBot = zbot;
  Sys->BiasEnergy = BiasPot(ztop, zbot, Sys->MeshBias.theta0, Sys->MeshBias.BiasStr) + WLEnergy(ztop - zbot, Sys);
}

double MeshBiasChangeWL(sys *Sys){
  double dzSuggested;
  if(Sys->choice == TOP){
    dzSuggested = (Sys->suggestedz - Sys->TopMesh[Sys->i0][Sys->j0].Pos.z)/Sys->MTot;
    Sys->MeshBias.thetadummy = Sys->MeshBias.thetaTop + dzSuggested;
    if( (Sys->MeshBias.thetadummy - Sys->MeshBias.thetaBot > Sys->WL.zmaxcut) || 
        (Sys->MeshBias.thetadummy - Sys->MeshBias.thetaBot < Sys->WL.zmincut)    ) return ALOT;
    
    return   BiasPot(Sys->MeshBias.thetadummy, Sys->MeshBias.thetaBot, Sys->MeshBias.theta0, Sys->MeshBias.BiasStr)
           - BiasPot(Sys->MeshBias.thetaTop  , Sys->MeshBias.thetaBot, Sys->MeshBias.theta0, Sys->MeshBias.BiasStr)
           + WLEnergy(Sys->MeshBias.thetadummy - Sys->MeshBias.thetaBot, Sys)
           - WLEnergy(Sys->MeshBias.thetaTop   - Sys->MeshBias.thetaBot, Sys);
  }
  else{
    dzSuggested = (Sys->suggestedz - Sys->BotMesh[Sys->i0][Sys->j0].Pos.z)/Sys->MTot;
    Sys->MeshBias.thetadummy = Sys->MeshBias.thetaBot + dzSuggested;
    if( (Sys->MeshBias.thetaTop - Sys->MeshBias.thetadummy > Sys->WL.zmaxcut) || 
        (Sys->MeshBias.thetaTop - Sys->MeshBias.thetadummy < Sys->WL.zmincut)    ) return ALOT;
    
    return   BiasPot(Sys->MeshBias.thetaTop, Sys->MeshBias.thetadummy, Sys->MeshBias.theta0, Sys->MeshBias.BiasStr)
           - BiasPot(Sys->MeshBias.thetaTop, Sys->MeshBias.thetaBot  , Sys->MeshBias.theta0, Sys->MeshBias.BiasStr)
           + WLEnergy(Sys->MeshBias.thetaTop - Sys->MeshBias.thetadummy, Sys)
           - WLEnergy(Sys->MeshBias.thetaTop - Sys->MeshBias.thetaBot  , Sys);
  }
}

double MeshBiasGlobalWL(sys *Sys){
  if(Sys->choice == TOP){
    Sys->MeshBias.thetadummy = Sys->MeshBias.thetaTop + Sys->suggestedz;
    if( (Sys->MeshBias.thetadummy - Sys->MeshBias.thetaBot > Sys->WL.zmaxcut) || 
        (Sys->MeshBias.thetadummy - Sys->MeshBias.thetaBot < Sys->WL.zmincut)    ) return ALOT;
    
    return   BiasPot(Sys->MeshBias.thetadummy, Sys->MeshBias.thetaBot, Sys->MeshBias.theta0, Sys->MeshBias.BiasStr)
           - BiasPot(Sys->MeshBias.thetaTop  , Sys->MeshBias.thetaBot, Sys->MeshBias.theta0, Sys->MeshBias.BiasStr)
           + WLEnergy(Sys->MeshBias.thetadummy - Sys->MeshBias.thetaBot, Sys)
           - WLEnergy(Sys->MeshBias.thetaTop   - Sys->MeshBias.thetaBot, Sys);
  }
  else{
    Sys->MeshBias.thetadummy = Sys->MeshBias.thetaBot + Sys->suggestedz;
    if( (Sys->MeshBias.thetaTop - Sys->MeshBias.thetadummy > Sys->WL.zmaxcut) || 
        (Sys->MeshBias.thetaTop - Sys->MeshBias.thetadummy < Sys->WL.zmincut)    ) return ALOT;
    
    return   BiasPot(Sys->MeshBias.thetaTop, Sys->MeshBias.thetadummy, Sys->MeshBias.theta0, Sys->MeshBias.BiasStr)
           - BiasPot(Sys->MeshBias.thetaTop, Sys->MeshBias.thetaBot  , Sys->MeshBias.theta0, Sys->MeshBias.BiasStr)
           + WLEnergy(Sys->MeshBias.thetaTop - Sys->MeshBias.thetadummy, Sys)
           - WLEnergy(Sys->MeshBias.thetaTop - Sys->MeshBias.thetaBot  , Sys);
  }
}

void StoreBiasChangeWL(sys *Sys, double DeltaEBias){
  if (Sys->choice == TOP) Sys->MeshBias.thetaTop  =  Sys->MeshBias.thetadummy;
  else                    Sys->MeshBias.thetaBot  =  Sys->MeshBias.thetadummy;
  Sys->BiasEnergy +=  DeltaEBias;
}

void UpdateWL(sys *Sys){
  double z = Sys->MeshBias.thetaTop - Sys->MeshBias.thetaBot;
  int imin = Return_WL_Index(z, Sys), i;
  double mean = 0.0, msqr = 0.0, min = Sys->WL.Hitsz[1][0];
    
  Sys->WL.Hitsz[1][imin]      += 1.0;
  Sys->WL.FreeEnergy[1][imin] += Sys->WL.eps;
  Sys->BiasEnergy             += Sys->WL.eps;  
  Sys->TotEnergy              += Sys->WL.eps;
  
  for(i=0;i<Sys->WL.binnums;i++){
    if(Sys->WL.Hitsz[1][i] < min) min = Sys->WL.Hitsz[1][i];
    msqr += Sys->WL.Hitsz[1][i]*Sys->WL.Hitsz[1][i]/Sys->WL.binnums;
    mean += Sys->WL.Hitsz[1][i]/Sys->WL.binnums;
  }
  msqr -= mean*mean;
  
  if( min > 0.8*mean ){
    for(i=0;i<Sys->WL.binnums;i++) Sys->WL.Hitsz[1][i] = 0.0;
    Sys->WL.eps *= 0.5;
  }

  int r = max(1, max(Sys->SweepNum/100, max(1, Sys->SweepNum/500)));
    
  if( (Sys->SweepIndex % r == 0) ){
    double meanF = 0.0, Q = 0.0;
    for(i=0;i<Sys->WL.binnums;i++){
      Q     += Sys->WL.Hitsz[1][i];
      meanF += Sys->WL.FreeEnergy[1][i]/Sys->WL.binnums;
    }
    Q *= Sys->WL.dz;Q=1./Q;
    for(i=0;i<Sys->WL.binnums;i++) Sys->WL.FreeEnergy[1][i] -= meanF;
    Sys->BiasEnergy -= meanF;Sys->TotEnergy -= meanF;

    FILE *fp = fopen(Sys->WL.PotentialFile, "w");
    fprintf(fp, "#z(1) FreeEnergy(2) Probs(3) current eps: %e min: %e mean: %e Norm: %e\n", 
            Sys->WL.eps, min, mean, Q
           );
    for(i=0;i<Sys->WL.binnums;i++){
      fprintf(fp, "%lf %lf %lf\n", 
	      Sys->WL.Hitsz[0][i], Sys->WL.FreeEnergy[1][i], Q*Sys->WL.Hitsz[1][i]
             );
    }
    fflush(fp);fclose(fp);
  }
}

//HarmonicPlusVeff Bias

int ErrorCheckHarmonicPlusVeffBias(sys *Sys){
  int Err=0;
  
  double En = PrintBiasPotentialHarmonicPlusVeff(Sys);
  
  double ztop=0.0, zbot=0.0;
  for(int i = 0;i<Sys->Nx;i++){
    for(int j = 0;j<Sys->Ny;j++){
      ztop += Sys->TopMesh[i][j].Pos.z/Sys->MTot;
      zbot += Sys->BotMesh[i][j].Pos.z/Sys->MTot;
    }
  }
  
  ztop-=Sys->MeshBias.theta0;
  zbot+=Sys->MeshBias.theta0;
  
  if(fabs(ztop-Sys->MeshBias.thetaTop)>1e-6){
    printf("displacement of top is %lf stored %lf\ntheta0: %lf\n", Sys->MeshBias.thetaTop + Sys->MeshBias.theta0, Sys->MeshBias.thetaTop, Sys->MeshBias.theta0);
    Err=6; 
  }
  if(fabs(zbot-Sys->MeshBias.thetaBot)>1e-6){
    printf("displacement of bot is %lf stored %lf\ntheta0: %lf\n", Sys->MeshBias.thetaBot - Sys->MeshBias.theta0, Sys->MeshBias.thetaBot,-Sys->MeshBias.theta0);
    Err=6; 
  }
  
  if(fabs(En-Sys->BiasEnergy)>1e-6){
    printf("Bias energy off computed %lf stored %lf\n", En, Sys->BiasEnergy);
    Err=6;
  }
  
  if(Err!=0){
    printf("\nMeshBiasInfo\n");
    printf("  thetaBot: %lf\n  thetaTop: %lf\n  thetaVal: %lf\n  theta0: %lf\n  VeffVal: %lf\n", 
           Sys->MeshBias.thetaBot, Sys->MeshBias.thetaTop, Sys->MeshBias.thetaTop - Sys->MeshBias.thetaBot + 2.0*Sys->MeshBias.theta0, 2.0*Sys->MeshBias.theta0,
           Sys->MTot*ReturnParticlePotentialLatticeSite(Sys, Sys->MeshBias.thetaTop - Sys->MeshBias.thetaBot + 2.0*Sys->MeshBias.theta0, 1.0, 1.0));
  }
  
  return Err;
}

double PrintBiasPotentialHarmonicPlusVeff(sys *Sys){
  double ztop=0.0, zbot=0.0;
  for(int i = 0;i<Sys->Nx;i++){
    for(int j = 0;j<Sys->Ny;j++){
      ztop += Sys->TopMesh[i][j].Pos.z/Sys->MTot;
      zbot += Sys->BotMesh[i][j].Pos.z/Sys->MTot;
    }
  }
  double thetaAve = ztop - zbot;
  ztop-=Sys->MeshBias.theta0;
  zbot+=Sys->MeshBias.theta0;
  return  Sys->MeshBias.BiasStr*(zbot*zbot + ztop*ztop)
         -Sys->MeshBias.VeffStr*ReturnParticlePotentialLatticeSite(Sys, thetaAve, 1.0, 1.0);
}

void ComputeBiasPotentialHarmonicPlusVeff(sys *Sys){
 double ztop=0.0, zbot=0.0;
  for(int i = 0;i<Sys->Nx;i++){
    for(int j = 0;j<Sys->Ny;j++){
      ztop += Sys->TopMesh[i][j].Pos.z/Sys->MTot;
      zbot += Sys->BotMesh[i][j].Pos.z/Sys->MTot;
    }
  }
  double thetaAve = ztop - zbot;
  ztop-=Sys->MeshBias.theta0;
  zbot+=Sys->MeshBias.theta0;
  
  Sys->MeshBias.thetaTop = ztop;
  Sys->MeshBias.thetaBot = zbot;
  Sys->BiasEnergy        =  Sys->MeshBias.BiasStr*(zbot*zbot + ztop*ztop) 
                          - Sys->MeshBias.VeffStr*ReturnParticlePotentialLatticeSite(Sys, thetaAve, 1.0, 1.0);
}

double MeshBiasChangeHarmonicPlusVeff(sys *Sys){
  double dzSuggested, partialtheta = 2.0*Sys->MeshBias.theta0;
  if(Sys->choice == TOP){
    dzSuggested = (Sys->suggestedz - Sys->TopMesh[Sys->i0][Sys->j0].Pos.z)/Sys->MTot;
    Sys->MeshBias.thetadummy = Sys->MeshBias.thetaTop + dzSuggested;
    partialtheta -= Sys->MeshBias.thetaBot;
    return  Sys->MeshBias.BiasStr*dzSuggested*(dzSuggested + 2.0*Sys->MeshBias.thetaTop)
          - Sys->MeshBias.VeffStr*(
              ReturnParticlePotentialLatticeSite(Sys, Sys->MeshBias.thetadummy + partialtheta, 1.0, 1.0)
             -ReturnParticlePotentialLatticeSite(Sys, Sys->MeshBias.thetaTop   + partialtheta, 1.0, 1.0)
          );
  }
  else{
    dzSuggested = (Sys->suggestedz - Sys->BotMesh[Sys->i0][Sys->j0].Pos.z)/Sys->MTot;
    Sys->MeshBias.thetadummy = Sys->MeshBias.thetaBot + dzSuggested;
    partialtheta += Sys->MeshBias.thetaTop;
    return  Sys->MeshBias.BiasStr*dzSuggested*(dzSuggested + 2.0*Sys->MeshBias.thetaBot)
          - Sys->MeshBias.VeffStr*(
              ReturnParticlePotentialLatticeSite(Sys, partialtheta - Sys->MeshBias.thetadummy, 1.0, 1.0)
             -ReturnParticlePotentialLatticeSite(Sys, partialtheta - Sys->MeshBias.thetaBot  , 1.0, 1.0)
          );
  }
}

double MeshBiasGlobalHarmonicPlusVeff(sys *Sys){
  double partialtheta = 2.0*Sys->MeshBias.theta0;
  
  if(Sys->choice == TOP){
    Sys->MeshBias.thetadummy = Sys->MeshBias.thetaTop + Sys->suggestedz;
    partialtheta -= Sys->MeshBias.thetaBot;
    return   Sys->MeshBias.BiasStr*Sys->suggestedz*(Sys->suggestedz + 2.0*Sys->MeshBias.thetaTop)
           - Sys->MeshBias.VeffStr*(
              ReturnParticlePotentialLatticeSite(Sys, Sys->MeshBias.thetadummy + partialtheta, 1.0, 1.0)
             -ReturnParticlePotentialLatticeSite(Sys, Sys->MeshBias.thetaTop   + partialtheta, 1.0, 1.0)
            );
  }
  else{
    Sys->MeshBias.thetadummy = Sys->MeshBias.thetaBot + Sys->suggestedz;
    partialtheta += Sys->MeshBias.thetaTop;
    return   Sys->MeshBias.BiasStr*Sys->suggestedz*(Sys->suggestedz + 2.0*Sys->MeshBias.thetaBot)
           - Sys->MeshBias.VeffStr*(
              ReturnParticlePotentialLatticeSite(Sys, partialtheta - Sys->MeshBias.thetadummy, 1.0, 1.0)
             -ReturnParticlePotentialLatticeSite(Sys, partialtheta - Sys->MeshBias.thetaBot  , 1.0, 1.0)
            );
  }
}

void StoreBiasChangeHarmonicPlusVeff(sys *Sys, double DeltaEBias){
  if (Sys->choice == BOT) Sys->MeshBias.thetaBot = Sys->MeshBias.thetadummy;
  else             Sys->MeshBias.thetaTop = Sys->MeshBias.thetadummy;
  Sys->BiasEnergy+=DeltaEBias;
}


