#include "utilities.h"
#include "meshbiasfuncs.h"
#include "mesh.h"

void ExtractValuesFromStringList( void *array, char *opt, char *str, int flavornums){
  int i=0;
  if( strcmp(opt, "INT") == 0 ){
    ((int *)array)[0]= atoi(strtok(str , " "));
    for (i=1;i<flavornums;i++){
      ((int *)array)[i]=atoi(strtok(NULL , " "));
    }
  }
  else if( strcmp(opt, "DOUBLE")== 0 ){
   ((double *)array)[0]= atof(strtok(str , " "));
   for (i=1;i<flavornums;i++){
     ((double *)array)[i]=atof(strtok(NULL , " "));
    }
  }
  else{printf("Could not determine what this OPT: %s is\n", opt);}
}

int IndexPBC(int index, int Nmax){
  if (index<0) return index + Nmax;
  else if (index>=Nmax) return index - Nmax;
  else return index;
}

int indexDiffPBC(int dindex, int Nmax){
  if      (dindex<=-Nmax/2) return dindex + Nmax;
  else if (dindex> Nmax/2) return dindex - Nmax;
  else return dindex;
}

int computeindexdistsqr(int di, int dj, int Nx, int Ny){
  int i = indexDiffPBC(di, Nx), j = indexDiffPBC(dj, Ny);return i*i + j*j;
}

int ReturndrIndex(int di, int dj, int Nx, int Ny, int *rsqrlookup, int lenrvals, int int_rcutoffsqr){
  int int_drsqr = computeindexdistsqr(di, dj, Nx, Ny);
  int t=-1;
  if(int_drsqr > int_rcutoffsqr) return -1;
  for(t=0;t<lenrvals;t++){
    if(int_drsqr == rsqrlookup[t]) break;
  }
  return t;
}

void CopyVec(vec *target, vec source){
 target->x=source.x; 
 target->y=source.y;
 target->z=source.z; 
}

double ReturnHeightStd(sys *Sys){
  double sigmasqr=0.0, kx, ky;
  int i,j;
  for(i=0;i<Sys->Nx;i++){
    kx=2.0*M_PI*(i-Sys->Nx/2)/Sys->Lx;
    for(j=0;j<Sys->Ny;j++){
      ky=2.0*M_PI*(j-Sys->Ny/2)/Sys->Ly;
    if(kx*kx + ky*ky>0.0)sigmasqr += 1.0/(Sys->A*Sys->kc*(kx*kx + ky*ky)*(kx*kx + ky*ky));
    }
  }
  return sqrt(sigmasqr);
}

double ReturnTotalArea(sys *Sys, int choice){
  MeshPoint **MyMesh = (choice == 0 ? Sys->BotMesh: Sys->TopMesh);
  double Area=0.0;
  for(int i=0;i<Sys->Nx;i++){
    for(int j=0;j<Sys->Ny;j++){
      Area+=MyMesh[i][j].nNorm;
    }
  }
  return Area*Sys->dx*Sys->dy;
}

void ReturnMeshHeights(double *zbot, double *ztop, sys *Sys){
  *zbot = *ztop = 0.0;
  for(int i = 0;i<Sys->Nx;i++){
    for(int j = 0;j<Sys->Ny;j++){
      *ztop += Sys->TopMesh[i][j].Pos.z/Sys->MTot;
      *zbot += Sys->BotMesh[i][j].Pos.z/Sys->MTot;
    }
  }
}

// void ComputeMeshHeights(double *zbot, double *ztop, sys *Sys){
//   *zbot = *ztop = 0.0;
//   for(int i = 0;i<Sys->Nx;i++){
//     for(int j = 0;j<Sys->Ny;j++){
//       *ztop += Sys->TopMesh[i][j].Pos.z/Sys->MTot;
//       *zbot += Sys->BotMesh[i][j].Pos.z/Sys->MTot;
//     }
//   }
// }

void PrintVec(vec target, char *note){printf("%s: %lf %lf %lf\n", note, target.x, target.y, target.z);}

void PrintConfiguration(FILE *fp, sys *Sys){
  int i, j;
  
  fprintf(fp, "%ld\nSweep no: %ld out of %ld\n", Sys->DoF, Sys->SweepIndex, Sys->SweepNum);
  for (i=0;i<Sys->Nx;i++){
   for (j=0;j<Sys->Ny;j++){
    fprintf(
      fp, "N %lf %lf %lf %lf %lf %lf %lf %lf\n", 
      Sys->TopMesh[i][j].Pos.x,  Sys->TopMesh[i][j].Pos.y,  Sys->TopMesh[i][j].Pos.z,
      Sys->TopMesh[i][j].nhat.x, Sys->TopMesh[i][j].nhat.y, Sys->TopMesh[i][j].nhat.z,
      Sys->TopMesh[i][j].lc_x, Sys->TopMesh[i][j].lc_y
    );
   }
  }
  
  for (i=0;i<Sys->Nx;i++){
   for (j=0;j<Sys->Ny;j++){
    fprintf(
      fp, "N %lf %lf %lf %lf %lf %lf %lf %lf\n", 
      Sys->BotMesh[i][j].Pos.x,  Sys->BotMesh[i][j].Pos.y,  Sys->BotMesh[i][j].Pos.z,
      Sys->BotMesh[i][j].nhat.x, Sys->BotMesh[i][j].nhat.y, Sys->BotMesh[i][j].nhat.z,
      Sys->BotMesh[i][j].lc_x, Sys->BotMesh[i][j].lc_y
    );
   }
  }
}

void PrintData(FILE *fp, sys *Sys, long currentMCstep){
  int i, shift = 17;
  double zbot, ztop;
  ReturnMeshHeights(&zbot, &ztop, Sys);
  
  if (currentMCstep==0){
    fprintf(fp, "#Sweep(1) TAcc(2) NAcc(3) SHAcc(4) TotE(5) MeshE(6) BiasE(7) Veff(8) ZBot(9) ZTop(10) TotAreaBot(11) TotAreaTop(12) TotPop(13) TotPopBot(14) TotPopTop(15) ThetaValTop(16) ThetaValBot(17) ");
    for(i = 0; i<Sys->NumPTypes; i++) fprintf(fp, "NBot%d(%d) ",     i+1, shift+1+i);
    for(i = 0; i<Sys->NumPTypes; i++) fprintf(fp, "NTop%d(%d) ",     i+1, shift + 1 +   Sys->NumPTypes + i);
    for(i = 0; i<Sys->NumPTypes; i++) fprintf(fp, "Bonds%d(%d) ",    i+1, shift + 1 + 2*Sys->NumPTypes + i);
    fprintf(fp, "\n");
  }
    
  fprintf( fp, "%ld %.2lf %.2lf %.2lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf ", 
    Sys->SweepIndex, 
    Sys->AcceptanceRate/(currentMCstep+1), 
    Sys->NodeAccept/(currentMCstep + 1 - Sys->ShiftAttempt),
    Sys->ShiftAccept/(Sys->ShiftAttempt+1),
    Sys->TotEnergy, Sys->MeshEnergy, Sys->BiasEnergy, Sys->TotVeff, zbot, ztop,
    ReturnTotalArea(Sys, 0), ReturnTotalArea(Sys, 1),
    Sys->TotCurPop[0] + Sys->TotCurPop[1], Sys->TotCurPop[0], Sys->TotCurPop[1], Sys->MeshBias.thetaTop, Sys->MeshBias.thetaBot
  );
  
  for(i = 0; i<Sys->NumPTypes; i++) fprintf(fp, "%lf ", Sys->CurParPop[0][i]);
  for(i = 0; i<Sys->NumPTypes; i++) fprintf(fp, "%lf ", Sys->CurParPop[1][i]);
  for(i = 0; i<Sys->NumPTypes; i++) fprintf(fp, "%lf ", Sys->BondNums[i]);
  fprintf(fp,"\n");
}

void PrintDataToScreen(FILE *fp, sys *Sys, long currentMCstep){
  double zbot, ztop;
  ReturnMeshHeights(&zbot, &ztop, Sys);
  
  if (currentMCstep==0) fprintf(fp, "Sweep(1) TAcc(2) NAcc(3) SHAcc(4) TotE(5) MeshE(6) BiasE(7) Veff(8) MemGap(9) TotPopBot(10) TotPopTop(11) ThetaValTop(12) ThetaValBot(13)\n");
  fprintf( fp, "%ld %.2lf %.2lf %.2lf %lf %lf %lf %lf %lf %.0lf %.0lf %lf %lf\n", 
    Sys->SweepIndex, 
    Sys->AcceptanceRate/(currentMCstep+1),
    Sys->NodeAccept/(currentMCstep + 1 - Sys->ShiftAttempt),
    Sys->ShiftAccept/(Sys->ShiftAttempt+1),
    Sys->TotEnergy, Sys->MeshEnergy, Sys->BiasEnergy, Sys->TotVeff, ztop - zbot, 
    Sys->TotCurPop[0], Sys->TotCurPop[1], Sys->MeshBias.thetaTop, Sys->MeshBias.thetaBot
  );
}


// void UpdateCorrelations(sys *Sys){
//   if(Sys->SweepIndex>0 && Sys->Hhist.binnums > 1){
//     ComputeCurrentCorrelations(Sys);
//     corrfuncs *CorrFuncs = &(Sys->CorrFuncs);
//     int samplenums = Sys->SweepIndex/Sys->datarate, rnums=CorrFuncs->lenrvals, statenums = CorrFuncs->statenums, compnums = CorrFuncs->compnums;
//     int r, p1, p2, p1prime, p2prime, tmp;
//     
//     for(r=0;r<rnums;r++){for(p1=0;p1<statenums;p1++){for(p2=0;p2<statenums;p2++)CorrFuncs->Ave_Correlations[r][p1][p2] += CorrFuncs->Cur_Correlations[r][p1][p2];}}
//     
//     for(p1=0;p1<statenums;p1++){CorrFuncs->Ave_MeanProbs[p1] += CorrFuncs->Cur_MeanProbs[p1];}
//     
// //     FILE *fp=fopen(Sys->Corrfile, "w");
// //     fprintf(fp, "#Samplenums: %d columns are rvals, then <P(i_bot, j_top) P(i'_bot, j'_top)>, note that P(%d, %d) refers to empty-empty. Note that entry %d onwards, we also have mean P(i_bot, j_top) \n", samplenums, compnums-1, compnums-1, statenums*statenums);
// //     fprintf(fp, "# rvals(1) ");
//     
//     tmp=2;
//     for(p1=0;p1<compnums;p1++){
//       for(p2=0;p2<compnums;p2++){
//         for(p1prime=0;p1prime<compnums;p1prime++){
//           for(p2prime=0;p2prime<compnums;p2prime++){
//             fprintf(fp, "<P(%d,%d)P(%d,%d)>(%d) ", p1, p2, p1prime, p2prime, tmp);
//             tmp++;
//           }
//         }
//       }
//     }
//     for(p1=0;p1<compnums;p1++){for(p2=0;p2<compnums;p2++){fprintf(fp, "Pbar(%d,%d)(%d) ", p1, p2, tmp);tmp++;}}
//     fprintf(fp, "\n");
//     
//     for(r=0;r<rnums;r++){
//       fprintf(fp, "%lf ", CorrFuncs->rvals[r]);
//       for(p1=0;p1<statenums;p1++){
//         for(p2=0;p2<statenums;p2++){
//           fprintf(fp, "%lf ", CorrFuncs->Ave_Correlations[r][p1][p2]/samplenums);
//         }
//       }
//       for(p1=0;p1<statenums;p1++) fprintf(fp, "%lf ", CorrFuncs->Ave_MeanProbs[p1]/samplenums);
//       fprintf(fp, "\n");
//     }
//     fclose(fp);
//   }
// }


// entry=0;
//             
//             for(p1=0;p1<CompNums;p1++){
//               for(p2=p1;p2<CompNums;p2++){
//                 c1 = c2 = 0.0;
//                 for(prime=0;prime<CompNums;prime++){
//                   c1 += BotMesh[i][j].ParticleProbs[p1][prime]           - AveProbs[p1][prime];
//                   c2 += BotMesh[iprime][jprime].ParticleProbs[p2][prime] - AveProbs[p2][prime];
//                 }
//                 CorrFuncs->Cur_Correlations[rindex][entry] += c1*c2;
//                 entry++;
//               }
//             }
//             for(p1=0;p1<CompNums;p1++){
//               for(p2=p1;p2<CompNums;p2++){
//                 c1 = c2 = 0.0;
//                 for(prime=0;prime<CompNums;prime++){
//                   c1 += TopMesh[i][j].ParticleProbs[prime][p1]           - AveProbs[prime][p1];
//                   c2 += TopMesh[iprime][jprime].ParticleProbs[prime][p2] - AveProbs[prime][p2];
//                 }
//                 CorrFuncs->Cur_Correlations[rindex][entry] += c1*c2;
//                 entry++;
//               }
//             }
//             
//             for(p1=0;p1<CompNums;p1++){
//               CorrFuncs->Cur_Correlations[rindex][entry] += (BotMesh[i][j].ParticleProbs[p1][p1]           - AveProbs[p1][p1])
//                                                            *(TopMesh[iprime][jprime].ParticleProbs[p1][p1] - AveProbs[p1][p1]);
//               entry++;
//             }
//             
//             for(p1=0;p1<CompNums-1;p1++){
//               CorrFuncs->Cur_Correlations[rindex][entry] += (BotMesh[i][j].ParticleProbs[p1][CompNums-1]           - AveProbs[p1][CompNums-1])
//                                                            *(TopMesh[iprime][jprime].ParticleProbs[p1][CompNums-1] - AveProbs[p1][CompNums-1]);
//               entry++;
//             }
//             for(p1=0;p1<CompNums-1;p1++){
//               CorrFuncs->Cur_Correlations[rindex][entry] += (BotMesh[i][j].ParticleProbs[CompNums-1][p1]           - AveProbs[CompNums-1][p1])
//                                                            *(TopMesh[iprime][jprime].ParticleProbs[CompNums-1][p1] - AveProbs[CompNums-1][p1]);
//               entry++;
//             }


