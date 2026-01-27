#include "parseparameters.h"
#include "membranesys.h"
#include "utilities.h"

// 
const char *__restrict__ ArgumentPrint(){
  return "\nMC properties:\n"
  "   MC Sweeps:       %d\n"
  "   Mesh disp:       %lf\n"
  "   Bias type:       %s\n"
  "   Bias strength:   %lf\n"
  "   Bias location:   %lf\n"
  "   InitialHeight:   %lf\n"
  "   RandSeed:        %ld\n"
  "   Boundary Cond:   %s\n"
  "\nMesh Properties:\n"
  "   Bending Modulus: %lf\n"
  "   Lx:              %lf\n"
  "   Ly:              %lf\n"
  "   Node nums in x:  %d\n"
  "   Node nums in y:  %d\n"
  "\nParticle Properties:\n"
  "   Ptype num:       %d\n"
  "   Lengths:         %s\n"
  "   Binding str:     %s\n"
  "   Binding ranges:  %s\n"
  "   Densities:       %s\n"
  "   PDiam:           %f\n"
  "\nOutput Options:\n"
  "   Tag:             %s\n"
  "   OutputDir:       %s\n"
  "   Framenums:       %d\n"
  "   Datanums:        %d\n"
  "   Histbnums:       %d\n"
  "\nRestart Info:\n"
  "   Inputfilename:   %s\n"
  "   Equil option:    %s\n"
  "\nError Check:        %s\n"
  "\n";
}

const char *__restrict__ WritePopBias(){
  return "Bond Bias Params:\n"
  "   Bond Strengths:  %s\n"
  "   Bond Locations:  %s\n\n";
}

//Replace all occurances of %s above with [^n]
const char *__restrict__ ArgumentRead(){
  return "\nMC properties:\n"
  "   MC Sweeps:       %d\n"
  "   Mesh disp:       %lf\n"
  "   Bias type    :   %[^\n]\n"
  "   Bias strength:   %lf\n"
  "   Bias location:   %lf\n"
  "   InitialHeight:   %lf\n"
  "   RandSeed:        %ld\n"
  "   Boundary Cond:   %[^\n]\n"
  "\nMesh Properties:\n"
  "   Bending Modulus: %lf\n"
  "   Lx:              %lf\n"
  "   Ly:              %lf\n"
  "   Node nums in x:  %d\n"
  "   Node nums in y:  %d\n"
  "\nParticle Properties:\n"
  "   Ptype num:       %d\n"
  "   Lengths:         %[^\n]"
  "   Binding str:     %[^\n]"
  "   Binding ranges:  %[^\n]"
  "   Densities:       %[^\n]"
  "   PDiam:           %lf\n"
  "\nOutput Options:\n"
  "   Tag:             %[^\n]"
  "   OutputDir:       %[^\n]"
  "   Framenums:       %d\n"
  "   Datanums:        %d\n"
  "   Histbnums:       %d\n"
  "\nRestart Info:\n"
  "   Inputfilename:   %[^\n]"
  "   Equil option:    %[^\n]"
  "\nError Check:        %[^\n]\n"
  "\n";
}

const char *__restrict__ ReadPopBias(){
  return "Bond Bias Params:\n"
  "   Bond Strengths:  %[^n]\n"
  "   Bond Locations:  %[^n]\n\n";
}

/* Our argp parser. */
static error_t parse_opt (int key, char *arg, struct argp_state *state){
/* Parse each option, one at a time. */
/* Get the input argument from argp_parse, which we
     know is a pointer to our arguments structure. */
  struct arguments *arguments = state->input;
  
  switch (key)
    {
    
    case 'a':sscanf(arg, "%lf", &arguments->Lx);break;
    case 'b':sscanf(arg, "%lf", &arguments->Ly);break;
    case 'c':sscanf(arg, "%d",  &arguments->Nx);break;
    case 'd':sscanf(arg, "%d",  &arguments->Ny);break;
    case 'e':sscanf(arg, "%lf", &arguments->z0);break;
    case 'f':sscanf(arg, "%lf", &arguments->MDisp);break;
    case 'g':sscanf(arg, "%lf", &arguments->kc);break;
    
    case 'h':sscanf(arg, "%d",  &arguments->NumPType);break;
    case 'i':strcpy(arguments->LPType,    arg);break;
    case 'j':strcpy(arguments->BindPType, arg);break;
    case 'k':strcpy(arguments->RanPType,  arg);break;
    
    case 'l':strcpy(arguments->Densities,  arg);break;
        
    case 'm':strcpy(arguments->OutputDir, arg);break;
    case 'n':strcpy(arguments->Tag,  arg);break;
    
    case 'o':sscanf(arg, "%d",  &arguments->framenums);break;
    case 'p':sscanf(arg, "%d",  &arguments->datnums);break;
    case 'q':sscanf(arg, "%d",  &arguments->Histbnums);break;
        
    case 'r':strcpy(arguments->Inputfile, arg);break;
    case 's':strcpy(arguments->Equil,     arg);break;
    
    case 't':sscanf(arg, "%d",  &arguments->MCSweeps);break;
    
    case 'u':sscanf(arg, "%lf", &arguments->BiasStr);break;
    case 'v':strcpy(arguments->BiasType, arg);break;
    case 'w':sscanf(arg, "%lf", &arguments->theta0);break;
    case 'x':sscanf(arg, "%d",  &arguments->WLbins);break;
    case 'y':sscanf(arg, "%lf", &arguments->WLzmin);break;
    case 'z':sscanf(arg, "%lf", &arguments->WLzmax);break;
     
    case 'A':strcpy(arguments->PopStr,  arg);break;
    case 'B':strcpy(arguments->PopLoc,  arg);break;
    
    case 'C':strcpy(arguments->ErrorCheck,  arg);break;
    case 'D':sscanf(arg, "%ld", &arguments->RandSeed);break;
    case 'E':strcpy(arguments->BC,  arg);break;
    case 'F':sscanf(arg, "%lf", &arguments->pDiam);break;
    
    
    default:
      return ARGP_ERR_UNKNOWN;
    }
  return 0;
}

static struct argp argp = { options, parse_opt, "", "" };

void InterpretInputCommands( struct arguments *arguments, int argc, char **argv){
  
//Default values
  
//Mesh properties
  arguments->Lx=120.0;
  arguments->Ly=120.0;
  arguments->Nx=40;
  arguments->Ny=40;
  arguments->kc=20.0;
  
//Particle properties
  
//   arguments->NumPType=0;

//   arguments->NumPType=1;
//   snprintf(arguments->LPType,    1024, "0.0 0.0");
//   snprintf(arguments->BindPType, 1024, "4.0");
//   snprintf(arguments->RanPType,  1024, "4.0");
//   snprintf(arguments->Densities, 1024, "1000.0 1000.0");

//   arguments->NumPType=2;
//   snprintf(arguments->LPType,    1024, "0.0 2.0 0.0 2.0");
//   snprintf(arguments->BindPType, 1024, "4.0 0.0");
//   snprintf(arguments->RanPType,  1024, "4.0 4.0");
//   snprintf(arguments->Densities,  1024, "800.0 10.0 200.0 100.0");
  
  arguments->NumPType = 3;
  snprintf(arguments->LPType,    1024, "5.5 17.5 40.0 5.5 17.5 40.0");
  snprintf(arguments->BindPType, 1024, "10.5 0.0 4.0");
  snprintf(arguments->RanPType,  1024, "4.0 0.0 4.0");
  snprintf(arguments->Densities, 1024, "200 0 8000 200 0 8000");
  
  arguments->pDiam = 4.0;
  
//MC properties
  arguments->MDisp    = 0.4;
  arguments->z0       = 60.0;
  arguments->BiasStr  = 1.0;
  arguments->theta0   = 82.0;
  snprintf(arguments->BiasType, 1024, "Minimum");
  
//   snprintf(arguments->BiasType, 1024, "Pop");
  snprintf(arguments->PopStr,   1024, "6.0 6.0");
  snprintf(arguments->PopLoc,   1024, "0.5 0.5");
  
  snprintf(arguments->BC,   1024, "Frame");
  
  arguments->WLzmin   = 85.0;
  arguments->WLzmax   = 135.0;
  arguments->WLbins   = 501;
  arguments->MCSweeps = 100000;
  arguments->RandSeed = 8431446;
  
// Output options  
  arguments->framenums = 100;
  arguments->datnums   = 1000;
  arguments->Histbnums = 0;
  sscanf("/home/jaffar/Desktop/werk/membranesims/Data/", "%s", arguments->OutputDir);
//   arguments->Tag[0]='\0';
  strncpy(arguments->Tag, "Test", 1024);
  
// Input options
  strncpy(arguments->Inputfile, "None", 1024);
  strncpy(arguments->Equil,  "No", 1024);
  
// ErrorChecking argument
  strncpy(arguments->ErrorCheck, "Naw", 1024);
    
//Parse Arguments
  argp_parse (&argp, argc, argv, 0, 0, arguments);  
}

void PrintInputCommands(FILE *fp, struct arguments *arguments){
  fprintf(fp, ArgumentPrint(),
              arguments->MCSweeps,
              arguments->MDisp,
              arguments->BiasType, arguments->BiasStr, arguments->theta0, arguments->z0,
              arguments->RandSeed, arguments->BC,
              arguments->kc, arguments->Lx, arguments->Ly,
              arguments->Nx, arguments->Ny,
              arguments->NumPType,
              arguments->LPType, arguments->BindPType, arguments->RanPType,
              arguments->Densities, arguments->pDiam,
              arguments->Tag, arguments->OutputDir,
              arguments->framenums, arguments->datnums, arguments->Histbnums,
              arguments->Inputfile, arguments->Equil,
              arguments->ErrorCheck
         );
  if(strcmp(arguments->BiasType, "Pop")==0) fprintf(fp, WritePopBias(), arguments->PopStr, arguments->PopLoc);
}

void ConfigurationRestart(FILE *fp, sys *Sys){
  int i, j;
  for (i=0;i<Sys->Nx;i++){
   for (j=0;j<Sys->Ny;j++){
    fprintf(
      fp, "N %lf %lf %lf %lf %lf %lf\n", 
      Sys->TopMesh[i][j].Pos.x,  Sys->TopMesh[i][j].Pos.y,  Sys->TopMesh[i][j].Pos.z,
      Sys->TopMesh[i][j].nhat.x, Sys->TopMesh[i][j].nhat.y, Sys->TopMesh[i][j].nhat.z
    );
   }
  }
  for (i=0;i<Sys->Nx;i++){
   for (j=0;j<Sys->Ny;j++){
    fprintf(
      fp, "N %lf %lf %lf %lf %lf %lf\n", 
      Sys->BotMesh[i][j].Pos.x,  Sys->BotMesh[i][j].Pos.y,  Sys->BotMesh[i][j].Pos.z,
      Sys->BotMesh[i][j].nhat.x, Sys->BotMesh[i][j].nhat.y, Sys->BotMesh[i][j].nhat.z
    );
   }
  }
}

void WriteRestart(sys *Sys, struct arguments *arguments){
  Sys->fileptr_restart = fopen(Sys->Restartfile, "w");
  PrintInputCommands(Sys->fileptr_restart, arguments);
  ConfigurationRestart(Sys->fileptr_restart, Sys);
  fclose(Sys->fileptr_restart);
}

void ReadConfiguration(FILE *fp, sys *Sys){
  char line[1024], *dummy;
  int i,j;
  dummy = fgets(line, 1024, fp);
  for (i=0;i<Sys->Nx;i++){
   for (j=0;j<Sys->Ny;j++){
    sscanf(
      line, "N %lf %lf %lf %lf %lf %lf\n", 
      &(Sys->TopMesh[i][j].Pos.x),  &(Sys->TopMesh[i][j].Pos.y),  &(Sys->TopMesh[i][j].Pos.z),
      &(Sys->TopMesh[i][j].nhat.x), &(Sys->TopMesh[i][j].nhat.y), &(Sys->TopMesh[i][j].nhat.z)
    );
    dummy = fgets(line, 1024, fp);
   }
  }
  
  for (i=0;i<Sys->Nx;i++){
   for (j=0;j<Sys->Ny;j++){
    sscanf(
      line, "N %lf %lf %lf %lf %lf %lf\n", 
      &(Sys->BotMesh[i][j].Pos.x),  &(Sys->BotMesh[i][j].Pos.y),  &(Sys->BotMesh[i][j].Pos.z),
      &(Sys->BotMesh[i][j].nhat.x), &(Sys->BotMesh[i][j].nhat.y), &(Sys->BotMesh[i][j].nhat.z)
    );
    dummy = fgets(line, 1024, fp);
   }
  }
  dummy=dummy;
  fclose(fp);
}

void ReadInputCommands(FILE *fp, struct arguments *arguments){
  int tmp=fscanf(fp, ArgumentRead(), 
                     &arguments->MCSweeps,
                     &arguments->MDisp,
                     &arguments->BiasType, &arguments->BiasStr, &arguments->theta0,
                     &arguments->z0,
                     &arguments->RandSeed, &arguments->BC,
                     &arguments->kc, &arguments->Lx, &arguments->Ly,
                     &arguments->Nx, &arguments->Ny,
                     &arguments->NumPType,
                     &arguments->LPType, &arguments->BindPType, &arguments->RanPType,
                     &arguments->Densities, &arguments->pDiam,
                     &arguments->Tag, &arguments->OutputDir,
                     &arguments->framenums, &arguments->datnums, &arguments->Histbnums,
                     &arguments->Inputfile, &arguments->Equil,
                     &arguments->ErrorCheck
                );
  if(strcmp(arguments->BiasType, "Pop")==0) tmp += fscanf(fp, ReadPopBias(), &arguments->PopStr, &arguments->PopLoc);
  tmp=tmp;
}
