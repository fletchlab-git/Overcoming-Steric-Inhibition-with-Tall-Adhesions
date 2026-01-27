#ifndef INIT_H
#define INIT_H

#include "parseparameters.h"
#include "membranesys.h"
#include "meshbiasfuncs.h"
#include "errorcheck.h"
#include "utilities.h"
#include "mesh.h"
#include "effectiveparticlepotential.h"



#include <string.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <dirent.h>

void SetOutPutFiles(sys *Sys, struct arguments *arguments );
void InitializeStructs(sys *Sys, struct arguments *arguments);
void ModifyInitialConfig(sys *Sys, double z0);
void WriteRestart(sys *Sys, struct arguments *arguments);
#endif
