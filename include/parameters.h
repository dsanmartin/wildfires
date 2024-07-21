#ifndef PARAMETERS_H
#define PARAMETERS_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h> 
#include "structures.h"
#include "utils.h"

Parameters read_parameters_file(const char *filePath); 

#endif