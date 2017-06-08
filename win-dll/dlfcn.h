/******************************************************************************
 *                                                                            *
 * dlfcn.h                                                                    *
 *                                                                            *
 * Copyright 2014 Ambinet Systems.                                            *
 *                                                                            *
 * Description:                                                               *
 *                                                                            *
 ******************************************************************************/
#ifndef _DL_FUNCTIONS_
#define _DL_FUNCTIONS_

#include <sys/types.h>

#define RTLD_LAZY   0x00001              /* Lazy function call binding.       */
#define RTLD_NOW    0x00002              /* Immediate function call binding.  */

void *dlopen(char const *moduleName, int mode);
int   dlclose(void *hModule);
void *dlsym(void *hModule, char const *symbolName);
char const *dlerror(void);

#endif
