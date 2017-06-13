/******************************************************************************
 *                                                                            *
 * dlfcn.c                                                                    *
 *                                                                            *
 * Copyright 2014 Ambinet Systems.                                            *
 *                                                                            *
 * Description:                                                               *
 *                                                                            *
 ******************************************************************************/

#include <windows.h>

#include <dlfcn.h>


static int s_err = 0;

/******************************************************************************/
static const char *err_string(void)
{
    switch (s_err) {
        case ERROR_SUCCESS:
            return NULL;
            break;
        case ERROR_MOD_NOT_FOUND:
            return "Module not found";
            break;
        case ERROR_PROC_NOT_FOUND:
            return "Symbol not found";
            break;
        case ERROR_BAD_EXE_FORMAT:
            return "Invalid image format";
            break;
        case ERROR_SHARING_VIOLATION:
            return "Sharing violation";
            break;
        default:
            return "Operation failed";
            break;
    }
}
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/


/******************************************************************************/
void *dlopen(char const *moduleName, int mode)
{
    HMODULE	hModule;

    if (mode & RTLD_LAZY)
        OutputDebugStringA("Library does not support RTLD_LAZY; using RTLD_NOW\n");

    hModule = LoadLibraryA(moduleName);

    if (hModule == NULL)
        s_err = (GetLastError());
    else
        s_err = (ERROR_SUCCESS);

    return hModule;
}
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/


/******************************************************************************/
int dlclose(void *hModule)
{
    if(!FreeLibrary((HMODULE)hModule)) {
        s_err = GetLastError();

        return 1;
    }

    s_err = (ERROR_SUCCESS);

    return 0;
}
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/


/******************************************************************************/
void *dlsym(void *hModule, char const *symbolName)
{
    void *symbol = GetProcAddress((HMODULE)hModule, symbolName);

    if (symbol == NULL)
        s_err = GetLastError();
    else
        s_err = ERROR_SUCCESS;

    return symbol;
}
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/


/******************************************************************************/
char const *dlerror(void)
{
    char const *err = err_string();

    s_err = 0;

    return err;
}
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
