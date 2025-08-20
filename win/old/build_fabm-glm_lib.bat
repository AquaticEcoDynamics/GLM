@rem Script to build the fabm library for glm

@echo off

if "%1%"=="" ( echo no configuration type given
               exit /B 5 )
if "%2%"=="" ( echo no platform type given
               exit /B 5 )

set Configuration=%1%
set Platform=%2%

echo %Platform%-%Configuration%

@rem Check env var VisualStudioVersion
@rem   vs2015 = "14.0"
@rem   vs2017 = "15.0"
@rem
@rem Odd thing to note - batch files include the quotes in variables
@rem so set GV without quotes otherwise they appear in the middle of
@rem the Generator variable.
@rem
@rem However, we need those quotes around the Generator variable
@rem where the cmake call is made
@rem


if "%VisualStudioVersion%"=="14.0" (
  set GV=14 2015
) else if "%VisualStudioVersion%"=="15.0" (
  set GV=15 2017
) else (
  echo Unknown Visual Studio version
)

if "%Platform%"=="x64" (
 set Generator="Visual Studio %GV% Win64"
) else (
 set Generator="Visual Studio %GV%"
)

set startdir=%cd%
set prevdir=%cd%

:loop1
   if EXIST ".\fabm-git\." ( goto :done1 )
   cd ..
   rem echo now in dir %cd%
   rem echo prev in dir %prevdir%
   if "%prevdir%"=="%cd%" (
      rem chdir did nothing - probably at the top of the tree
      echo Cannot find source tree
      chdir %startdir%
      exit /B 6
   )
   set prevdir=%cd%
   goto :loop1

:done1

chdir fabm-git
set FABM_BASE=%cd%
echo Base directories: %FABM_BASE%

echo Current Directory %startdir%
chdir %startdir%
if not EXIST ".\%Platform%-%Configuration%\." ( mkdir "%Platform%-%Configuration%" )
set build_dir=%startdir%\%Platform%-%Configuration%\fabm-glm
echo Build directory: %build_dir%
if not EXIST ".\%build_dir%\." ( mkdir "%build_dir%" )
chdir "%build_dir%"

set install_prefix=%startdir%\%Platform%-%Configuration%
echo Install directory: %install_prefix%

set FABM_HOST=glm

@echo on

cmake "%FABM_BASE%\src" ^
      -G %Generator% ^
      -DFABM_EMBED_VERSION=on ^
      -DFABM_HOST=glm ^
      -DCMAKE_Fortran_COMPILER=ifort ^
      -DCMAKE_VS_PLATFORM_NAME=%Platform% ^
      -DCMAKE_CONFIGURATION_TYPES="Debug;Release" ^
      -DCMAKE_INSTALL_PREFIX="%install_prefix%"

@rem pause

@rem cmake --build . --clean-first --config "%Configuration%" --target INSTALL
cmake --build . --clean-first --config %Configuration% --target INSTALL

@rem pause

@chdir %startdir%

