#ifndef __pctPhysicalConstants_h
#define __pctPhysicalConstants_h

// Prevent macro redefinition warning caused by SystemOfUnits.h defining pascal
// which is already a keyword in minwindef.h
#ifdef pascal
#  undef pascal
#endif
#include "CLHEP/Units/PhysicalConstants.h"

#endif
