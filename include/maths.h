/*
 * $Id$
 * 
 *       This source code is part of
 * 
 *        G   R   O   M   A   C   S
 * 
 * GROningen MAchine for Chemical Simulations
 * 
 *               VERSION 2.0
 * 
 * Copyright (c) 1991-1999
 * BIOSON Research Institute, Dept. of Biophysical Chemistry
 * University of Groningen, The Netherlands
 * 
 * Please refer to:
 * GROMACS: A message-passing parallel molecular dynamics implementation
 * H.J.C. Berendsen, D. van der Spoel and R. van Drunen
 * Comp. Phys. Comm. 91, 43-56 (1995)
 * 
 * Also check out our WWW page:
 * http://md.chem.rug.nl/~gmx
 * or e-mail to:
 * gromacs@chem.rug.nl
 * 
 * And Hey:
 * Good ROcking Metal Altar for Chronical Sinners
 */

#ifndef _maths_h
#define _maths_h

static char *SRCID_maths_h = "$Id$";

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#ifdef HAVE_IDENT
#ident	"@(#) maths.h 1.11 11/24/92"
#endif /* HAVE_IDENT */

#include <math.h>
#include "typedefs.h"

#ifndef M_PI
#define	M_PI		3.14159265358979323846
#endif

#ifndef M_PI_2
#define	M_PI_2		1.57079632679489661923
#endif

#ifndef M_2PI
#define	M_2PI		6.28318530718
#endif

#ifndef M_SQRT2
#define M_SQRT2 sqrt(2.0)
#endif

#ifdef CPLUSPLUS
extern "C" {
#endif

extern	int		gmx_nint(real a);
extern  real            sign(real x,real y);

#ifdef CPLUSPLUS
}
#endif

#endif	/* _maths_h */
