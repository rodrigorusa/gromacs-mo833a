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
 * Green Red Orange Magenta Azure Cyan Skyblue
 */

#ifndef _rbin_h
#define _rbin_h

static char *SRCID_rbin_h = "$Id$";

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#ifdef HAVE_IDENT
#ident	"@(#) rbin.h 1.8 2/2/97"
#endif /* HAVE_IDENT */
#include "sysstuff.h"
#include "typedefs.h"
#include "network.h"
	
typedef struct {
  int    nreal;
  int    maxreal;
  double *rbuf;
} t_bin;

extern t_bin *mk_bin(void);
/* Create a real bin */

extern void reset_bin(t_bin *b);
/* Reset number of entries to zero */

extern int add_binr(FILE *log,t_bin *b,int nr,real r[]);
extern int add_bind(FILE *log,t_bin *b,int nr,double r[]);
/* Add reals to the bin. Returns index */

extern void sum_bin(t_bin *b,t_commrec *cr);
/* Globally sum the reals in the bin */

extern void extract_binr(t_bin *b,int index,int nr,real r[]);
extern void extract_bind(t_bin *b,int index,int nr,double r[]);
/* Extract values from the bin, starting from index (see add_bin) */

#endif	/* _rbin_h */
