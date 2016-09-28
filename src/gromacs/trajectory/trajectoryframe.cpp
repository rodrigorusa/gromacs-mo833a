/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2016, by the GROMACS development team, led by
 * Mark Abraham, David van der Spoel, Berk Hess, and Erik Lindahl,
 * and including many others, as listed in the AUTHORS file in the
 * top-level source directory and at http://www.gromacs.org.
 *
 * GROMACS is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public License
 * as published by the Free Software Foundation; either version 2.1
 * of the License, or (at your option) any later version.
 *
 * GROMACS is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with GROMACS; if not, see
 * http://www.gnu.org/licenses, or write to the Free Software Foundation,
 * Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA.
 *
 * If you want to redistribute modifications to GROMACS, please
 * consider that scientific software is very special. Version
 * control is crucial - bugs must be traceable. We will be happy to
 * consider code for inclusion in the official distribution, but
 * derived work must not be called official GROMACS. Details are found
 * in the README & COPYING files - if they are missing, get the
 * official version at http://www.gromacs.org.
 *
 * To help us fund GROMACS development, we humbly ask that you cite
 * the research papers on the package. Check out http://www.gromacs.org.
 */
#include "gmxpre.h"

#include "trajectoryframe.h"

#include <cstdio>

#include <algorithm>

#include "gromacs/math/veccompare.h"
#include "gromacs/topology/atoms.h"
#include "gromacs/utility/compare.h"

void comp_frame(FILE *fp, t_trxframe *fr1, t_trxframe *fr2,
                gmx_bool bRMSD, real ftol, real abstol)
{
    fprintf(fp, "\n");
    cmp_int(fp, "not_ok", -1, fr1->not_ok, fr2->not_ok);
    cmp_int(fp, "natoms", -1, fr1->natoms, fr2->natoms);
    if (cmp_bool(fp, "bTitle", -1, fr1->bTitle, fr2->bTitle))
    {
        cmp_str(fp, "title", -1, fr1->title, fr2->title);
    }
    if (cmp_bool(fp, "bStep", -1, fr1->bStep, fr2->bStep))
    {
        cmp_int(fp, "step", -1, fr1->step, fr2->step);
    }
    cmp_int(fp, "step", -1, fr1->step, fr2->step);
    if (cmp_bool(fp, "bTime", -1, fr1->bTime, fr2->bTime))
    {
        cmp_real(fp, "time", -1, fr1->time, fr2->time, ftol, abstol);
    }
    if (cmp_bool(fp, "bLambda", -1, fr1->bLambda, fr2->bLambda))
    {
        cmp_real(fp, "lambda", -1, fr1->lambda, fr2->lambda, ftol, abstol);
    }
    if (cmp_bool(fp, "bAtoms", -1, fr1->bAtoms, fr2->bAtoms))
    {
        cmp_atoms(fp, fr1->atoms, fr2->atoms, ftol, abstol);
    }
    if (cmp_bool(fp, "bPrec", -1, fr1->bPrec, fr2->bPrec))
    {
        cmp_real(fp, "prec", -1, fr1->prec, fr2->prec, ftol, abstol);
    }
    if (cmp_bool(fp, "bX", -1, fr1->bX, fr2->bX))
    {
        cmp_rvecs(fp, "x", std::min(fr1->natoms, fr2->natoms), fr1->x, fr2->x, bRMSD, ftol, abstol);
    }
    if (cmp_bool(fp, "bV", -1, fr1->bV, fr2->bV))
    {
        cmp_rvecs(fp, "v", std::min(fr1->natoms, fr2->natoms), fr1->v, fr2->v, bRMSD, ftol, abstol);
    }
    if (cmp_bool(fp, "bF", -1, fr1->bF, fr2->bF))
    {
        cmp_rvecs(fp, "f", std::min(fr1->natoms, fr2->natoms), fr1->f, fr2->f, bRMSD, ftol, abstol);
    }
    if (cmp_bool(fp, "bBox", -1, fr1->bBox, fr2->bBox))
    {
        cmp_rvecs(fp, "box", 3, fr1->box, fr2->box, FALSE, ftol, abstol);
    }
}