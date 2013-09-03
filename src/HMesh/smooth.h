/**
 * @file smooth.h
 * @brief Functions for mesh smoothing.
 */

/* ----------------------------------------------------------------------- *
 * This file is part of GEL, http://www.imm.dtu.dk/GEL
 * Copyright (C) the authors and DTU Informatics
 * For license and list of authors, see ../../doc/intro.pdf
 * ----------------------------------------------------------------------- */

#ifndef __HMESH_SMOOTH_H
#define __HMESH_SMOOTH_H

namespace HMesh
{
    class Manifold;

    /// Simple laplacian smoothing with an optional weight.
    void laplacian_smooth(HMesh::Manifold& m, float t=1.0f);

    /// Taubin smoothing is similar to laplacian smoothing but reduces shrinkage
    void taubin_smooth(HMesh::Manifold& m, int iter);

    /// Fuzzy vector median smoothing is effective when it comes to preserving sharp edges.
    void fvm_smooth(HMesh::Manifold& m, int iter);
}
#endif