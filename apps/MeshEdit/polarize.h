//
//  polarize.h
//  GEL
//
//  Created by J. Andreas Bærentzen on 18/03/12.
//  Copyright 2012 __MyCompanyName__. All rights reserved.
//

#include <CGLA/Vec3d.h>
#include <HMesh/Manifold.h>
#include <HMesh/AttributeVector.h>

#ifndef __POLARIZE_H__
#define __POLARIZE_H__

void polarize_mesh(HMesh::Manifold& m, HMesh::VertexAttributeVector<double>& fun, double vmin, double vmax, const int divisions, HMesh::VertexAttributeVector<double>& parametrization);

void make_height_fun(const HMesh::Manifold& m, HMesh::VertexAttributeVector<double>& fun, 
                     double& vmin, double &vmax);

#endif
