/* ----------------------------------------------------------------------- *
 * This file is part of GEL, http://www.imm.dtu.dk/GEL
 * Copyright (C) the authors and DTU Informatics
 * For license and list of authors, see ../../doc/intro.pdf
 * ----------------------------------------------------------------------- */

#include "dual.h"

#include <vector>
#include "../CGLA/Vec3d.h"

#include "Manifold.h"
#include "AttributeVector.h"

namespace HMesh
{
    using namespace std;
    using namespace CGLA;
    
    void dual(Manifold& m)
    {
        // make sure every face knows its number
        int i = 0;
        
        FaceAttributeVector<int> ftouched;
        for(FaceIDIterator f = m.faces_begin(); f != m.faces_end(); ++f, ++i)
            ftouched[*f] = i;
        
        vector<Vec3d> vertices;
        vertices.resize(m.no_faces());
        vector<int> faces;
        vector<int> indices;
        
        // Create new vertices. Each face becomes a vertex whose position
        // is the centre of the face
        i = 0;
        for(FaceIDIterator f = m.faces_begin(); f != m.faces_end(); ++f, ++i)
            vertices[i] = centre(m, *f);
            
            // Create new faces. Each vertex is a new face with N=valency of vertex
            // edges.
            for(VertexIDIterator v = m.vertices_begin(); v!= m.vertices_end(); ++v){
                if(boundary(m, *v))
                    continue;
                
                Walker w = m.walker(*v);
                vector<int> index_tmp;
                for(; !w.full_circle(); w = w.circulate_vertex_ccw())
                    indices.push_back(ftouched[w.face()]);
                
                // Insert face valency in the face vector.
                faces.push_back(w.no_steps());
                
            }
            
            // Clear the manifold before new geometry is inserted.
            m.clear();
            
            // And build
            m.build(    vertices.size(), 
                    reinterpret_cast<double*>(&vertices[0]),
                    faces.size(), 
                    &faces[0],
                    &indices[0]);
        }
    }