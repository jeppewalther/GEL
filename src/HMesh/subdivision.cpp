/* ----------------------------------------------------------------------- *
 * This file is part of GEL, http://www.imm.dtu.dk/GEL
 * Copyright (C) the authors and DTU Informatics
 * For license and list of authors, see ../../doc/intro.pdf
 * ----------------------------------------------------------------------- */

#include "subdivision.h"

#include <vector>
#include "../CGLA/Vec3d.h"

#include "Manifold.h"
#include "AttributeVector.h"

namespace HMesh
{
    using namespace std;
    using namespace CGLA;

    void cc_split(Manifold& m_in, Manifold& m_out)
    {
        const int Invalid = -1;

        vector<Vec3d> new_points;
        new_points.reserve(m_in.no_vertices());

        VertexAttributeVector<int> vtouched(m_in.allocated_vertices(), Invalid);
        HalfEdgeAttributeVector<int> htouched(m_in.allocated_halfedges(), Invalid);

        int npsize = 0;
        for(VertexIDIterator v = m_in.vertices_begin(); v != m_in.vertices_end(); ++v){       
            vtouched[*v] = npsize;
            new_points.push_back(m_in.pos(*v));
            ++npsize;
        }

        for(HalfEdgeIDIterator h = m_in.halfedges_begin(); h != m_in.halfedges_end(); ++h){
            if(htouched[*h] != Invalid)
                continue;

            Walker w = m_in.walker(*h);
            htouched[*h] = htouched[w.opp().halfedge()] = npsize;
            new_points.push_back((m_in.pos(w.vertex()) + m_in.pos(w.opp().vertex())) * 0.5f);
            ++npsize;

        }
        vector<int> indices;
        vector<int> faces;

        for(FaceIDIterator f = m_in.faces_begin(); f != m_in.faces_end(); ++f){           
            for(Walker w = m_in.walker(*f); !w.full_circle(); w = w.circulate_face_cw()){
                indices.push_back(npsize);
                indices.push_back(htouched[w.halfedge()]);
                indices.push_back(vtouched[w.vertex()]);
                indices.push_back(htouched[w.next().halfedge()]);
                faces.push_back(4);
            }
            new_points.push_back(centre(m_in, *f));
            ++npsize;
        }

        m_out.clear();
        m_out.build(npsize, reinterpret_cast<double*>(&new_points[0]), faces.size(), &faces[0], &indices[0]);
    }
    
    void cc_smooth(Manifold& m)
    {
        VertexAttributeVector<Vec3d> new_vertices(m.no_vertices(), Vec3d(0));
        for(FaceIDIterator fi = m.faces_begin(); fi != m.faces_end(); ++fi)
        {				
            FaceID f = *fi;
            Walker w = m.walker(f);
            for(; !w.full_circle(); w = w.next())
            {
                VertexID v = w.vertex();
                float val = valency(m, v);
                float A = (1.0f-3.0f/val)	* (1.0f/val);
                float B = sqr(1.0f/val);
                Walker w2 = m.walker(f);
                for(; !w2.full_circle(); w2 = w2.next())
                {
                    VertexID v2 = w2.vertex();
                    if(v==v2)
                        new_vertices[v] += A * m.pos(v2);
                    else
                        new_vertices[v] += B * m.pos(v2);
                }
                
            }
        }
        for(VertexIDIterator vi = m.vertices_begin(); vi != m.vertices_end(); ++vi)
            m.pos(*vi) = new_vertices[*vi];
    }


}
