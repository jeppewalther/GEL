/* ----------------------------------------------------------------------- *
 * This file is part of GEL, http://www.imm.dtu.dk/GEL
 * Copyright (C) the authors and DTU Informatics
 * For license and list of authors, see ../../doc/intro.pdf
 * ----------------------------------------------------------------------- */

#include "refine_edges.h"

#include <vector>
#include <iterator>

#include "Manifold.h"
#include "AttributeVector.h"

namespace HMesh
{
    using namespace std;

    float average_edge_length(const Manifold& m)
    {
        float lsum = 0;
        for(HalfEdgeIDIterator h = m.halfedges_begin(); h != m.halfedges_end(); ++h)
            lsum += length(m, *h);
        return lsum / m.no_halfedges();
    }

    int refine_edges(Manifold& m, float t)
    {
        int work = 0;
       
        vector<HalfEdgeID> hedges;
        hedges.reserve(m.no_halfedges());
        copy(m.halfedges_begin(), m.halfedges_end(), back_inserter(hedges));

        HalfEdgeAttributeVector<int> touched(m.allocated_halfedges(), 0);

        cout << "Refining edges";
        for(vector<HalfEdgeID>::iterator h = hedges.begin(); h != hedges.end(); ++h){
            Walker w = m.walker(*h);

            if(!m.in_use(*h) || w.face() == InvalidFaceID || length(m, *h) < t || touched[*h])
                continue;

            ++work;
            if( work % 10000 == 0)
                cout << ".";

            touched[w.opp().halfedge()] = 1;
//            VertexID v = m.split_edge(*h);

//            Walker wv = m.walker(v);

//            FaceID f1 = wv.opp().face();
           // if(f1 != InvalidFaceID)
             //   m.split_face_by_vertex(f1);

//            FaceID f2 =  wv.prev().face();
           // if(f2 != InvalidFaceID)
            //    m.split_face_by_vertex(f2);

        }
        cout << endl;
        return work;
    }

}
