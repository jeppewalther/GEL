/* ----------------------------------------------------------------------- *
 * This file is part of GEL, http://www.imm.dtu.dk/GEL
 * Copyright (C) the authors and DTU Informatics
 * For license and list of authors, see ../../doc/intro.pdf
 * ----------------------------------------------------------------------- */

#include "smooth.h"

#include <vector>
#include <algorithm>
#include "../CGLA/Mat3x3f.h"
#include "../CGLA/Vec3d.h"

#include "Manifold.h"
#include "AttributeVector.h"

namespace HMesh
{
    using namespace std;
    using namespace CGLA;

    Vec3d laplacian(const Manifold& m, VertexID v)
    {
        Vec3d avg_pos(0);
        int n = 0;

        for(Walker w = m.walker(v); !w.full_circle(); w = w.circulate_vertex_cw()){
            avg_pos += m.pos(w.vertex());
            ++n;
        }
        return avg_pos / n - m.pos(v);
    }

    void laplacian_smooth(Manifold& m, float t)
    {
        VertexAttributeVector<Vec3d> pos(m.allocated_vertices());

        for(VertexIDIterator v = m.vertices_begin(); v != m.vertices_end(); ++v){
            if(!boundary(m, *v))
                pos[*v] =  t * laplacian(m, *v) + m.pos(*v);
        }

        for(VertexIDIterator v = m.vertices_begin(); v != m.vertices_end(); ++v){
            if(!boundary(m, *v))
                m.pos(*v) = pos[*v];
        }
    }

    void taubin_smooth(Manifold& m, int max_iter)
    {
        for(int iter = 0; iter < max_iter; ++iter) {
            VertexAttributeVector<Vec3d> lap(m.allocated_vertices());

            for(VertexIDIterator v = m.vertices_begin(); v != m.vertices_end(); ++v){
                if(!boundary(m, *v))
                    lap[*v] =  laplacian(m, *v);
            }

            for(VertexIDIterator v = m.vertices_begin(); v != m.vertices_end(); ++v){
                if(!boundary(m, *v))
                    m.pos(*v) += (iter%2 == 0) ? 0.5  * lap[*v] : -0.52 * lap[*v];
            }
        }
    }

    void face_neighbourhood(Manifold& m, FaceAttributeVector<int>& touched, FaceID f, vector<Vec3d>& nbrs)
    {	
        nbrs.push_back(normal(m, f));
        for(Walker wf = m.walker(f); !wf.full_circle(); wf = wf.circulate_face_cw()){
            for(Walker wv = m.walker(wf.vertex()); !wv.full_circle(); wv = wv.circulate_vertex_cw()){
                FaceID fn = wv.face();
                if(fn != InvalidFaceID && touched[fn] != touched[f]){
                    Vec3d norm = normal(m, fn);
                    if(!isnan(sqr_length(norm)))
                        nbrs.push_back(norm);
                    else
                        cout << "bad normal detected" << endl;
                    touched[fn] = touched[f];
                }
            }
        }
    }



    Vec3d filtered_normal(Manifold& m,  FaceAttributeVector<int>& touched, FaceID f)
    {
        const float sigma = .1f;

        vector<Vec3d> norms;
        face_neighbourhood(m, touched, f, norms);
        float min_dist_sum=1e32f;
        long int median=-1;

        for(size_t i=0;i<norms.size();++i)
        {
            float dist_sum = 0;
            for(size_t j=0;j<norms.size(); ++j)
                dist_sum += 1.0f - dot(norms[i], norms[j]);
            if(dist_sum < min_dist_sum)
            {
                min_dist_sum = dist_sum;
                median = i;
            }
        }
        assert(median != -1);
        Vec3d median_norm = norms[median];
        Vec3d avg_norm(0);
        for(size_t i=0;i<norms.size();++i)
        {
            float w = exp((dot(median_norm, norms[i])-1)/sigma);
            if(w<1e-2) w = 0;
            avg_norm += w*norms[i];
        }
        return normalize(avg_norm);
    }

    void fvm_smooth(HMesh::Manifold& m, int max_iter)
    {
        for(int iter = 0;iter<max_iter; ++iter)
        {
            FaceAttributeVector<int> touched(m.allocated_faces(), -1);
            FaceAttributeVector<Vec3d> filtered_norms(m.allocated_faces());

            int i = 0;
            for(FaceIDIterator f = m.faces_begin(); f != m.faces_end(); ++f,++i){
                touched[*f] = i;
                filtered_norms[*f] = filtered_normal(m, touched, *f);
            }

            VertexAttributeVector<Vec3d> vertex_positions(m.allocated_vertices());

            for(int iter=0;iter<20;++iter)
            {
                for(VertexIDIterator v = m.vertices_begin(); v != m.vertices_end(); ++v){
                    Vec3d move(0);
                    for(Walker w = m.walker(*v); !w.full_circle(); w = w.circulate_vertex_cw()){
                        FaceID f1 = w.face();
                        FaceID f2 = w.opp().face();
                        Vec3d dir = m.pos(w.vertex()) - m.pos(*v);

                        if(f1 != InvalidFaceID){
                            Vec3d n1 = filtered_norms[f1];
                            move += 0.05 * n1 * dot(n1, dir);
                        }
                        if(f2 != InvalidFaceID)
                        {
                            Vec3d n2 = filtered_norms[f2];
                            move += 0.05 * n2 * dot(n2, dir);
                        }
                    }
                    vertex_positions[*v] = m.pos(*v) + move;
                }
                for(VertexIDIterator v = m.vertices_begin(); v != m.vertices_end(); ++v)
                    m.pos(*v) = vertex_positions[*v];
            }
        }
    }


}
