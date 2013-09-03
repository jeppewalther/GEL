//
//  polarize.cpp
//  GEL
//
//  Created by J. Andreas Bærentzen on 18/03/12.
//  Copyright 2012 __MyCompanyName__. All rights reserved.
//
#include <queue>

#include "polarize.h"
#include <HMesh/triangulate.h>
#include <HMesh/curvature.h>
#include <HMesh/quadric_simplify.h>
#include <HMesh/mesh_optimization.h>
#include <HMesh/smooth.h>

#include <complex>

using namespace CGLA;
using namespace std;
using namespace HMesh;

void shortest_edge_triangulate_face(Manifold& m, FaceID f0, VertexAttributeVector<int>& level_set_id_vertex)
{
    queue<FaceID> face_queue;
    
    face_queue.push(f0);
    
    while(!face_queue.empty())
    {
        FaceID f = face_queue.front();
        face_queue.pop();
        
        // Create a vector of vertices.
        vector<VertexID> verts;
        for(Walker w = m.walker(f); !w.full_circle(); w = w.circulate_face_ccw())
        {
            FaceID fa = w.face();
            FaceID fb = f;
            assert(fa==fb);
            verts.push_back(w.vertex());
        }
        // If there are just three we are done.
        if(verts.size() == 3) continue;
        
        // Find vertex pairs that may be connected.
        vector<pair<int,int> > vpairs;
        const int N = verts.size();
        for(int i = 0; i < N - 2; ++i){
            for(int j = i + 2; j < N; ++j){
                if(verts[i] != verts[j] && 
                   !connected(m, verts[i], verts[j]) &&
                   (level_set_id_vertex[verts[i]] == 0 || level_set_id_vertex[verts[i]] != level_set_id_vertex[verts[j]])
                   )
                    vpairs.push_back(pair<int,int>(i, j));
            }
        }
        if(vpairs.empty()){
            cout << "Warning: could not triangulate a face." 
            << "Probably a vertex appears more than one time in other vertex's one-ring" << endl;
            continue;
        }
        
        /* For all vertex pairs, find the edge lengths. Combine the
         vertices forming the shortest edge. */
        
        float min_len=FLT_MAX;
        int min_k = -1;
        for(size_t k = 0; k < vpairs.size(); ++k){
            int i = vpairs[k].first;
            int j = vpairs[k].second;
            float len = sqr_length(m.pos(verts[i]) - m.pos(verts[j]));
            
            if(len<min_len){
                min_len = len;
                min_k = k;
            }
        }
        assert(min_k != -1);
        
        {
            // Split faces along edge whose midpoint is closest to isovalue
            int i = vpairs[min_k].first;
            int j = vpairs[min_k].second;
            FaceID f_new = m.split_face_by_edge(f, verts[i], verts[j]);
            
            if(no_edges(m, f)>3)
                face_queue.push(f);
            if(no_edges(m, f_new)>3)
                face_queue.push(f_new);
        }
        
    }
}


void shortest_edge_split_face(Manifold& m, FaceID f0, VertexAttributeVector<int>& level_set_id_vertex)
{
    queue<FaceID> face_queue;
    
    face_queue.push(f0);
    
    while(!face_queue.empty())
    {
        FaceID f = face_queue.front();
        face_queue.pop();
        
        // Create a vector of vertices.
        vector<VertexID> verts;
        for(Walker w = m.walker(f); !w.full_circle(); w = w.circulate_face_ccw())
        {
            verts.push_back(w.vertex());
        }
        
        // Find vertex pairs that may be connected.
        vector<pair<int,int> > vpairs;
        const int N = verts.size();
        for(int i = 0; i < N ; ++i){
            for(int j = 3; j < N-2; ++j){
                int jj = (j+i)%N;
                if(!connected(m, verts[i], verts[jj]) &&
                   (level_set_id_vertex[verts[i]] != level_set_id_vertex[verts[jj]]) &&
                   (level_set_id_vertex[verts[(i+1)%N]] == level_set_id_vertex[verts[i]]) &&
                   (level_set_id_vertex[verts[(i+N-1)%N]] == level_set_id_vertex[verts[i]]) &&
                   (level_set_id_vertex[verts[(jj+1)%N]] == level_set_id_vertex[verts[jj]]) &&
                   (level_set_id_vertex[verts[(jj+N-1)%N]] == level_set_id_vertex[verts[jj]]))
                    vpairs.push_back(pair<int,int>(i, jj));
            }
        }
        if(vpairs.empty()) continue;
        
        
        /* For all vertex pairs, find the edge lengths. Combine the
         vertices forming the shortest edge. */
        
        float min_len=FLT_MAX;
        int min_k = -1;
        for(size_t k = 0; k < vpairs.size(); ++k){
            int i = vpairs[k].first;
            int j = vpairs[k].second;
            float len = sqr_length(m.pos(verts[i]) - m.pos(verts[j]));
            
            if(len<min_len){
                min_len = len;
                min_k = k;
            }
        }
        assert(min_k != -1);
        
        {
            // Split faces along edge whose midpoint is closest to isovalue
            int i = vpairs[min_k].first;
            int j = vpairs[min_k].second;
            FaceID f_new = m.split_face_by_edge(f, verts[i], verts[j]);
            if(no_edges(m, f)>5)
                face_queue.push(f);
            if(no_edges(m, f_new)>5)
                face_queue.push(f_new);
        }
        
    }
}



struct EdgeQElem {
    float priority;
    HalfEdgeID he;
    EdgeQElem(float _priority, HalfEdgeID _he): priority(_priority), he(_he) {}
};

bool operator<(const EdgeQElem& l, const EdgeQElem& r)
{
    return l.priority < r.priority;
}

class FunctionalDifference: public EnergyFun
{
    VertexAttributeVector<float>& fun;
    VertexAttributeVector<int>& status;
public:
    FunctionalDifference(VertexAttributeVector<float>& _fun, VertexAttributeVector<int>& _status): fun(_fun), status(_status) {}
    virtual double delta_energy(const Manifold& m, HalfEdgeID h) const 
    {
        Walker w = m.walker(h);
        if(status[w.vertex()] == 1 && status[w.opp().vertex()]==1)
            return DBL_MAX;
        return sqr_length(m.pos(w.next().vertex())-m.pos(w.opp().next().vertex()))/(1e-6+abs(fun[w.next().vertex()]-fun[w.opp().next().vertex()])) - sqr_length(m.pos(w.vertex())-m.pos(w.opp().vertex()))/(1e-6+abs(fun[w.vertex()]-fun[w.opp().vertex()]));
    }
};

class TriangleQualityValence: public EnergyFun
{
    VertexAttributeVector<int>& idv;
    MinAngleEnergy mae;
    ValencyEnergy vae;
public:
    TriangleQualityValence(VertexAttributeVector<int>& _idv): idv(_idv), mae(-1) {}
    virtual double delta_energy(const Manifold& m, HalfEdgeID h) const 
    {
        Walker w = m.walker(h);
        if(idv[w.next().vertex()] == idv[w.opp().next().vertex()] || 
           idv[w.vertex()] == idv[w.opp().vertex()])
            return DBL_MAX;
        
        VertexID v1 = w.opp().vertex();
        VertexID v2 = w.vertex();
        VertexID vo1 = w.next().vertex();
        VertexID vo2 = w.opp().next().vertex();
        
        int val1  = valency(m, v1);
        int val2  = valency(m, v2);
        int valo1 = valency(m, vo1);
        int valo2 = valency(m, vo2);
        
        // The optimal valency is four for a boundary vertex
        // and six elsewhere.
        int t1 = boundary(m, v1) ? 4 : 6;
        int t2 = boundary(m, v2) ? 4 : 6;
        int to1 = boundary(m, vo1) ? 4 : 6;
        int to2 = boundary(m, vo2) ? 4 : 6;
        
        int before = 
        max(max(sqr(val1-t1),sqr(val2-t2)), max(sqr(valo1-to1), sqr(valo2-to2)));
        int after = 
        max(max(sqr(valo1+1-to1),sqr(val1-1-t1)), max(sqr(val2-1-t2),sqr(valo2+1-to2)));
        
        return static_cast<double>(after-before);
        
        //        return vae.delta_energy(m,h);
        float la = length(m.pos(w.next().vertex())-m.pos(w.opp().next().vertex()));
        float lb = length(m.pos(w.vertex())-m.pos(w.opp().vertex()));
        return la-lb;
        return mae.delta_energy(m,h);
    }
};

class TriangleQuality: public EnergyFun
{
    VertexAttributeVector<int>& idv;
    MinAngleEnergy mae;
    ValencyEnergy vae;
public:
    TriangleQuality(VertexAttributeVector<int>& _idv): idv(_idv), mae(-1) {}
    virtual double delta_energy(const Manifold& m, HalfEdgeID h) const 
    {
        Walker w = m.walker(h);
        if(idv[w.next().vertex()] == idv[w.opp().next().vertex()] || 
           idv[w.vertex()] == idv[w.opp().vertex()])
            return DBL_MAX;
        return mae.delta_energy(m,h);
    }
};

Vec3d grad(HMesh::Manifold& m, HMesh::VertexAttributeVector<double>& fun, HMesh::FaceID f)
{
    if(no_edges(m,f) != 3)
        return Vec3d(0.0);
    
    
    Vec3d n = normal(m, f);
    
    Vec3d gsum(0.0);
    for(Walker w = m.walker(f); !w.full_circle(); w = w.next())
    {
        Vec3d gdir = normalize(cross(n, m.pos(w.vertex()) - m.pos(w.opp().vertex())));
        float l = dot(m.pos(w.next().vertex())-m.pos(w.vertex()), gdir);
        gdir *= fun[w.next().vertex()]/l;
        gsum += gdir;
    }
    return gsum;
}

complex<double> complex_slerp(double t, const complex<double>& a, const complex<double>&b)
{
    double omega = (arg(a/b));
    cout << omega << endl;
    return (a*sin((1-t)*omega) + b*sin(t*omega))/sin(omega);
}

double mod_2pi(double x)
{
    return x-floor(x/(2.0 * M_PI))*2.0 * M_PI;
}

void extend_fun2(HMesh::Manifold& m, HMesh::HalfEdgeID h,
                 HMesh::VertexAttributeVector<double>& fun, 
                 HMesh::VertexAttributeVector<double>& fun2)
{
    Walker w = m.walker(h);
    FaceID f = w.face();
    Vec3d a = m.pos(w.opp().vertex());
    Vec3d b = m.pos(w.vertex());
    Vec3d c = m.pos(w.next().vertex());
    Vec3d g = grad(m, fun, f);
    Vec3d n = normal(m, f);
    Vec3d N = -normalize(cross(g, n));
    float dot_nba = dot(N,b-a);
    float v;
    if(dot_nba > 1e-5)
    {
        float t = dot(N, c-a)/dot_nba;
        double aval = fun2[w.opp().vertex()];
        double bval = fun2[w.vertex()];
        double val = (1-t)*aval + t*bval;
        v = mod_2pi(val);
        cout << t << " , " << v << endl; 
    }
    else
        v =  fun2[w.vertex()];
    fun2[w.next().vertex()] = v;
}

inline bool same_level(float a, float b) {return abs(a-b) < 0.00001;}

Vec3d laplacian(const Manifold& m, VertexID v)
{
    Vec3d avg_pos(0);
    float asum = 0;
    
    for(Walker w = m.walker(v); !w.full_circle(); w = w.circulate_vertex_cw()){
        float  a = barycentric_area(m, w.vertex());
        avg_pos += a * m.pos(w.vertex());
        asum += a;
    }
    return avg_pos / asum - m.pos(v);
}

void aw_laplacian_smooth(Manifold& m, float t)
{
    VertexAttributeVector<Vec3d> pos(m.allocated_vertices());
    
    for(VertexIDIterator v = m.vertices_begin(); v != m.vertices_end(); ++v){
        if(!boundary(m, *v))
        {
            Vec3d n = normal(m, *v);
            Vec3d l = laplacian(m, *v);
            if(sqr_length(n) > 0.8)
                l -= n * dot(n,l);
            pos[*v] =  t * l + m.pos(*v);
        }
    }
    
    for(VertexIDIterator v = m.vertices_begin(); v != m.vertices_end(); ++v){
        if(!boundary(m, *v))
            m.pos(*v) = pos[*v];
    }
}

void polarize_mesh(Manifold& m, VertexAttributeVector<double>& fun, double vmin, double vmax, const int divisions, VertexAttributeVector<double>& parametrization)
{
    vmax -= 0.01 * (vmax-vmin);
    float interval = (vmax-vmin)/divisions;
    
    VertexAttributeVector<int> status(m.allocated_vertices(), 0);
    
    
    // ----------------------------------------
    cout << "Tracing level set curves" << endl;
    
    vector<HalfEdgeID> hidvec;
    for(HalfEdgeIDIterator hid = m.halfedges_begin(); hid != m.halfedges_end(); ++hid)
    {
        Walker w = m.walker(*hid);
        if(fun[w.vertex()] > fun[w.opp().vertex()])
            hidvec.push_back(*hid);
    }
    
    for(size_t i = 0; i<hidvec.size(); ++i)
    {
        Walker w = m.walker(hidvec[i]);
        
        float b = (fun[w.vertex()]- vmin)/interval;
        float a = (fun[w.opp().vertex()] - vmin)/interval;
        float floor_b = floor(b);
        float floor_a = floor(a);
        
        Vec3d pb = m.pos(w.vertex());
        for(int j=floor_b; j>floor_a; --j)
        {
            float t = (j-a) / (b-a);
            Vec3d p = t * pb + (1.0-t) * m.pos(w.opp().vertex());
            VertexID v_new = m.split_edge(w.halfedge());
            w = w.prev();
            status[v_new] = 1;
            fun[v_new] = j * interval + vmin;
            m.pos(v_new) = p;
        }
    }
    
    bool did_work;
    do
    {
        did_work = false;
        
        for(FaceIDIterator fid = m.faces_begin(); fid != m.faces_end(); ++fid)
            for(Walker w = m.walker(*fid);!w.full_circle(); w = w.next())
                if(status[w.vertex()] == 1 && !(status[w.next().vertex()]==1 && same_level(fun[w.vertex()],fun[w.next().vertex()]))
                   && !(status[w.prev().vertex()]==1 && same_level(fun[w.vertex()],fun[w.prev().vertex()])))
                {
                    Walker w0 = w;
                    w = w.next().next();
                    do
                    {
                        if(status[w.vertex()] == 1 && w.next().halfedge() != w0.halfedge() &&
                           same_level(fun[w0.vertex()],fun[w.vertex()]))
                        {
                            m.split_face_by_edge(*fid, w0.vertex(), w.vertex());
                            did_work = true;
                            break;
                        }
                        w = w.next();
                    }
                    while(!w.full_circle());
                    break;
                }
    }
    while(did_work);
    
    
    
    
    // ----------------------------
    cout << "Numbering the level sets" << endl;
    
    float max_length = 0;
    double max_length_fun = 0;
    int max_length_id =-1;
    HalfEdgeID max_length_h;
    HalfEdgeAttributeVector<int> level_set_id(m.allocated_halfedges(), 0);
    VertexAttributeVector<int> level_set_id_vertex(m.allocated_vertices(), 0);
    int no_id=1;
    for(HalfEdgeIDIterator hid = m.halfedges_begin(); hid != m.halfedges_end(); ++hid)
    {
        Walker w = m.walker(*hid);
        if(status[w.vertex()] == 1 && status[w.opp().vertex()] == 1 &&
           same_level(fun[w.vertex()], fun[w.opp().vertex()]) &&
           level_set_id[w.halfedge()] == 0)
        {
            float level_set_length = 0;
            while(level_set_id[w.halfedge()] != no_id)
            {
                level_set_length += length(m,w.halfedge());
                level_set_id[w.halfedge()] = no_id;
                level_set_id[w.opp().halfedge()] = no_id;
                level_set_id_vertex[w.vertex()] = no_id;
                w = w.next();
                while(status[w.vertex()] != 1 || !same_level(fun[w.vertex()], fun[w.opp().vertex()]))
                    w = w.circulate_vertex_cw();
            }            
            if(level_set_length > max_length)
            {
                max_length = level_set_length;
                max_length_id = no_id;
                max_length_h = w.halfedge();
                max_length_fun = fun[w.vertex()];
            }
            ++no_id;
        }
    }
    cout << max_length << " " << max_length_id << endl;
    cout << "Number of level sets : " << (no_id-1);
    // ----------------------------
    
    shortest_edge_triangulate(m);
    
    
    // -------------
    cout << "Parametrize level sets " << endl;
    
    queue<HalfEdgeID> hq;
    HalfEdgeAttributeVector<int> touched(m.no_halfedges(), 0);
    Walker w = m.walker(max_length_h);
    float param = 0;
    do
    {
        parametrization[w.opp().vertex()] = 2.0 * M_PI * param / max_length;
        param += length(m, w.halfedge());
        w = w.next();
        while(level_set_id[w.halfedge()] != max_length_id)
            w = w.circulate_vertex_cw();
        hq.push(w.halfedge());
        hq.push(w.opp().halfedge());
        touched[w.halfedge()] = 1;
        touched[w.opp().halfedge()] = 1;
    }            
    while(w.halfedge() != max_length_h);
    
    while(!hq.empty())
    {
        HalfEdgeID h = hq.front();
        hq.pop();
        if(!touched[w.next().halfedge()])
        {
            extend_fun2(m, h, fun, parametrization);
            touched[w.next().halfedge()] = 1;
            touched[w.prev().halfedge()] = 1;
            Walker w = m.walker(h);
            if(!touched[w.next().opp().halfedge()])
            {
                touched[w.next().opp().halfedge()] = 1;
                hq.push(w.next().opp().halfedge());
            }
            if(!touched[w.prev().opp().halfedge()])
            {
                touched[w.prev().opp().halfedge()] = 1;
                hq.push(w.prev().opp().halfedge());
                
            }
        }
    }
    
    return; 
    
    // ----------------------------
    cout << "Remove vertices not on level set curves" << endl;
    
    vector<VertexID> vid_vec;
    for(VertexIDIterator vid = m.vertices_begin(); vid != m.vertices_end(); ++vid)
        if(status[*vid]==0)
            vid_vec.push_back(*vid);
    
    random_shuffle(vid_vec.begin(), vid_vec.end());
    for (size_t i=0; i<vid_vec.size(); ++i) {
        FaceID f = m.merge_one_ring(vid_vec[i]);
        if(f != InvalidFaceID)
            shortest_edge_triangulate_face(m, f, level_set_id_vertex);
        else
            cout << "vertex not removed " << valency(m, vid_vec[i]) << endl; 
    }
    
    for(FaceIDIterator fid = m.faces_begin(); fid != m.faces_end(); ++fid)
        if(no_edges(m, *fid) > 3)
            shortest_edge_triangulate_face(m, *fid, level_set_id_vertex);
    
    
    VertexAttributeVector<Vec3d> recalled_positions;
    for(VertexIDIterator vid = m.vertices_begin(); vid != m.vertices_end(); ++vid)
        recalled_positions[*vid] = m.pos(*vid);
    
    
    TriangleQuality tq_energy(level_set_id_vertex);
    priority_queue_optimization(m, tq_energy);
    
    
    
    //// ----------------------------
    cout << "smooth level set curves" << endl;
    
    
    for(int iter=0;iter< 14;++iter)
    {
        TriangleQualityValence tqv_energy(level_set_id_vertex);
        TriangleQuality tq_energy(level_set_id_vertex);
        priority_queue_optimization(m, tqv_energy);
        priority_queue_optimization(m, tq_energy);
        
        VertexAttributeVector<Vec3d> new_pos(m.allocated_vertices(), Vec3d(0));
        VertexAttributeVector<float> new_pos_wsum(m.allocated_vertices(), 0.0);
        for(HalfEdgeIDIterator hid = m.halfedges_begin(); hid != m.halfedges_end(); ++hid)
        {
            Walker w = m.walker(*hid);
            if(level_set_id[w.halfedge()] != 0)
            {
                float weight = 1.0;//sqr(valency(m,w.opp().vertex())-2);
                new_pos[w.vertex()] += weight*m.pos(w.opp().vertex());
                new_pos_wsum[w.vertex()] += weight;
            }
        }
        for(VertexIDIterator vid = m.vertices_begin(); vid != m.vertices_end(); ++vid)
        {
            float weight = 1.0;//sqr(valency(m,*vid)-2);
            new_pos[*vid] += weight*m.pos(*vid);
            new_pos_wsum[*vid] += weight;
            m.pos(*vid) = 0.75*m.pos(*vid) + 0.25 *new_pos[*vid] / new_pos_wsum[*vid];
        }
        
        priority_queue_optimization(m, tqv_energy);
        priority_queue_optimization(m, tq_energy);
        vector<HalfEdgeID> hidvec;
        for(HalfEdgeIDIterator hid = m.halfedges_begin(); hid != m.halfedges_end(); ++hid)
        {
            Walker w = m.walker(*hid);
            if(w.halfedge() < w.opp().halfedge() && 
               level_set_id[w.halfedge()] != 0 &&
               (
                ((level_set_id[w.opp().halfedge()] == level_set_id[w.opp().next().halfedge()] ||
                  level_set_id[w.halfedge()] == level_set_id[w.next().halfedge()]) &&
                 valency(m, w.vertex())+valency(m,w.opp().vertex())>10) ||
                
                valency(m, w.vertex())+valency(m,w.opp().vertex())>14
                )
               )
                hidvec.push_back(w.halfedge());
        }
        
        for(size_t i=0;i<hidvec.size(); ++i)
        {
            Walker w = m.walker(hidvec[i]);
            VertexID v = m.split_edge(hidvec[i]);
            level_set_id_vertex[v] = level_set_id[w.halfedge()];
            level_set_id[w.prev().halfedge()] = level_set_id[w.halfedge()];
            level_set_id[w.prev().opp().halfedge()] = level_set_id[w.halfedge()];
            shortest_edge_triangulate_face(m, w.face(), level_set_id_vertex);
            shortest_edge_triangulate_face(m, w.opp().face(), level_set_id_vertex);
        }
        
        priority_queue_optimization(m, tqv_energy);
        priority_queue_optimization(m, tq_energy);
        
        for(HalfEdgeIDIterator hid = m.halfedges_begin(); hid != m.halfedges_end(); ++hid)
        {
            Walker w = m.walker(*hid);
            if(level_set_id[w.halfedge()] != 0 &&
               !(level_set_id[w.opp().halfedge()] == level_set_id[w.opp().next().halfedge()] ||
                 level_set_id[w.halfedge()] == level_set_id[w.next().halfedge()])  &&
               valency(m, w.vertex())<5 &&
               valency(m, w.opp().vertex())<5 &&
               precond_collapse_edge(m, w.halfedge()))
            {
                m.collapse_edge(w.halfedge(), true);
                did_work = true;
            }
        }
        
        
    }
    
    return;
    
    priority_queue<EdgeQElem> edge_queue;
    for(HalfEdgeIDIterator hid = m.halfedges_begin(); hid != m.halfedges_end(); ++hid)
    {
        Walker w = m.walker(*hid);
        if(w.halfedge()<w.opp().halfedge() &&
           level_set_id[w.halfedge()] == 0)
        {
            Vec3d v = normalize(m.pos(w.vertex()) - m.pos(w.opp().vertex()));
            
            float weight = 0;
            Walker wo = m.walker(w.opp().vertex());
            for(; !w.full_circle(); w = w.circulate_vertex_ccw())
                if(level_set_id[w.halfedge()] != 0)
                {
                    Vec3d e = normalize(m.pos(w.vertex()) - m.pos(w.opp().vertex()));
                    weight += abs(dot(v,e));
                }
            for(; !wo.full_circle(); wo = wo.circulate_vertex_ccw())
                if(level_set_id[wo.halfedge()] != 0)
                {
                    Vec3d e = normalize(m.pos(wo.vertex()) - m.pos(wo.opp().vertex()));
                    weight += abs(dot(v,e));
                }
            edge_queue.push(EdgeQElem(weight, *hid));        
        }
    }
    
    while(!edge_queue.empty())
    {
        HalfEdgeID h = edge_queue.top().he;
        edge_queue.pop();
        
        Walker w = m.walker(h);
        if(no_edges(m, w.face())==3 && no_edges(m,w.opp().face())==3 &&
           !(level_set_id[w.next().halfedge()] == level_set_id[w.opp().prev().halfedge()] ||
             level_set_id[w.prev().halfedge()] == level_set_id[w.opp().next().halfedge()]))
            m.merge_faces(w.face(), w.halfedge());
    }
    
    
    return;
    
    
    do{
        did_work = false;
        for(HalfEdgeIDIterator hid = m.halfedges_begin(); hid != m.halfedges_end(); ++hid)
        {
            Walker w = m.walker(*hid);
            
            if(level_set_id[w.halfedge()] != 0 && no_edges(m, w.face())==3 &&
               precond_collapse_edge(m, w.halfedge()))
            {
                m.collapse_edge(w.halfedge(), true);
                did_work = true;
            }
        }
    }
    while(did_work);
    
    
    return;
    
    
    int k=0;
    do {
        ++k;
        did_work = false;
        for(HalfEdgeIDIterator hid = m.halfedges_begin(); hid != m.halfedges_end(); ++hid)
        {
            Walker w0 = m.walker(*hid);
            
            if(level_set_id[w0.halfedge()] != 0 &&
               (level_set_id[w0.next().halfedge()] == 0 && level_set_id[w0.prev().halfedge()] == 0))
            {
                
                Walker w = w0;
                bool do_split = false;
                for(;!w.full_circle(); w = w.next())
                {
                    if(level_set_id[w.halfedge()] != 0 &&
                       (level_set_id[w.next().halfedge()] == level_set_id[w.halfedge()]))
                        do_split = true;
                }
                if(do_split)
                {
                    VertexID v = m.split_edge(w0.halfedge());
                    level_set_id_vertex[v] = level_set_id[w0.halfedge()];
                    level_set_id[w0.prev().halfedge()] = level_set_id[w0.halfedge()];
                    level_set_id[w0.prev().opp().halfedge()] = level_set_id[w0.halfedge()];
                    did_work = true;
                }
            }
            
            for(FaceIDIterator fid = m.faces_begin(); fid != m.faces_end(); ++fid)
                if(no_edges(m,*fid) >= 6)
                {
                    shortest_edge_split_face(m, *fid, level_set_id_vertex);
                    did_work = true;
                }
            
        }
    } while (did_work && k<1);
    
    
}

void make_height_fun(const HMesh::Manifold& m, HMesh::VertexAttributeVector<double>& fun,
                     double& vmin, double& vmax)
{
    VertexIDIterator vid = m.vertices_begin();
    vmin = vmax = m.pos(*vid)[2];
    for(; vid != m.vertices_end(); ++vid)
    {
        double v = m.pos(*vid)[1];
        fun[*vid] = v;
        vmin = min(v, vmin);
        vmax = max(v, vmax);
    }
}

//    //-------------------------
//    cout << "Remove short level set edges" << endl;
//    float avglen=0;
//    int n=0;
//    for(HalfEdgeIDIterator hid = m.halfedges_begin(); hid != m.halfedges_end(); ++hid)
//        {
//            Walker w = m.walker(*hid);
//            if(level_set_id[w.halfedge()] != 0)
//            {
//                avglen += length(m, w.halfedge());
//                ++n;
//            }
//        }
//    avglen /= n;
//    for(HalfEdgeIDIterator hid = m.halfedges_begin(); hid != m.halfedges_end(); ++hid)
//        if (length(m,*hid)<0.25*avglen && precond_collapse_edge(m, *hid)) {
//            m.collapse_edge(*hid);
//        }    


//        vector<HalfEdgeID> hidvec;
//        for(HalfEdgeIDIterator hid = m.halfedges_begin(); hid != m.halfedges_end(); ++hid)
//        {
//            Walker w = m.walker(*hid);
//            if(w.halfedge() < w.opp().halfedge() && 
//               level_set_id[w.halfedge()] != 0 &&
//               valency(m, w.vertex())+valency(m,w.opp().vertex())>12)
//                hidvec.push_back(w.halfedge());
//        }
//        
//        for(int i=0;i<hidvec.size(); ++i)
//        {
//            Walker w = m.walker(hidvec[i]);
//            VertexID v = m.split_edge(hidvec[i]);
//            level_set_id_vertex[v] = level_set_id[w.halfedge()];
//            level_set_id[w.prev().halfedge()] = level_set_id[w.halfedge()];
//            level_set_id[w.prev().opp().halfedge()] = level_set_id[w.halfedge()];
//            shortest_edge_triangulate_face(m, w.face(), level_set_id_vertex);
//            shortest_edge_triangulate_face(m, w.opp().face(), level_set_id_vertex);
//        }
//        
//        //     priority_queue_optimization(m, tq_energy);




//// ----------------------------
//cout << "smooth level set curves" << endl;
//
//for(int iter=0;iter<2;++iter)
//{
//    VertexAttributeVector<Vec3d> new_pos(m.allocated_vertices(), Vec3d(0));
//    for(HalfEdgeIDIterator hid = m.halfedges_begin(); hid != m.halfedges_end(); ++hid)
//    {
//        Walker w = m.walker(*hid);
//        if(status[w.vertex()] == 1 && status[w.opp().vertex()] == 1 &&
//           same_level(fun[w.vertex()], fun[w.opp().vertex()]))
//        {
//            new_pos[w.vertex()] += m.pos(w.vertex()) + m.pos(w.opp().vertex());
//        }
//    }
//    for(VertexIDIterator vid = m.vertices_begin(); vid != m.vertices_end(); ++vid)
//        if(status[*vid] == 1)
//            m.pos(*vid) = new_pos[*vid] / (4.0);
//}



//return;
//
//for (int iter=0; iter<1; ++iter) {
//    cout << __FILE__ << __LINE__ << endl;
//    cout << __FILE__ << __LINE__ << endl;
//    TriangleQuality tq_energy(level_set_id, level_set_id_vertex);
//    priority_queue_optimization(m, tq_energy);
//    //        FunctionalDifference energy(fun, status);
//    //        priority_queue_optimization(m, energy);
//    
//    // ----------------------------
//    cout << "Remove vertices not on level set curves" << endl;
//    
//    vector<VertexID> vid_vec;
//    for(VertexIDIterator vid = m.vertices_begin(); vid != m.vertices_end(); ++vid)
//    if(status[*vid]==0)
//    vid_vec.push_back(*vid);
//    
//    for (int i=0; i<vid_vec.size(); ++i) {
//        FaceID f = m.merge_one_ring(vid_vec[i]);
//    }
//    cout << "-" << endl;
//    }    
// --------------------------
// Triangulate polygons by connecting only vertices on different curves.

//    shortest_edge_triangulate(m);
//    TriangleQuality tq_energy(level_set_id);
//    priority_queue_optimization(m, tq_energy);


//        k=0;
//        do {
//            ++k;
//            did_work = false;
//            priority_queue<EdgeQElem> edge_queue;
//            
//            for(HalfEdgeIDIterator hid = m.halfedges_begin(); hid != m.halfedges_end(); ++hid)
//            {
//                Walker w = m.walker(*hid);
//                if(status[w.opp().vertex()] == 0 && status[w.vertex()] == 1)
//                {
//                    float weight = (abs(fun[w.vertex()] - fun[w.opp().vertex()]*length(m, w.halfedge())) + 1e-6);            
//                    edge_queue.push(EdgeQElem(weight, w.halfedge()));
//                }
//                
//                
//                while(!edge_queue.empty())
//                {
//                    HalfEdgeID he = edge_queue.top().he;
//                    edge_queue.pop();
//                    if(m.in_use(he))
//                    {
//                        if(precond_collapse_edge(m,he))
//                        {
//                            m.collapse_edge(he);
//                            did_work = true;
//                       }
//                    }
//                }
//            } 
//        }while (did_work && k < 100);
//            
//            cout << "k=" << k << endl;
//    priority_queue_optimization(m, energy);


//priority_queue<EdgeQElem> edge_queue;
//HalfEdgeAttributeVector<int> time_stamp(m.allocated_halfedges(), 0);
//for(HalfEdgeIDIterator hid = m.halfedges_begin(); hid != m.halfedges_end(); ++hid)
//if(level_set_id[*hid]==0)
//{
//    Walker w = m.walker(*hid);
//    if(w.halfedge()<w.opp().halfedge() && !(status[w.vertex()]==1 && status[w.opp().vertex()]==1))
//        edge_queue.push(EdgeQElem(-(sqr(fun[w.vertex()]-fun[w.opp().vertex()]))*sqr_length(m.pos(w.vertex())-m.pos(w.opp().vertex())),*hid,0));
//        }
//
//shortest_edge_triangulate(m);
//int k=0;
//while(!edge_queue.empty())
//{
//    HalfEdgeID hid = edge_queue.top().he;
//    if(m.in_use(hid) && time_stamp[hid]== edge_queue.top().time_stamp)
//    {
//        Walker w = m.walker(hid);
//        Walker wa = w.next();
//        Walker wb = w.prev().opp();
//        
//        if(m.merge_faces(w.face(), hid))
//        {
//            cout << ".";
//            if(valency(m, wa.opp().vertex())==2 && precond_collapse_edge(m, wa.halfedge()))
//            {
//                HalfEdgeID h = wa.halfedge();
//                wa = wa.prev();
//                m.collapse_edge(h, false);
//                ++time_stamp[wa.halfedge()];
//                ++time_stamp[wa.opp().halfedge()];
//                if(level_set_id[wa.halfedge()]==0 && !(status[wa.vertex()]==1 && status[wa.opp().vertex()]==1))
//                    edge_queue.push(EdgeQElem(-(sqr(fun[wa.vertex()]-fun[wa.opp().vertex()]))*sqr_length(m.pos(wa.vertex())-m.pos(wa.opp().vertex())),wa.halfedge(),time_stamp[wa.halfedge()]));
//            }
//            if(valency(m, wb.opp().vertex())==2 && precond_collapse_edge(m, wb.halfedge()))
//            {
//                HalfEdgeID h = wb.halfedge();
//                wb = wb.prev();
//                m.collapse_edge(h, false);
//                ++time_stamp[wb.halfedge()];
//                ++time_stamp[wb.opp().halfedge()];
//                if(level_set_id[wb.halfedge()]==0 && !(status[wb.vertex()]==1 && status[wb.opp().vertex()]==1))
//                    edge_queue.push(EdgeQElem(-(sqr(fun[wb.vertex()]-fun[wb.opp().vertex()]))*sqr_length(m.pos(wb.vertex())-m.pos(wb.opp().vertex())),wb.halfedge(),time_stamp[wb.halfedge()]));
//            }
//            
//            }
//            
//            }
//            
//            edge_queue.pop();
//            }
//

//   // bool did_work;
//    do{
//        did_work = false;
//    for(HalfEdgeIDIterator hid = m.halfedges_begin(); hid != m.halfedges_end(); ++hid)
//    {
//        Walker w = m.walker(*hid);
//        if(level_set_id[w.halfedge()] != 0 &&
//           valency(m, w.vertex())==3 &&
//           valency(m, w.opp().vertex())==3 &&
//           ((level_set_id[w.next().halfedge()] == 0 &&level_set_id[w.opp().next().halfedge()] == 0) ||
//           (level_set_id[w.prev().halfedge()] == 0 &&level_set_id[w.opp().prev().halfedge()] == 0)) &&
//           precond_collapse_edge(m, w.halfedge()))
//        {
//            cout << "collapsing!!!" << endl;
//            m.collapse_edge(w.halfedge(), true);
//            did_work = true;
//        }
//    }
//    }while (did_work);


//    for(VertexIDIterator vid = m.vertices_begin(); vid != m.vertices_end(); ++vid)
//    {
//        bool is_max = true;
//        bool is_min = true;
//        Walker w = m.walker(*vid);
//        float f = fun[*vid];
//        for(;!w.full_circle(); w = w.circulate_vertex_ccw())
//        {
//            if(fun[w.vertex()] < f)
//                is_min = false;
//            if(fun[w.vertex()] > f)
//                is_max = false;
//            
//        }
//        if(is_max || is_min)
//            status[*vid] = 2;
//    }



