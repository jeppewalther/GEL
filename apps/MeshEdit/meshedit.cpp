/*
 *  MeshEdit is a small application which allows you to load and edit a mesh.
 *  The mesh will be stored in GEL's half edge based Manifold data structure.
 *  A number of editing operations are supported. Most of these are accessible from the 
 *  console that pops up when you hit 'esc'.
 *
 *  Created by J. Andreas Bærentzen on 15/08/08.
 *  Copyright 2008 __MyCompanyName__. All rights reserved.
 *
 */

#include <string>
#include <iostream>
#include <vector>
#include <algorithm>

#include <GL/glew.h>

#include <GLGraphics/Console.h>


#include <CGLA/eigensolution.h>
#include <CGLA/Vec2d.h>
#include <CGLA/Vec3d.h>
#include <CGLA/Mat3x3d.h>
#include <CGLA/Mat2x2d.h>
#include <CGLA/Mat2x3d.h>
#include <CGLA/Mat4x4d.h>

#include <LinAlg/Matrix.h>
#include <LinAlg/Vector.h>
#include <LinAlg/LapackFunc.h>

#include <GLGraphics/gel_glut.h>

#include <HMesh/Manifold.h>
#include <HMesh/AttributeVector.h>
#include <HMesh/mesh_optimization.h>
#include <HMesh/curvature.h>
#include <HMesh/triangulate.h>
#include <HMesh/flatten.h>
#include <HMesh/dual.h>
#include <HMesh/load.h>
#include <HMesh/quadric_simplify.h>
#include <HMesh/smooth.h>
#include <HMesh/x3d_save.h>
#include <HMesh/obj_save.h>
#include <HMesh/off_save.h>
#include <HMesh/mesh_optimization.h>
#include <HMesh/triangulate.h>
#include <HMesh/cleanup.h>
#include <HMesh/cleanup.h>
#include <HMesh/refine_edges.h>
#include <HMesh/subdivision.h>

#include <Util/Timer.h>
#include <Util/ArgExtracter.h>

#include "polarize.h"
#include "harmonics.h"
#include "VisObj.h"

using namespace std;
using namespace HMesh;
using namespace Geometry;
using namespace GLGraphics;
using namespace CGLA;
using namespace Util;
using namespace LinAlg;

// Single global instance so glut can get access
Console theConsole;
bool console_visible = false;


inline VisObj& get_vis_obj(int i)
{
    static VisObj vo[9];
    return vo[i];
}

Console::variable<int> active(0);

inline VisObj& avo()
{
    return get_vis_obj(active);
}

inline Manifold& active_mesh()
{
    return avo().mesh();
}

inline GLViewController& active_view_control()
{
    return avo().view_control();
}



////////////////////////////////////////////////////////////////////////////////
bool MyConsoleHelp(const std::vector<std::string> & args)
{
    theConsole.printf("");
    theConsole.printf("----------------- HELP -----------------");
    theConsole.printf("Press ESC key to open and close console");
    theConsole.printf("Press TAB to see the available commands and functions");
    theConsole.printf("Functions are shown in green and variables in yellow");
    theConsole.printf("Setting a value: [command] = value");
    theConsole.printf("Getting a value: [command]");
    theConsole.printf("Functions: [function] [arg1] [arg2] ...");
    theConsole.printf("Entering arg1=? or arg1=help will give a description.");
    theConsole.printf("History: Up and Down arrow keys move through history.");
    theConsole.printf("Tab Completion: TAB does tab completion and makes suggestions.");
    theConsole.printf("");
    theConsole.printf("Keyboard commands (when console is not active):");
    theConsole.printf("w   : switch to display.render_mode = wireframe");
    theConsole.printf("i   : switch to display.render_mode = isophotes");
    theConsole.printf("r   : switch to display.render_mode = reflection");
    theConsole.printf("m   : switch to display.render_mode = metallic");
    theConsole.printf("g   : switch to display.render_mode = glazed");
    theConsole.printf("n   : switch to display.render_mode = normal");
    theConsole.printf("h   : switch to display.render_mode = harmonics");
    theConsole.printf("f   : toggle smooth/flat shading");
    theConsole.printf("1-9 : switch between active meshes.");
    theConsole.printf("d   : (display.render_mode = harmonics) diffuse light on and off");
    theConsole.printf("h   : (display.render_mode = harmonics) highlight on and off ");
    theConsole.printf("+/- : (display.render_mode = harmonics) which eigenvector to show");
    theConsole.printf("q   : quit program");
    theConsole.printf("ESC : open console");
    theConsole.printf("");
    theConsole.printf("Mouse: Left button rotates, middle zooms, right pans");
    theConsole.printf("----------------- HELP -----------------");
    theConsole.printf("");
    return true;
}

bool wantshelp(const std::vector<std::string> & args)
{
    if(args.size() == 0) 
        return false;
	
    string str = args[0];
	
    if(str=="help" || str=="HELP" || str=="Help" || str=="?") 
        return true;
	
    return false;
}

/// Function that aligns two meshes.
void console_align(const std::vector<std::string> & args)
{
    if(wantshelp(args)) {
        theConsole.printf("usage: align <dest> <src>");
        theConsole.printf("This function aligns dest mesh with src");
        theConsole.printf("In practice the GLViewController of src is copied to dst.");
        theConsole.printf("both arguments are mandatory and must be numbers between 1 and 9.");
        theConsole.printf("Note that results might be unexpexted if the meshes are not on the same scale");
    }
	
    int dest = 0;
	
    if(args.size()>0){
        istringstream a0(args[0]);
        a0 >> dest;
        --dest;
		
        if(dest <0 || dest>8)
        {
            theConsole.printf("dest mesh out of range (1-9)");
            return;
        }
    }
    else
    {
        theConsole.printf("neither source nor destination mesh?!");
        return;
    }
	
    int src = 0;
    if(args.size()>1){
        istringstream a1(args[1]);
        a1 >> src;
        --src;
		
        if(src <0 || src>8)
        {
            theConsole.printf("src mesh out of range (1-9)");
            return;
        }
    }
    else
    {
        theConsole.printf("no src mesh?");
        return;
    }
    get_vis_obj(dest).view_control() = get_vis_obj(src).view_control();
}

void console_polarize(const std::vector<std::string> & args)
{
    if(wantshelp(args)) {
        theConsole.printf("usage: polarize");
        return;
    }
    int divisions = 10;
	
    if(args.size() > 0){
        istringstream a0(args[0]);
        a0 >> divisions;
    }

    avo().save_old();

	double vmin, vmax;
    VertexAttributeVector<double> fun;
    VertexAttributeVector<double> par;
    make_height_fun(active_mesh(), fun, vmin, vmax);
    polarize_mesh(active_mesh(), fun, vmin, vmax, divisions, par);
}

void transform_mesh(Manifold& mani, const Mat4x4d& m)
{
    for(VertexIDIterator vid = mani.vertices_begin(); vid != mani.vertices_end(); ++vid)
        mani.pos(*vid) = m.mul_3D_point(mani.pos(*vid));
}

void console_scale(const std::vector<std::string> & args)
{
    if(wantshelp(args)) {
        theConsole.printf("usage: scale sx sy sz");
        return;
    }
    
    Vec3d s;

    if(args.size() > 0){
        istringstream a0(args[0]);
        a0 >> s[0];
    }
    if(args.size() > 1){
        istringstream a0(args[0]);
        a0 >> s[1];
    }
    if(args.size() > 2){
        istringstream a0(args[0]);
        a0 >> s[2];
    }
    
    avo().save_old();
    transform_mesh(avo().mesh(),scaling_Mat4x4d(s));
    avo().refit();
}


void console_flatten(const std::vector<std::string> & args)
{
    if(wantshelp(args)) {
        theConsole.printf("usage: flatten <floater|harmonic|barycentric>");
        theConsole.printf("This function flattens a meshs with a simple boundary. It is mostly for showing mesh");
        theConsole.printf("parametrization methods. The current mesh MUST have a SINGLE boundary loop");
        theConsole.printf("This loop is mapped to the unit circle in a regular fashion (equal angle intervals).");
        theConsole.printf("All non boundary vertices are placed at the origin. Then the system is relaxed iteratively");
        theConsole.printf("using the weight scheme given as argument.");
        return;
    }
	
    avo().save_old();
	
    WeightScheme ws = BARYCENTRIC_W;
    if(args.size()>0){
        if(args[0] == "floater")
            ws = FLOATER_W;
        else if(args[0] == "harmonic")
            ws = HARMONIC_W;
        else if(args[0] == "lscm")
            ws = LSCM_W;
    }
    else
        return;
	
    flatten(active_mesh(), ws);
	
    return;
}

void console_save(const std::vector<std::string> & args)
{
    if(wantshelp(args)) {
        theConsole.printf("usage: save <name.x3d|name.obj> ");
		
        return;
    }
    const string& file_name = args[0];
    if(args.size() == 1){
        if(file_name.substr(file_name.length()-4,file_name.length())==".obj"){
            obj_save(file_name, active_mesh());
			
            return;
        }
        else if(file_name.substr(file_name.length()-4,file_name.length())==".off"){
            off_save(file_name, active_mesh());
			
            return;
        }
        else if(file_name.substr(file_name.length()-4,file_name.length())==".x3d"){
            x3d_save(file_name, active_mesh());
			
            return;
        }
        theConsole.printf("unknown format");
        return; 
    }
    theConsole.printf("usage: save <name.x3d|name.obj> ");
}


void console_refine_edges(const std::vector<std::string> & args)
{
    if(wantshelp(args)) {
        theConsole.printf("usage: refine.split_edges <length>");
        theConsole.printf("splits edges longer than <length>; default is 0.5 times average length");
        return;
    }
	
    avo().save_old();
	
    float thresh = 0.5f;
	
    if(args.size() > 0){
        istringstream a0(args[0]);
        a0 >> thresh;
    }
	
    float avg_length = average_edge_length(active_mesh());
	
    refine_edges(active_mesh(), thresh * avg_length);
	
    return;
	
}

void console_refine_faces(const std::vector<std::string> & args)
{
    if(wantshelp(args)) {
        theConsole.printf("usage: refine.split_faces ");
        theConsole.printf("usage:  Takes no arguments. Inserts a vertex at the centre of each face.");
		
        return;
    }
    avo().save_old();
	
    triangulate_by_vertex_face_split(active_mesh());
	
    return;
	
}

void console_cc_subdivide(const std::vector<std::string> & args)
{
    if(wantshelp(args)) {
        theConsole.printf("usage: refine.catmull_clark ");
        theConsole.printf("Does one step of Catmull-Clark subdivision");
		
        return;
    }
    avo().save_old();
	
    cc_split(active_mesh(),active_mesh());
    cc_smooth(active_mesh());
	
    return;
}

void console_doosabin_subdivide(const std::vector<std::string> & args)
{
    if(wantshelp(args)) {
        theConsole.printf("usage: refine.doo_sabin ");
        theConsole.printf("Does one step of Doo-Sabin Subdivision");
		
        return;
    }
    avo().save_old();
	
    cc_split(active_mesh(),active_mesh());
    dual(active_mesh());
    
    return;
}

void console_dual(const std::vector<std::string> & args)
{
    if(wantshelp(args)) 
    {
        theConsole.printf("usage: dual ");
        theConsole.printf("Produces the dual by converting each face to a vertex placed at the barycenter.");
        return;
    }
    avo().save_old();
	
    dual(active_mesh());
	
    return;
}


void console_minimize_curvature(const std::vector<std::string> & args)
{
    if(wantshelp(args)) 
    {
        theConsole.printf("usage: optimize.minimize_curvature <anneal>");
        theConsole.printf("Flip edges to minimize mean curvature.");
        theConsole.printf("If anneal is true, simulated annealing (slow) is used rather than a greedy scheme");
        return;
    }
    avo().save_old();
	
    bool anneal=false;
    if(args.size() > 0)
    {
        istringstream a0(args[0]);
        a0 >> anneal;
    }
	
    minimize_curvature(active_mesh(), anneal);
    avo().post_create_display_list();
    return;
}

void console_minimize_dihedral(const std::vector<std::string> & args)
{
    if(wantshelp(args))
    {
        theConsole.printf("usage: optimize.minimize_dihedral <iter> <anneal> <use_alpha> <gamma> ");
        theConsole.printf("Flip edges to minimize dihedral angles.");
        theConsole.printf("Iter is the max number of iterations. anneal tells us whether to use ");
        theConsole.printf("simulated annealing and not greedy optimization. use_alpha (default=true) ");
        theConsole.printf("means to use angle and not cosine of anglegamma (default=4) is the power ");
        theConsole.printf("to which we raise the dihedral angle");
        return;
    }
    avo().save_old();
	
    int iter = 1000;
    if(args.size() > 0)
    {
        istringstream a0(args[0]);
        a0 >> iter;
    }
	
    bool anneal = false;
    if(args.size() > 1)
    {
        istringstream a0(args[1]);
        a0 >> anneal;
    }
	
    bool use_alpha = true;
    if(args.size() > 2)
    {
        istringstream a0(args[2]);
        a0 >> use_alpha;
    }
	
    float gamma = 4.0f;
    if(args.size() > 3)
    {
        istringstream a0(args[3]);
        a0 >> gamma;
    }
	
	
    minimize_dihedral_angle(active_mesh(), iter, anneal, use_alpha, gamma);
    return;
}

void console_maximize_min_angle(const std::vector<std::string> & args)
{
    if(wantshelp(args)) 
    {
        theConsole.printf("usage: optimize.maximize_min_angle <thresh> <anneal>");
        theConsole.printf("Flip edges to maximize min angle - to make mesh more Delaunay.");
        theConsole.printf("If the dot product of the normals between adjacent faces < thresh");
        theConsole.printf("no flip will be made. anneal selects simulated annealing rather ");
        theConsole.printf("nthan greedy optimization.");
        return;
    }
    avo().save_old();
	
    float thresh = 0.0f;
    if(args.size() > 0)
    {
        istringstream a0(args[0]);
        a0 >> thresh;
    }
    bool anneal = false;
    if(args.size() > 1)
    {
        istringstream a0(args[1]);
        a0 >> anneal;
    }
    maximize_min_angle(active_mesh(),thresh,anneal);
    return;
}


void console_optimize_valency(const std::vector<std::string> & args)
{
    if(wantshelp(args)) 
    {
        theConsole.printf("usage: optimize.valency <anneal> ");
        theConsole.printf("Optimizes valency for triangle meshes. Anneal selects simulated annealing rather than greedy optim.");
        return;
    }
    avo().save_old();
	
    bool anneal = false;
    if(args.size() > 0)
    {
        istringstream a0(args[0]);
        a0 >> anneal;
    }
    optimize_valency(active_mesh(), anneal);
    return;
}

void console_analyze(const std::vector<std::string> & args)
{
    if(wantshelp(args)) 
    {
        theConsole.printf("usage:  harmonics.analyze");
        theConsole.printf("Creates the Laplace Beltrami operator for the mesh and finds all eigensolutions.");
        theConsole.printf("It also projects the vertices onto the eigenvectors - thus transforming the mesh");
        theConsole.printf("to this basis.");
        theConsole.printf("Note that this will stall the computer for a large mesh - as long as we use Lapack.");
        return;
    }
    avo().harmonics_analyze_mesh(theConsole);
    return;
}


void console_partial_reconstruct(const std::vector<std::string> & args)
{
    if(args.size() != 3)
        theConsole.printf("usage: haramonics.partial_reconstruct <e0> <e1> <s>");
	
    if(wantshelp(args)) {
        theConsole.printf("Reconstruct from projections onto eigenvectors. The two first arguments indicate");
        theConsole.printf("the eigenvector interval that we reconstruct from. The last argument is the ");
        theConsole.printf("scaling factor. Thus, for a vertex, v, the formula for computing the position, p, is:");
        theConsole.printf("for (i=e0; i<=e1;++i) p += proj[i] * Q[i][v] * s;");
        theConsole.printf("where proj[i] is the 3D vector containing the x, y, and z projections of the mesh onto");
        theConsole.printf("eigenvector i. Q[i][v] is the v'th coordinate of the i'th eigenvector.");
        theConsole.printf("Note that if vertex coordinates are not first reset, the result is probably unexpected.");
    }
    avo().save_old();
	
    if(args.size() != 3)
        return;
	
    int E0,E1;
    float scale;
    istringstream a0(args[0]);
    a0 >> E0;
    istringstream a1(args[1]);
    a1 >> E1;
    istringstream a2(args[2]);
    a2 >> scale;
    avo().harmonics_partial_reconstruct(E0,E1,scale);
    return;
}

void console_reset_shape(const std::vector<std::string> & args)
{
    if(wantshelp(args)) 
    {
        theConsole.printf("usage: harmonics.reset_shape ");
        theConsole.printf("Simply sets all vertices to 0,0,0. Call this before doing partial_reconstruct");
        theConsole.printf("unless you know what you are doing.");
        return;
    }
    avo().save_old();
    avo().harmonics_reset_shape();
    return;
}


void console_close_holes(const std::vector<std::string> & args)
{
    if(wantshelp(args)) 
    {
        theConsole.printf("usage: cleanup.close_holes");
        theConsole.printf("This function closes holes. It simply follows the loop of halfvectors which");
        theConsole.printf("enclose the hole and add a face to which they all point.");
        return;
    }
    avo().save_old();
	
    close_holes(active_mesh());
    return;
}

void console_reload(const std::vector<std::string> & args)
{
    if(wantshelp(args)) 
    {
        theConsole.printf("usage:  load <file>");
        theConsole.printf("(Re)loads the current file if no argument is given, but");
        theConsole.printf("if an argument is given, then that becomes the current file");
        return;
    }
    avo().save_old();
	
    if(!avo().reload(args.size() > 0 ? args[0]:""))
        theConsole.printf("failed to load");
	
    return;
}


void console_add_mesh(const std::vector<std::string> & args)
{
    if(wantshelp(args)) 
    {
        theConsole.printf("usage:  add_mesh <file>");
        theConsole.printf("Loads the file but without clearing the mesh. Thus, the loaded mesh is added to the");
        theConsole.printf("current model.");
        return;
    }
    avo().save_old();
	
    if(!avo().add_mesh(args.size() > 0 ? args[0]:""))
        theConsole.printf("failed to load");
	
    return;
}
void console_valid(const std::vector<std::string> & args)
{
    if(wantshelp(args)) 
    {
        theConsole.printf("usage:  validity");
        theConsole.printf("Tests validity of Manifold");
        return;
    }
	if(valid(active_mesh()))
		theConsole.printf("Mesh is valid");
	else
		theConsole.printf("Mesh is invalid - check console output");
	return;
}

void console_info(const std::vector<std::string> & args)
{
    if(wantshelp(args)) 
    {
        theConsole.printf("usage:  info");
        theConsole.printf("Provides information about mesh.");
        return;
    }
    Vec3d p0, p7;
    bbox(active_mesh(), p0, p7);
    stringstream bbox_corners;
    bbox_corners << p0 << " - " << p7 << endl;
	theConsole.printf("Bounding box corners : %s", bbox_corners.str().c_str());
    map<int,int> val_hist;
    
    for(VertexIDIterator vi = active_mesh().vertices_begin(); vi != active_mesh().vertices_end(); ++vi)
    {
        int val = valency(active_mesh(), *vi);
        if(val_hist.find(val) == val_hist.end())
            val_hist[val] = 0;
        ++val_hist[val];
    }
    
    theConsole.printf("Valency histogam");
    for(map<int,int>::iterator iter = val_hist.begin(); iter != val_hist.end(); ++iter)
    {
        stringstream vhl;
        vhl << iter->first << ", " << iter->second;
        theConsole.printf("%d, %d", iter->first, iter->second);
    }

	theConsole.printf("Mesh contains %d faces", active_mesh().no_faces());
	theConsole.printf("Mesh contains %d halfedges", active_mesh().no_halfedges());
	theConsole.printf("Mesh contains %d vertices", active_mesh().no_vertices());
	return;
}


void console_simplify(const std::vector<std::string> & args)
{
    if(wantshelp(args)) 
    {
        theConsole.printf("usage: simplify <fraction> ");
        theConsole.printf("Performs Garland Heckbert (quadric based) mesh simplification.");
        theConsole.printf("The only argument is the fraction of vertices to keep.");
        return;
    }
    avo().save_old();
	
    float keep_fraction;
    if(args.size() == 0)
    {
        theConsole.print("you must specify fraction of vertices to keep");
        return;
    }
    istringstream a0(args[0]);
    a0 >> keep_fraction;
	
    Vec3d p0, p7;
    bbox(active_mesh(), p0, p7);
    Vec3d d = p7-p0;
    float s = 1.0/d.max_coord();
    Vec3d pcentre = (p7+p0)/2.0;
    for(VertexIDIterator vi = active_mesh().vertices_begin(); vi != active_mesh().vertices_end(); ++vi){
        active_mesh().pos(*vi) = (active_mesh().pos(*vi) - pcentre) * s;
    }
    cout << "Timing the Garland Heckbert (quadric based) mesh simplication..." << endl;
    Timer timer;
    timer.start();
	
    //simplify
    quadric_simplify(active_mesh(),keep_fraction,0.0001f,true);
	
    cout << "Simplification complete, process time: " << timer.get_secs() << " seconds" << endl;
	
    //clean up the mesh, a lot of edges were just collapsed 
    active_mesh().cleanup();
	
    for(VertexIDIterator vi = active_mesh().vertices_begin(); vi != active_mesh().vertices_end(); ++vi)
        active_mesh().pos(*vi) = active_mesh().pos(*vi)*d.max_coord() + pcentre;
    return;
}

void console_vertex_noise(const std::vector<std::string> & args)
{
    if(wantshelp(args)) 
    {
        theConsole.printf("usage: noise.perturb_vertices <amplitude>");
        theConsole.printf("adds a random vector to each vertex. A random vector in the unit cube is generated and");
        theConsole.printf("to ensure an isotropic distribution, vectors outside the unit ball are discarded.");
        theConsole.printf("The vector is multiplied by the average edge length and then by the amplitude specified.");
        theConsole.printf("If no amplitude is specified, the default (0.5) is used.");
        return;
    }
    avo().save_old();
	
    float avg_length = average_edge_length(active_mesh());
	
    float noise_amplitude = 0.5f;
    if(args.size() > 0) {
        istringstream a0(args[0]);
        a0 >> noise_amplitude;
    }
	
    gel_srand(0);
    for(VertexIDIterator vi = active_mesh().vertices_begin(); vi != active_mesh().vertices_end(); ++vi){
        Vec3d v;
        do{
            v = Vec3d(gel_rand(),gel_rand(),gel_rand());
            v /= (float)(GEL_RAND_MAX);
        } 
        while(sqr_length(v) > 1.0);
		
        v -= Vec3d(0.5);
        v *= 2.0;
        v *= noise_amplitude;
        v *= avg_length;
        active_mesh().pos(*vi) += v;
    }		
    return;
}

void console_perpendicular_vertex_noise(const std::vector<std::string> & args)
{
    if(wantshelp(args)) {
        theConsole.printf("usage: noise.perturb_vertices_perpendicular <amplitude>");
        theConsole.printf("adds the normal times a random scalar times amplitude times");
        theConsole.printf("times average edge length to the vertex. (default amplitude=0.5)");
        return;
    }
    avo().save_old();
	
    float avg_length = average_edge_length(active_mesh());
	
    float noise_amplitude = 0.5;
    if(args.size() > 0) 
    {
        istringstream a0(args[0]);
        a0 >> noise_amplitude;
    }
	
    VertexAttributeVector<Vec3d> normals(active_mesh().allocated_vertices());
    for(VertexIDIterator vi = active_mesh().vertices_begin(); vi != active_mesh().vertices_end(); ++vi)
        normals[*vi] = normal(active_mesh(), *vi);
	
    gel_srand(0);
    for(VertexIDIterator vi = active_mesh().vertices_begin(); vi != active_mesh().vertices_end(); ++vi)
    {
        float rval = 0.5-gel_rand() / float(GEL_RAND_MAX);
        active_mesh().pos(*vi) += normals[*vi]*rval*noise_amplitude*avg_length*2.0;
    }
    return;
}

void console_noisy_flips(const std::vector<std::string> & args)
{
    if(wantshelp(args)){
        theConsole.printf("usage:  noise.perturb_topology <iter>");
        theConsole.printf("Perform random flips. iter (default=1) is the number of iterations.");
        theConsole.printf("mostly for making nasty synthetic test cases.");
        return;
    }
    avo().save_old();
	
    int iter = 1;
    if(args.size() > 0){
        istringstream a0(args[0]);
        a0 >> iter;
    }
	
    randomize_mesh(active_mesh(),  iter);
    return;
}

void console_laplacian_smooth(const std::vector<std::string> & args)
{
    if(wantshelp(args)) {
        theConsole.printf("usage:  smooth.laplacian <weight> <iter>");
        theConsole.printf("Perform Laplacian smoothing. weight is the scaling factor for the Laplacian.");
        theConsole.printf("default weight = 1.0. Default number of iterations = 1");
        return;
    }
    avo().save_old();
	
    float t=1.0;
    if(args.size() > 0){
        istringstream a0(args[0]);
        a0 >> t;
    }
    int iter = 1;
    if(args.size()>1){
        istringstream a0(args[1]);
        a0 >> iter;
    }
    /// Simple laplacian smoothing with an optional weight.
    for(int i=0;i<iter;++i) 
        laplacian_smooth(active_mesh(), t);
    return;
}


void console_mean_curvature_smooth(const std::vector<std::string> & args){
    if(wantshelp(args)) {
        theConsole.printf("usage:  smooth.mean_curvature <weight> <iter>");
        theConsole.printf("Perform mean curvature smoothing. weight is the scaling factor for the");
        theConsole.printf("mean curvature vector which has been normalized by dividing by edge lengths");
        theConsole.printf("this allows for larger steps as suggested by Desbrun et al.");
        theConsole.printf("default weight = 1.0. Default number of iterations = 1");
        return;
    }
    avo().save_old();
	
    double t=1.0;
    if(args.size() > 0){
        istringstream a0(args[0]);
        a0 >> t;
    }
    int iter=1;
    if(args.size() > 1){
        istringstream a0(args[1]);
        a0 >> iter;
    }	
    VertexAttributeVector<Vec3d> new_pos(active_mesh().allocated_vertices());
    for(int j = 0; j < iter; ++j){
        for(VertexIDIterator v = active_mesh().vertices_begin(); v != active_mesh().vertices_end(); ++v) {
            Vec3d m;
            double w_sum;
            unnormalized_mean_curvature_normal(active_mesh(), *v, m, w_sum);
            new_pos[*v] = Vec3d(active_mesh().pos(*v))  + (t * m/w_sum);
        }
        for(VertexIDIterator v = active_mesh().vertices_begin(); v != active_mesh().vertices_end(); ++v)
            active_mesh().pos(*v) = new_pos[*v];
    }
    return;
}

void console_taubin_smooth(const std::vector<std::string> & args)
{
    if(wantshelp(args)){
        theConsole.printf("usage:  smooth.taubin <iter>");
        theConsole.printf("Perform Taubin smoothing. iter (default=1) is the number of iterations.");
        return;
    }
    avo().save_old();
	
    int iter = 1;
    if(args.size() > 0){
        istringstream a0(args[0]);
        a0 >> iter;
    }
    /// Taubin smoothing is similar to laplacian smoothing but reduces shrinkage
    taubin_smooth(active_mesh(),  iter);
	
    return;
}

void console_fvm_smooth(const std::vector<std::string> & args)
{	
    if(wantshelp(args)){
        theConsole.printf("usage: smooth.fuzzy_vector_median <iter>");
        theConsole.printf("Smooth normals using fuzzy vector median smoothing. iter (default=1) is the number of iterations");
        theConsole.printf("This function does a very good job of preserving sharp edges.");
        return;
    }
    avo().save_old();
	
    int iter=1;
    if(args.size() > 0){
        istringstream a0(args[0]);
        a0 >> iter;
    }
    // Fuzzy vector median smoothing is effective when it comes to preserving sharp edges. 
    fvm_smooth(active_mesh(),  iter);
	
    return;
}

void console_triangulate(const std::vector<std::string> & args)
{	
    if(wantshelp(args)) {
        theConsole.printf("usage:  triangulate");
        theConsole.printf("This function triangulates all non triangular faces of the mesh.");
        theConsole.printf("you may want to call it after hole closing. For a polygon it simply connects");
        theConsole.printf("the two closest vertices in a recursive manner until only triangles remain");
        return;
    }
    avo().save_old();
	
    shortest_edge_triangulate(active_mesh());
    active_mesh().cleanup();
	valid(active_mesh());
    return;
}

void console_remove_faces(const std::vector<std::string> & args)
{
    avo().save_old();
    
    gel_srand(0);

//    for (FaceIDIterator f= active_mesh().faces_begin(); f != active_mesh().faces_end(); ++f) {
//        if(gel_rand() < 0.5 * GEL_RAND_MAX)
//        {
//            active_mesh().remove_face(*f);
//        }
//    }

//    for (VertexIDIterator v= active_mesh().vertices_begin(); v != active_mesh().vertices_end(); ++v) {
//        if(gel_rand() < 0.005 * GEL_RAND_MAX)
//        {
//            active_mesh().remove_vertex(*v);
//        }
//    }
    for (HalfEdgeIDIterator h= active_mesh().halfedges_begin(); h != active_mesh().halfedges_end(); ++h) {
        if(gel_rand() < 0.005 * GEL_RAND_MAX)
        {
            active_mesh().remove_edge(*h);
        }
    }

    active_mesh().cleanup();
    valid(active_mesh());
	
    return;
}


void console_remove_caps(const std::vector<std::string> & args)
{	
    if(wantshelp(args)) {
        theConsole.printf("usage:  cleanup.remove_caps thresh");
        theConsole.printf("Remove caps (triangles with one very big angle). The thresh argument is the fraction of PI to");
        theConsole.printf("use as threshold for big angle. Default is 0.85. Caps are removed by flipping.");
        return;
    }
    avo().save_old();
	
    float t = 0.85f;
    if(args.size() > 0){
        istringstream a0(args[0]);
        a0 >> t;
    }
    remove_caps(active_mesh(), static_cast<float>(M_PI) *t);
    active_mesh().cleanup();
	
    return;
}

void console_remove_needles(const std::vector<std::string> & args)
{	
    if(wantshelp(args)){
        theConsole.printf("usage: cleanup.remove_needles <thresh>");
        theConsole.printf("Removes very short edges by collapse. thresh is multiplied by the average edge length");
        theConsole.printf("to get the length shorter than which we collapse. Default = 0.1");
        return;
    }
    avo().save_old();
	
    float thresh = 0.1f;
    if(args.size() > 0){
        istringstream a0(args[0]);
        a0 >> thresh;
    }
    float avg_length = average_edge_length(active_mesh());
    remove_needles(active_mesh(), thresh * avg_length);
    active_mesh().cleanup();
	
    return;
}

void console_undo(const std::vector<std::string> & args)
{	
    if(wantshelp(args)) {
        theConsole.printf("usage:  undo");
        theConsole.printf("This function undoes one operation. Repeated undo does nothing");
        return;
    }
    avo().restore_old();
    avo().refit();
    return;
}


void reshape(int W, int H)
{
    active_view_control().reshape(W,H);
}

Console::variable<string> display_render_mode("normal");
Console::variable<int> display_smooth_shading(true);
Console::variable<float> display_gamma(2.2);

void display()
{
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	

    glPushMatrix();

    avo().display(display_render_mode, theConsole, display_smooth_shading, display_gamma);
	
    glPopMatrix();
	
    if(console_visible)
    {
        glUseProgram(0);
        theConsole.display();
	}
    
    glutSwapBuffers();
}

void animate() 
{	
    //usleep( (int)1e4 );
    active_view_control().try_spin();
    glutPostRedisplay();
}


void mouse(int button, int state, int x, int y) 
{
    Vec2i pos(x,y);
    if (state==GLUT_DOWN) 
    {
        if (button==GLUT_LEFT_BUTTON && glutGetModifiers() == 0)
            active_view_control().grab_ball(ROTATE_ACTION,pos);
        else if (button==GLUT_MIDDLE_BUTTON || glutGetModifiers() == GLUT_ACTIVE_CTRL) 
            active_view_control().grab_ball(ZOOM_ACTION,pos);
        else if (button==GLUT_RIGHT_BUTTON || glutGetModifiers() == GLUT_ACTIVE_ALT)
            active_view_control().grab_ball(PAN_ACTION,pos);
    }
    else if (state==GLUT_UP)
        active_view_control().release_ball();
}

void motion(int x, int y) {
    Vec2i pos(x,y);
    active_view_control().roll_ball(Vec2i(x,y));
}


void keyboard_spec(int key, int x, int y)
{
    if (console_visible)
        theConsole.special(key);
    glutPostRedisplay();
}


void keyboard(unsigned char key, int x, int y) 
{
    //toggle console with ESC
    if (key == 27)
    {
        console_visible = !console_visible;
        glutPostRedisplay();
        return;
    }
    
    if (console_visible)
    {
        theConsole.keyboard(key);
        if(key == 13)
        {
            avo().post_create_display_list();
            glutPostRedisplay();
        }
        return;
    }
    else {		
		
        switch(key) {
			case 'q': exit(0);
			case '\033':
                console_visible = false;
				break;
			case '1':
			case '2':
			case '3':
			case '4':
			case '5':
			case '6':
			case '7':
			case '8':
			case '9':
				active = key - '1'; break;
			case 'f': display_smooth_shading = !display_smooth_shading; break;
			case 'w':
				display_render_mode = "wire"; break;
			case 'n':
				display_render_mode = "normal"; break;
			case 'i':
				display_render_mode = "isophotes"; break;
			case 'r':
				display_render_mode = "reflection"; break;
			case 'h':
				display_render_mode = "harmonics"; break;
			case 't':
				display_render_mode = "toon"; break;
			case 'g':
				display_render_mode = "glazed"; break;
			case 'a':
				display_render_mode = "ambient_occlusion"; break;
			case 'c':
				display_render_mode = "copper"; break;
			case 'C':
				display_render_mode = "curvature_lines"; break;
			case 'M':
				display_render_mode = "mean_curvature"; break;
			case 'G':
				display_render_mode = "gaussian_curvature"; break;
        }
		
        if(string(display_render_mode).substr(0,3) == "har")
            avo().harmonics_parse_key(key);
		
        if(key != '\033') avo().post_create_display_list();
    }
    
    glutPostRedisplay();
}

void init_glut(int argc, char** argv)
{  
    glutInitDisplayMode(GLUT_RGBA|GLUT_DOUBLE|GLUT_DEPTH|GLUT_ALPHA);
    glutInitWindowSize(WINX, WINY);
    glutInit(&argc, argv);
    glutCreateWindow("MeshEdit");
    glutDisplayFunc(display);
    glutKeyboardFunc(keyboard);
    glutSpecialFunc(keyboard_spec);
    glutReshapeFunc(reshape);
    glutMouseFunc(mouse);
    glutMotionFunc(motion);
    glutIdleFunc(animate);
}
void init_gl()
{
    glewInit();
    glEnable(GL_LIGHTING);
    glEnable(GL_LIGHT0);
    glLightModeli(GL_LIGHT_MODEL_TWO_SIDE, 1);
	
    // Set the value of a uniform
    //glUniform2f(glGetUniformLocation(prog_P0,"WIN_SCALE"), win_size_x/2.0, win_size_y/2.0);
	
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
    glClearColor(1,1,1, 0.f);
    glColor4f(1.0f, 1.0f, 1.0f, 0.f);
    float material[4] = {1,1,1,1};
    glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, material);
    glEnable(GL_DEPTH_TEST);
	
    theConsole.reg_cmdN("harmonics.reset_shape", console_reset_shape, "");
    theConsole.reg_cmdN("harmonics.analyze", console_analyze, "");
    theConsole.reg_cmdN("harmonics.partial_reconstruct", console_partial_reconstruct,"");
    theConsole.reg_cmdN("simplify", console_simplify,"");
    theConsole.reg_cmdN("smooth.mean_curvature", console_mean_curvature_smooth,"");
    theConsole.reg_cmdN("smooth.laplacian", console_laplacian_smooth,"");
    theConsole.reg_cmdN("smooth.taubin", console_taubin_smooth,"");
    theConsole.reg_cmdN("smooth.fuzzy_vector_median", console_fvm_smooth,"");
	
    theConsole.reg_cmdN("optimize.valency", console_optimize_valency,"");
    theConsole.reg_cmdN("optimize.minimize_dihedral_angles", console_minimize_dihedral,"");
    theConsole.reg_cmdN("optimize.minimize_curvature", console_minimize_curvature,"");
    theConsole.reg_cmdN("optimize.maximize_min_angle", console_maximize_min_angle,"");
    theConsole.reg_cmdN("cleanup.close_holes", console_close_holes,"");
    theConsole.reg_cmdN("load_mesh", console_reload,"");
    theConsole.reg_cmdN("add_mesh", console_add_mesh,"");
	
    theConsole.reg_cmdN("cleanup.remove_caps", console_remove_caps,"");
    theConsole.reg_cmdN("cleanup.remove_needles", console_remove_needles,"");
    theConsole.reg_cmdN("triangulate", console_triangulate,"");
    theConsole.reg_cmdN("refine.split_edges", console_refine_edges,"");
    theConsole.reg_cmdN("refine.split_faces", console_refine_faces,"");
    theConsole.reg_cmdN("refine.catmull_clark", console_cc_subdivide,"");
    theConsole.reg_cmdN("refine.doo_sabin", console_doosabin_subdivide,"");
    theConsole.reg_cmdN("save_mesh", console_save,"");
    theConsole.reg_cmdN("noise.perturb_vertices", console_vertex_noise,"");
    theConsole.reg_cmdN("noise.perturb_vertices_perpendicular", console_perpendicular_vertex_noise,"");
    theConsole.reg_cmdN("noise.perturb_topology", console_noisy_flips,"");

    theConsole.reg_cmdN("remove_faces", console_remove_faces,"");

    theConsole.reg_cmdN("dual", console_dual,"");
    theConsole.reg_cmdN("flatten", console_flatten,"");
	
    theConsole.reg_cmdN("align", console_align,"");
	
    theConsole.reg_cmdN("undo", console_undo,"");
	
	theConsole.reg_cmdN("validity", console_valid,"");
	theConsole.reg_cmdN("info", console_info,"");

    theConsole.reg_cmdN("polarize", console_polarize ,"");
    
    theConsole.reg_cmdN("transform.scale", console_scale, "Scale mesh");
    
    active.reg(theConsole, "active_mesh", "The active mesh");
    display_render_mode.reg(theConsole, "display.render_mode", "Display render mode");
    display_smooth_shading.reg(theConsole, "display.smooth_shading", "1 for smooth shading 0 for flat");
    display_gamma.reg(theConsole, "display.gamma", "The gamma setting for the display");

}

int main(int argc, char** argv)
{
    ArgExtracter ae(argc, argv);
	
    init_glut(argc, argv);
    init_gl();
	
    theConsole.print("Welcome to MeshEdit");
    theConsole.newline();
   
    if(argc>1){		
        vector<string> files;
		ae.get_all_args(files);
		for(size_t i=1;i<files.size();++i)
			get_vis_obj(i-1).reload(files[i]);
    }
    glutMainLoop();
    return 0;
}





