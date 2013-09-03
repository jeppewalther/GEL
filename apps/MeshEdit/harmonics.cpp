/*
 *  harmonics.cpp
 *  GEL
 *
 *  Created by J. Andreas Bærentzen on 01/09/08.
 *  Copyright 2008 __MyCompanyName__. All rights reserved.
 *
 */

#include <iostream>

#include "harmonics.h"

#if USE_SPARSE_MATRIX
#include <arlgsym.h>
#endif

#include <CGLA/Vec3d.h>
#include <CGLA/Mat2x2d.h>
#include <LinAlg/Matrix.h>
#include <LinAlg/Vector.h>
#include <LinAlg/LapackFunc.h>

#include <GL/glew.h>
#include <GLGraphics/glsl_shader.h>

#include <HMesh/mesh_optimization.h>
#include <HMesh/curvature.h>
#include <HMesh/triangulate.h>
#include <HMesh/load.h>
#include <HMesh/x3d_save.h>

#include <GLGraphics/Console.h>

//#include "CSCMatrixBuilder.h"

using namespace CGLA;
using namespace std;
using namespace HMesh;
using namespace Geometry;
using namespace GLGraphics;
using namespace LinAlg;

namespace
{
	
	string vss =
	"#version 120\n"
	"#extension GL_EXT_gpu_shader4 : enable\n"
	"	\n"
	"	\n"
	"	attribute float eigenvalue;\n"
	"	attribute float eigenvalue2;\n"
	"	varying vec3 normal;\n"
	"	varying float eig;\n"
	"	varying float eig2;\n"
	"	\n"
	"	void main(void)\n"
	"	{\n"
	"		gl_Position =  ftransform();\n"
	"		normal = normalize(gl_NormalMatrix * gl_Normal);\n"
	"		eig = eigenvalue;\n"
	"		eig2 = eigenvalue2;\n"
	"	}\n";
	
	string fss = 	
	"#version 120\n"
	"#extension GL_EXT_gpu_shader4 : enable\n"
	"	\n"
	"	varying vec3 normal;\n"
	"	varying float eig;\n"
	"	varying float eig2;\n"
	"	uniform float eig_max;\n"
	"	uniform float eig_max2;\n"
	"	uniform bool do_highlight;\n"
	"	uniform bool do_diffuse;\n"
	"	const vec3 light_dir = vec3(0,0,1);\n"
	"	\n"
	" float basef(float x) {return max(0.0,min(1.0,2.0-4.0*abs(x)));\n}" 
	"	void main()\n"
	"	{\n"
	"		float dot_ln = max(0.0,dot(light_dir, normal));\n"
	"		\n"
	"		float eig_norm = eig/eig_max;\n"
	"		float stripe_signal = 250 * eig_norm;\n"
	//"		vec4 stripe_col = abs(stripe_signal) < 3.14 ? vec4(1,1,0,0) : vec4(.4,.4,.4,0);\n"
	"		vec4 stripe_col = -vec4(.4,.4,.4,0);\n"
	"		\n"
	"       float alpha = (1.0-eig_norm) * 2.0 * 3.1415926;\n"
	"       float offs = 2.0*3.1415/3.0;\n"
	"		gl_FragColor = vec4(0,0,1,0)*basef(eig_norm)+vec4(0,1,0,0)*basef(eig_norm-0.5)+vec4(1,0,0,0)* basef(eig_norm-1.0);\n"
	"		if(do_diffuse)   gl_FragColor *= dot_ln;\n"
	"		if(do_highlight) gl_FragColor += dot_ln*dot_ln*dot_ln*dot_ln*dot_ln*dot_ln*dot_ln*vec4(.5,.5,.5,0);\n"
	"		gl_FragColor -= vec4(.4,.4,.4,.4)*smoothstep(0.2,0.6,cos(stripe_signal));\n"
	"	}\n";
}



#if USE_SPARSE_MATRIX

double calculate_laplace_entry(VertexCirculator vc,Vec3d vertex)
{
	//* this piece of code is taken from the orginal make_laplace_operator...
	HalfEdgeIter h = vc.get_halfedge();
	Vec3d nbr(h->vert->pos);
	Vec3d left(h->next->vert->pos);
	Vec3d right(h->opp->prev->opp->vert->pos);
	
	double d_left = dot(normalize(nbr-left),normalize(vertex-left));
	double d_right = dot(normalize(nbr-right),normalize(vertex-right));
	double a_left  = acos(min(1.0, max(-1.0, d_left))); // w positive...?
	double a_right = acos(min(1.0, max(-1.0, d_right)));
	
	double w = 1.0/tan(a_left) + 1.0/tan(a_right);
	return -w;
}

void Harmonics::make_laplace_operator_sparse()
{
	bool is_generalized=true;
	CSCMatrixBuilder<double> mb_M;
	//+GENERALIZED
	CSCMatrixBuilder<double> mb_S;
	// S has the same type as V
	S.Resize(mani.no_vertices());
	
	for(VertexIter v = mani.vhandles_begin(); v != mani.vhandles_end(); ++v)
		if(!is_boundary(v))
		{
			int i = v->touched;
			double area_i = voronoi_area(v);
			Vec3d vertex(v->pos);
			double a_sum = 0;
			for(VertexCirculator vc(v); !vc.end(); ++vc)
			{
				int j = vc.get_vertex()->touched;
				double entry = calculate_laplace_entry(vc,vertex);
				if(!is_generalized)
				{
					double area_j = voronoi_area(vc.get_vertex());
					entry /= sqrt(area_i*area_j);
				}
				if(j > i)
					mb_M.insert_entry(j,entry);
				a_sum += entry;
			}
			//cout << a_sum << " ";
			mb_M.insert_entry(i,-a_sum);
			mb_M.next_column();
			
			if(!is_generalized)
				area_i = 1; // if standard S is an identity matrix;
			mb_S.insert_entry(i,area_i);
			mb_S.next_column_nonsort();
			S[i] = area_i;
		}     
	
	cout << "solving generalized problem... i = " << mani.no_vertices()<< endl;
	
	
	//+STANDARD
	//ArpackPP::ARluSymStdEig<double> dprob(50L, mb_M.get_Matrix(), "SA");
	//+GENERALIZED (shifted inv mode)
	ARluSymGenEig<double> dprob('S',maximum_eigenvalue, mb_M.get_Matrix(), mb_S.get_Matrix(), 0.0);
	//ArpackPP::ARluSymGenEig<double> dprob(number_of_eigenvalues, mb_M.get_Matrix(), mb_S.get_Matrix());
	
	dprob.FindEigenvectors();
	int conv = dprob.ConvergedEigenvalues();
	cout << conv << " eigenvectors found" << endl;
	Q.Resize(conv, dprob.GetN());
	V.Resize(conv);
	
	qnorm.Resize(conv);
	for(int i = 0; i < conv; i++)
	{
		V[i] = dprob.Eigenvalue(i);
		qnorm[i] = 0;
		for(int j = 0; j < dprob.GetN(); j++)
		{
			Q[i][j] = dprob.Eigenvector(i,j);
			qnorm[i] += Q[i][j]*Q[i][j];
		}
		
		cout << "(" << i << ":" << sqrt(V[i]/V[1]
										) << ")" << V[i]/V[1] << "V= "<< V[i] << " exp " << exp(-0.01*V[i]/V[1]) << " Qnorm " << qnorm[i]<< endl;
	}
	
	cout  << endl;
}

Harmonics::Harmonics(Manifold& _mani):mani(_mani)
{
	assert(is_initialized);
	
	maximum_eigenvalue = min(size_t(500),mani.no_vertices());
	
	triangulate(mani);
	mani.enumerate_vertices();
	make_laplace_operator_sparse();
	
	if(maximum_eigenvalue == -1)
	{
		cout << "not found" << endl;
		return;
	}
	
	proj.resize(maximum_eigenvalue);
	
	max_eig_values.resize(maximum_eigenvalue, 1e-10f);
	
	cout << endl << "Proj" << endl;
	for(int es=0; es<maximum_eigenvalue; ++es)  //o(n^2)
	{
		proj[es] = Vec3d(0.0);
		for(VertexIter v=mani.vhandles_begin(); v != mani.vhandles_end(); ++v)
		{
			proj[es] +=  Vec3d(v->pos) * Q[es][v->touched] * S[v->touched];
			max_eig_values[es] = max(max_eig_values[es], static_cast<float>(abs(Q[es][v->touched])));     
		}     
	}
	cout << endl;
}

#else
void Harmonics::make_laplace_operator()
{
	Q.Resize(mani.no_vertices(), mani.no_vertices());
	
	for(VertexIDIterator v = mani.vertices_begin(); v != mani.vertices_end(); ++v)
		if(!boundary(mani, *v)){
			int i = vtouched[*v];
			double area_i = voronoi_area(mani, *v);
			Vec3d vertex(mani.pos(*v));
			Vec3d curv_normal(0);
			double a_sum = 0;
            for(Walker wv = mani.walker(*v); !wv.full_circle(); wv = wv.circulate_vertex_cw())
			{
				int j = vtouched[wv.vertex()];
				double area_j = voronoi_area(mani, wv.vertex());
				
				Vec3d nbr(mani.pos(wv.vertex()));
				Vec3d left(mani.pos(wv.next().vertex()));
				Vec3d right(mani.pos(wv.opp().prev().opp().vertex()));
				
				double d_left = dot(normalize(nbr-left),
									normalize(vertex-left));
				double d_right = dot(normalize(nbr-right),
									 normalize(vertex-right));
				double a_left  = acos(min(1.0, max(-1.0, d_left)));
				double a_right = acos(min(1.0, max(-1.0, d_right)));
				
				double w = 1.0/tan(a_left) + 1.0/tan(a_right);
				
				Q[i][j] = -w/sqrt(area_i*area_j);						
				//Q[i][j] = -1;						
				a_sum += Q[i][j];
			}
			Q[i][i] = -a_sum;
		}
	EigenSolutionsSym(Q,V);
}

Harmonics::Harmonics(HMesh::Manifold& _mani, GLGraphics::Console& cs):mani(_mani), vtouched(_mani.allocated_vertices(), 0)
{
    string shader_path = "/Users/jab/GEL/apps/MeshEdit/";
	GLuint vs = create_glsl_shader(GL_VERTEX_SHADER, vss);
	GLuint fs = create_glsl_shader(GL_FRAGMENT_SHADER, fss);
	
	// Create the program
	prog_P0 = glCreateProgram();
	
	// Attach all shaders
	if(vs) glAttachShader(prog_P0, vs);
	if(fs) glAttachShader(prog_P0, fs);
	
	// Link the program object and print out the info log
	glLinkProgram(prog_P0);
	print_glsl_program_log(prog_P0);
	
	// Install program object as part of current state
	glUseProgram(0);
	
    display_harmonics_diffuse.reg(cs, "display.harmonics.diffuse", "");
    display_harmonics_time.reg(cs, "display.harmonics.time", "");
    display_harmonics_highlight.reg(cs, "display.harmonics.highlight", "");
	
	triangulate_by_edge_face_split(mani);

    int i = 0;
    for(VertexIDIterator v = mani.vertices_begin(); v != mani.vertices_end(); ++v, ++i)
        vtouched[*v] = i;
	maximum_eigenvalue = mani.no_vertices()-1;
	make_laplace_operator();
	
	proj.resize(maximum_eigenvalue);
	max_eig_values.resize(maximum_eigenvalue, 1e-10f);
	
	cout << "Projection magnitude" << endl;
	for(int es=0; es<maximum_eigenvalue; ++es)
	{
		proj[es] = Vec3d(0.0);
		for(VertexIDIterator v = mani.vertices_begin(); v != mani.vertices_end(); ++v)
		{
			proj[es] +=  Vec3d(mani.pos(*v)) * Q[es][vtouched[*v]];
			max_eig_values[es] = max(max_eig_values[es], static_cast<float>(abs(Q[es][vtouched[*v]])));
		}
	}
}

#endif

void Harmonics::add_frequency(int es, float scale)
{
	if(es<maximum_eigenvalue)
		for(VertexIDIterator v = mani.vertices_begin(); v != mani.vertices_end(); ++v){
			Vec3d p = Vec3d(proj[es]);
			double Qval = Q[es][vtouched[*v]];
			
			mani.pos(*v) += p * Qval * scale; 	
		}
}

void Harmonics::reset_shape()
{
	for(VertexIDIterator v = mani.vertices_begin(); v != mani.vertices_end(); ++v)
		mani.pos(*v) = Vec3d(0);	
}
void Harmonics::partial_reconstruct(int E0, int E1, float scale)
{
	for(int es=E0;es<=E1;++es)
		add_frequency(es, scale);
}

void Harmonics::parse_key(unsigned char key)
{
	switch(key) {
		case '+': 
			display_harmonics_time = display_harmonics_time+0.001; 
			break;
		case '-': 
			display_harmonics_time = display_harmonics_time-0.001; 
			break;
		case 'd':	
			display_harmonics_diffuse = !display_harmonics_diffuse; 
			break;
		case 'h':
			display_harmonics_highlight = !display_harmonics_highlight;
			break;			
	}
	
}


void Harmonics::draw_adf()
{
	vector<double> F(mani.no_vertices(),0);
	double F_max = 0;
	for(VertexIDIterator v = mani.vertices_begin(); v != mani.vertices_end(); ++v){
		int v_idx = vtouched[*v];
		for(int e = 1; e < V.Length(); ++e)
			F[v_idx] += sqr(Q[e][v_idx]) * exp(-(display_harmonics_time)*V[e]/V[1]);
		F_max = max(F[v_idx], F_max);
	}
	cout << "F max" <<  F_max << endl;
	glUseProgram(prog_P0);
	glUniform1f(glGetUniformLocation(prog_P0,"eig_max"),F_max);//2*M_PI);
	glUniform1i(glGetUniformLocation(prog_P0,"do_diffuse"),display_harmonics_diffuse);
   	glUniform1i(glGetUniformLocation(prog_P0,"do_highlight"),display_harmonics_highlight);
	GLuint attrib = glGetAttribLocationARB(prog_P0, "eigenvalue");
	
	glFrontFace(GL_CW);
	for(FaceIDIterator f = mani.faces_begin(); f != mani.faces_end(); ++f){
		glBegin(GL_TRIANGLES);
        for(Walker w = mani.walker(*f); !w.full_circle(); w = w.circulate_face_cw()){
			int i = vtouched[w.vertex()];
			glVertexAttrib1f(attrib,F[i]);
			glNormal3dv(normal(mani, w.vertex()).get());
			glVertex3dv(mani.pos(w.vertex()).get());
		}
		glEnd();
	}
	glFrontFace(GL_CCW);
	glUseProgram(0);
}

