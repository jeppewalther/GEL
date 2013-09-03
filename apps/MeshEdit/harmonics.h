/*
 *  harmonics.h
 *  GEL
 *
 *  Created by J. Andreas Bærentzen on 01/09/08.
 *  Copyright 2008 __MyCompanyName__. All rights reserved.
 *
 */

#ifndef __MESHEDIT_HARMONICS_H__
#define __MESHEDIT_HARMONICS_H__

#include <CGLA/Vec3d.h>
#include <HMesh/Manifold.h>
#include <HMesh/AttributeVector.h>
#include <LinAlg/Matrix.h>
#include <LinAlg/Vector.h>
#include <GLGraphics/ManifoldRenderer.h>
#include <GLGraphics/Console.h>

#define USE_SPARSE_MATRIX 0

class Harmonics
{
	HMesh::Manifold& mani;
    HMesh::VertexAttributeVector<int> vtouched;
	
	int maximum_eigenvalue;
    
    bool is_initialized;
    GLuint prog_P0;
	
	std::vector<CGLA::Vec3d> proj;
	
	LinAlg::CMatrix Q;
	LinAlg::CVector qnorm;
	LinAlg::CVector V;
	LinAlg::CVector S;
    
    GLGraphics::Console::variable<float> display_harmonics_time;
    GLGraphics::Console::variable<int> display_harmonics_diffuse;
    GLGraphics::Console::variable<int> display_harmonics_highlight;
    
	std::vector<float> max_eig_values;
    
	void make_laplace_operator();
	void make_laplace_operator_sparse();
    
public:
	
	/// Initial analysis of harmonics
	Harmonics(HMesh::Manifold& mani, GLGraphics::Console& cs);
	
	/// Add a frequency to mesh reconstruction
	void add_frequency(int f, float scale = 1.0f);
	
	/// Reset the shape to use 0 eigenvalues
	void reset_shape();
	
	/// Do a partial reconstruct with an interval of eigenvalues
	void partial_reconstruct(int E0, int E1, float scale=1.0f);
	
	/// Parse keystrokes that would influence the interactive display
	void parse_key(unsigned char key);
	
	/// Draw with eigenvalues
	void draw_adf();
	
};


class HarmonicsRenderer: public GLGraphics::ManifoldRenderer
{
    
public:
    HarmonicsRenderer(Harmonics* h)
    {
        glNewList(display_list,GL_COMPILE);
        if(h) h->draw_adf();
        glEndList();
    }
};


#endif

