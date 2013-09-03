/*
 *  VisObj.h
 *  GEL
 *
 *  Created by J. Andreas Bærentzen on 20/09/08.
 *  Copyright 2008 __MyCompanyName__. All rights reserved.
 *
 */
#ifndef __MESHEDIT_VISOBJ_H__
#define __MESHEDIT_VISOBJ_H__


#include <string>
#include <GL/glew.h>
#include <HMesh/Manifold.h>
#include <CGLA/Vec3d.h>
#include <GLGraphics/draw.h>
#include <GLGraphics/Console.h>
#include <GLGraphics/GLViewController.h>
#include <GLGraphics/ManifoldRenderer.h>
#include "harmonics.h"

extern int WINX;
extern int WINY;

class VisObj
	{
		std::string file;
		GLGraphics::GLViewController view_ctrl;
		bool create_display_list;
		HMesh::Manifold mani;
		HMesh::Manifold old_mani;
		
		Harmonics* harmonics;
        GLGraphics::ManifoldRenderer* renderer;
		CGLA::Vec3d bsphere_center;
		float bsphere_radius;
	public:
		VisObj():
        file(""), view_ctrl(WINX,WINY, CGLA::Vec3f(0), 1.0), create_display_list(true), harmonics(0) {}
		
		float get_bsphere_radius() const { return bsphere_radius;}
		
		HMesh::Manifold& mesh() {return mani;}
		
		void save_old() {old_mani = mani;}
		void restore_old() {mani = old_mani;}
		
		GLGraphics::GLViewController& view_control() {return view_ctrl;}
        
        void refit();
		
		bool reload(std::string _file);

		bool add_mesh(std::string _file);
		
		void display(const std::string& display_method , GLGraphics::Console& cs, bool smooth, float gamma);
		
		void post_create_display_list()
		{
			create_display_list = true;
		}
		
		void harmonics_analyze_mesh(GLGraphics::Console& cs)
		{
			delete harmonics;
			harmonics = new Harmonics(mani, cs);
		}
		
		void harmonics_reset_shape()
		{
			if(harmonics)
				harmonics->reset_shape();
		}
		
		void harmonics_parse_key(unsigned char key)
		{
			harmonics->parse_key(key);
		}
		
		void harmonics_partial_reconstruct(int eig0, int eig1, float scale)
		{
			if(harmonics)
				harmonics->partial_reconstruct(eig0, eig1, scale);
		}
		
	};

#endif
