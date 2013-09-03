/*
 *  VisObj.cpp
 *  GEL
 *
 *  Created by J. Andreas Bærentzen on 20/09/08.
 *  Copyright 2008 __MyCompanyName__. All rights reserved.
 *
 */

#include "VisObj.h"
#include "polarize.h"

#include <GLGraphics/Console.h>
#include <HMesh/Manifold.h>
#include <HMesh/AttributeVector.h>
#include <HMesh/load.h>
#include <HMesh/curvature.h>

#include <CGLA/Mat3x3d.h>
#include <CGLA/Vec3d.h>
#include <CGLA/Vec4d.h>

using namespace std;
using namespace CGLA;
using namespace HMesh;
using namespace GLGraphics;

int WINX=800, WINY=800;

void VisObj::refit()
{
	bsphere(mani, bsphere_center, bsphere_radius);
    
	view_ctrl.set_centre(Vec3f(bsphere_center));
	view_ctrl.set_eye_dist(2*bsphere_radius);
}

bool VisObj::reload(string _file)
{
	if(_file != "") file = _file;
	mani.clear();
	if(!load(file, mani))
		return false;
    refit();
	return true;
}

bool VisObj::add_mesh(string file)
{
	if(!load(file, mani))
		return false;
	bsphere(mani, bsphere_center, bsphere_radius);
	view_ctrl.set_centre(Vec3f(bsphere_center));
	view_ctrl.set_eye_dist(2*bsphere_radius);
	return true;
}



void VisObj::display(const std::string& display_method , Console& cs, bool smooth, float gamma)
{
	if(create_display_list){
		create_display_list = false;
        
		delete renderer;
        static Console::variable<float> r(-1);
        r.reg(cs,"display.bsphere_radius","Radius of the bounding sphere of present object");
        if(r<=0)
            r = bsphere_radius;
        
		string short_name = display_method.substr(0,3);
		if(short_name== "wir")
			renderer = new WireframeRenderer(mani, smooth);
        
		else if(short_name == "har")
			renderer = new HarmonicsRenderer(harmonics);
        
		else if(short_name == "iso")
			renderer = new IsophoteLineRenderer(mani, smooth);
        
		else if(short_name == "ref")
			renderer = new ReflectionLineRenderer(mani, smooth);
        
		else if(short_name == "gla")
			renderer = new GlazedRenderer(mani, smooth, bsphere_radius);
        
		else if(short_name == "too")
			renderer = new ToonRenderer(mani, smooth);
        
		else if(short_name == "cur"){
            static Console::variable<string> line_direction("min");
            static Console::variable<string> method("tensors");
            static Console::variable<int> smoothing_iter(1);
            
            line_direction.reg(cs,"display.curvature_lines.direction", "");
            method.reg(cs, "display.curvature_lines.method", "");
            smoothing_iter.reg(cs, "display.curvature_lines.smoothing_iter", "");
			
			VertexAttributeVector<Mat3x3d> curvature_tensors(mani.allocated_vertices());
			VertexAttributeVector<Vec3d> min_curv_direction(mani.allocated_vertices());
			VertexAttributeVector<Vec3d> max_curv_direction(mani.allocated_vertices());
			string _line_direction = line_direction;
			VertexAttributeVector<Vec3d>& lines = (_line_direction == "min") ? min_curv_direction : max_curv_direction;
			VertexAttributeVector<double> curvature(mani.allocated_vertices());
			
            if(string(method) == "tensors")
            {
                curvature_tensors_from_edges(mani, curvature_tensors);
                for(int i=0;i<smoothing_iter; ++i)
                    smooth_curvature_tensors(mani,curvature_tensors);
                
                curvature_from_tensors(mani, curvature_tensors,
                                       min_curv_direction,
                                       max_curv_direction,
                                       curvature);
            }
            else
                curvature_paraboloids(mani,
                                      min_curv_direction,
                                      max_curv_direction,
                                      curvature);
			renderer = new LineFieldRenderer(mani, smooth, lines, r);
		}
		else if(short_name == "gau"){
            static Console::variable<float> smoothing(2.0f);
            smoothing.reg(cs, "display.gaussian_curvature_renderer.smoothing", "");
			VertexAttributeVector<double> scalars(mani.allocated_vertices());
			gaussian_curvature_angle_defects(mani, scalars, smoothing);
			double max_G = 0;
            
            for(VertexIDIterator v = mani.vertices_begin(); v != mani.vertices_end(); ++v)
				max_G = max(abs(scalars[*v]), max_G);
            
			renderer = new ScalarFieldRenderer(mani, smooth, scalars, max_G, gamma);
			
		}
		else if(short_name == "mea"){
            static Console::variable<int> smoothing(2);
            smoothing.reg(cs, "display.mean_curvature_renderer.smoothing", "");
            
			VertexAttributeVector<double> scalars(mani.allocated_vertices());
			mean_curvatures(mani, scalars, smoothing);
			double max_G = 0;
			double mean = 0;
            
            for(VertexIDIterator v = mani.vertices_begin(); v != mani.vertices_end(); ++v){
				max_G = max(abs(scalars[*v]), max_G);
				mean += scalars[*v];
			}
            
			renderer = new ScalarFieldRenderer(mani, smooth, scalars, max_G, gamma);
		}
		else if(short_name == "amb"){
            static Console::variable<int> smoothing(2);
            smoothing.reg(cs, "display.ambient_occlusion_renderer.smoothing", "");
            
			VertexAttributeVector<double> scalars(mani.allocated_vertices());
			mean_curvatures(mani, scalars, smoothing);
			double max_G = 0;
            
            for(VertexIDIterator v = mani.vertices_begin(); v != mani.vertices_end(); ++v)
				max_G = max(abs(scalars[*v]), max_G);
            
			renderer = new AmbientOcclusionRenderer(mani, smooth, scalars, max_G);
			
		}		else
			renderer = new NormalRenderer(mani, smooth);
		
	}
	view_ctrl.set_gl_modelview();
	renderer->draw();
}
