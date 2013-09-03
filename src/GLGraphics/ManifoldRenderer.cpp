/* ----------------------------------------------------------------------- *
 * This file is part of GEL, http://www.imm.dtu.dk/GEL
 * Copyright (C) the authors and DTU Informatics
 * For license and list of authors, see ../../doc/intro.pdf
 * ----------------------------------------------------------------------- */

#include "ManifoldRenderer.h"

#include <algorithm>
#include <string>
#include <cstdlib>
#include "../CGLA/Mat3x3d.h"
#include "../GLGraphics/glsl_shader.h"
#include "../HMesh/Manifold.h"
#include "../HMesh/AttributeVector.h"
#include "../HMesh/curvature.h"

using namespace CGLA;
using namespace HMesh;
using namespace std;

namespace GLGraphics
{
    
    GLuint get_noise_texture_id()
    {
        static GLuint texname=0;
        static bool was_here = false;
        
        if(!was_here)
        {
            was_here = true;
            int width = 32;
            int height = 32;
            int depth = 32;
            vector<unsigned char> texels(width*height*depth);
            for (int i = 0; i < width*height*depth; ++i)
            {
                int intensity = 255.0 * (float(gel_rand()) / GEL_RAND_MAX);
                texels[i] = (unsigned char) intensity;
            }
            
            glGenTextures(1, &texname);	
            glBindTexture(GL_TEXTURE_3D, texname);	
            glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
            glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
            glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_WRAP_S, GL_REPEAT);
            glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_WRAP_T, GL_REPEAT);
            glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_WRAP_R, GL_REPEAT);
            glTexImage3D(GL_TEXTURE_3D, 0, GL_INTENSITY8, width, height, depth, 0, GL_RED, GL_UNSIGNED_BYTE, &texels[0]);
        }
        
        return texname;
    }
    
    
    int WireframeRenderer::maximum_face_valency(const Manifold& m)
    {
        int max_val = 0;
        for(FaceIDIterator f = m.faces_begin(); f != m.faces_end(); ++f)
            max_val = max(max_val, no_edges(m, *f));
        return max_val;
    }
    
    WireframeRenderer::WireframeRenderer(HMesh::Manifold& m, bool smooth): idbuff_renderer(0)
    {
        if(GLEW_EXT_geometry_shader4 && maximum_face_valency(m) > 3)
        {
            GLint viewp[4];
            glGetIntegerv(GL_VIEWPORT,viewp);
            idbuff_renderer = new IDBufferWireframeRenderer(viewp[2], viewp[3], m);
        }
        else
        {
            glNewList(display_list,GL_COMPILE);
            if(GLEW_EXT_geometry_shader4)
                draw_triangles_in_wireframe(m,smooth, Vec3f(1,0,0));				
            else
                draw_wireframe_oldfashioned(m,smooth, Vec3f(1,0,0));
            glEndList();
        }
    }
    
    void WireframeRenderer::draw()
    {
        if(idbuff_renderer)
        {
            glEnable(GL_LIGHTING);
            idbuff_renderer->draw(Vec3f(1,0,0),Vec3f(1));
            glDisable(GL_LIGHTING);
        }
        else
            glCallList(display_list);
    }
    
    void SimpleShaderRenderer::init_shaders(const std::string& vss, 
                                            const std::string& fss)
    {
        vs = create_glsl_shader(GL_VERTEX_SHADER, vss);
        print_glsl_program_log(vs);
        
        fs = create_glsl_shader(GL_FRAGMENT_SHADER, fss);
        print_glsl_program_log(fs);
        
        prog = glCreateProgram();
        
        if(vs) glAttachShader(prog, vs);
        if(fs) glAttachShader(prog, fs);
        
        glLinkProgram(prog);
        print_glsl_program_log(prog);
        
    }
    
    void SimpleShaderRenderer::compile_display_list(const Manifold& m, bool smooth)
    {
        glNewList(display_list,GL_COMPILE);
        GLGraphics::draw(m, smooth);
        glEndList();	
    }
    
    void SimpleShaderRenderer::draw()
    {
        GLint old_prog;
        glGetIntegerv(GL_CURRENT_PROGRAM, &old_prog);
        glUseProgram(prog);
        glCallList(display_list);
        glUseProgram(old_prog);
    }
    
    const string NormalRenderer::vss = 
    "varying vec3 _n;\n"
    "varying vec3 v;\n"
    "\n"
    "void main(void)\n"
    "{\n"
    "	gl_Position = ftransform();\n"
    "	v = vec3(gl_ModelViewMatrix * gl_Vertex);\n"
    "	_n = normalize(gl_NormalMatrix * gl_Normal);\n"
    "}\n";
    
    const string NormalRenderer::fss = 
    "varying vec3 _n;\n"
    "varying vec3 v;\n"
    "\n"
    "void main(void)\n"
    "{\n"
    "   vec3 n = normalize(_n);\n"
    "	vec3 l = normalize(-v);\n"
    "	vec3 e = l;\n"
    "	vec3 r = normalize(2.0*dot(l, n)*n - l);\n"
    "	\n"
    "	vec4 a = vec4(0.0,0.1,.3,1.0);\n"
    "   float dot_ln = abs(dot(l, n));\n"
    "	vec4 d = vec4(0.7) * dot_ln;\n"
    "	vec4 s = vec4(0.3)*smoothstep(0.98,0.9999,dot(r, e));\n"
    "	\n"
    "	gl_FragColor =  d+s;\n"
    "}\n";
    
    
    const string ReflectionLineRenderer::vss = 
    "varying vec3 _n;\n"
    "varying vec3 v;\n"
    "\n"
    "void main(void)\n"
    "{\n"
    "	gl_Position = ftransform();\n"
    "	v = vec3(gl_ModelViewMatrix * gl_Vertex);\n"
    "	_n = normalize(gl_NormalMatrix * gl_Normal);\n"
    "}\n";
    
    
    const string ReflectionLineRenderer::fss = 
    "uniform float detail;\n"
    "\n"
    "varying vec3 _n;\n"
    "varying vec3 v;\n"
    "\n"
    "void main(void)\n"
    "{\n"
    "   vec3 n = normalize(_n);\n"
    "	// calculate the reflection\n"
    "	vec3 r = normalize(2.0*dot(-v, n)*n + v);\n"
    "	vec3 viewer_lightdir = vec3(0, 0, 1.0);\n"
    "   float diff  = dot(n,viewer_lightdir);\n"
    "	\n"
    "	vec2 r2 = normalize(vec2(r[0], r[2]));\n"
    "	vec2 x = vec2(1, 0);\n"
    "	float angle = acos(dot(r2, x));\n"
    "	\n"
    "	// decide if we hit a white or black ring, based on y value\n"
    "	gl_FragColor = diff * vec4(1.0) + smoothstep(0.8, 1.0,cos(13.0*angle)) * vec4(-1.0);\n"
    "}\n";
    
    const string IsophoteLineRenderer::vss = 
    "varying vec3 _n;\n"
    "varying vec3 v;\n"
    "\n"
    "void main(void)\n"
    "{\n"
    "	gl_Position = ftransform();\n"
    "	v = vec3(gl_ModelViewMatrix * gl_Vertex);\n"
    "	_n = normalize(gl_NormalMatrix * gl_Normal);\n"
    "}\n";
    
    
    const string IsophoteLineRenderer::fss = 
    "uniform float detail;\n"
    "\n"
    "varying vec3 _n;\n"
    "varying vec3 v;\n"
    "\n"
    "void main(void)\n"
    "{\n"
    "   vec3 n = normalize(_n);\n"
    "	vec3 viewer_lightdir = vec3(0, 0, 1.0);\n"
    "	vec3 isophote_lightdir = viewer_lightdir;\n"
    "	float angle = acos(dot(n, isophote_lightdir));\n"
    "   float diff  = dot(n,viewer_lightdir);\n"
    "	\n"
    "	// decide if we hit a white or black ring, based on y value\n"
    "	gl_FragColor = diff * vec4(1.0) + smoothstep(0.8, 1.0,cos(20.0*angle)) * vec4(-1.0);\n"
    "}\n";
    
    const string ToonRenderer::vss = 
    "varying vec3 _n;\n"
    "varying vec3 v;\n"
    "\n"
    "void main(void)\n"
    "{\n"
    "	gl_Position = ftransform();\n"
    "	v = vec3(gl_ModelViewMatrix * gl_Vertex);\n"
    "	_n = normalize(gl_NormalMatrix * gl_Normal);\n"
    "}\n";
    
    const string ToonRenderer::fss = 
    "varying vec3 _n;\n"
    "varying vec3 v;\n"
    "\n"
    "void main(void)\n"
    "{\n"
    "   vec3 n = normalize(_n);\n"
    "	vec3 l = normalize(-v);\n"
    "	vec3 e = l;\n"
    "	vec3 r = normalize(2.0*dot(l, n)*n - l);\n"
    "	\n"
    "	vec4 a = vec4(0.0,0.1,.3,1.0);\n"
    "   float dot_ln = abs(dot(l, n));\n"
    "	vec4 d = vec4(0.7,0.7,0.0,1.0) * 0.25 * (smoothstep(0.23,0.25,dot_ln)+smoothstep(0.45,0.47,dot_ln)+smoothstep(0.7,0.72,dot_ln)+smoothstep(0.9,0.92,dot_ln));\n"
    "	vec4 s = vec4(0.5,0.3,0.4,1.0)*smoothstep(0.96,0.98,dot(r, e));\n"
    "	\n"
    "	gl_FragColor =  d+s;\n"
    "}\n";
    
    
    void GlazedRenderer::draw()
    {
        GLint old_prog;
        glGetIntegerv(GL_CURRENT_PROGRAM, &old_prog);
        glUseProgram(prog);
        glBindTexture(GL_TEXTURE_3D, get_noise_texture_id());
        glUniform1iARB(glGetUniformLocationARB(prog, "noise_tex"),0);
        glUniform1fARB(glGetUniformLocationARB(prog, "noise_scale"),12.0/bsphere_rad);
        glCallList(display_list);
        glUseProgram(old_prog);
    }
    
    
    
    const string GlazedRenderer::vss = 
    "varying vec3 _n;\n"
    "varying vec3 v;\n"
    "varying vec3 v_obj;\n"
    "\n"
    "void main(void)\n"
    "{\n"
    "	gl_Position = ftransform();\n"
    "   v_obj = gl_Vertex.xyz;\n"
    "	v = vec3(gl_ModelViewMatrix * gl_Vertex);\n"
    "	_n = normalize(gl_NormalMatrix * gl_Normal);\n"
    "}\n"
    "\n";
    
    const string GlazedRenderer::fss =
    "uniform sampler3D noise_tex;\n"
    "uniform float noise_scale;\n"
    "varying vec3 _n;\n"
    "varying vec3 v;\n"
    "varying vec3 v_obj;\n"
    "\n"
    "vec4 glazed_shader(vec4 mat_col,  vec4 light_col, vec3 light_dir)\n"
    "{\n"
    "   vec3 n = normalize(_n);\n"
    "	vec3 e = normalize(-v);\n"
    "	vec3 r = normalize(2.0*dot(e, n)*n - e);\n"
    "	float d = max(0.05,dot(light_dir, n));\n"
    "	vec4 diff = mat_col * light_col *d; 	\n"
    "	vec4 refl = smoothstep(0.7,0.75,dot(r,light_dir)) * light_col;\n"
    "	return 0.15*refl + diff;\n"
    "}\n"
    "\n"
    "void main(void)\n"
    "{\n"
    "	vec4 mat_col = vec4(0.9,1.0,0.4,1.0) + vec4(-0.1,-0.1,0.12,0.0) * texture3D(noise_tex, noise_scale*v_obj).x\n"
    " + vec4(0.05) * texture3D(noise_tex, 500.0*v_obj).x;\n"
    "	\n"
    "	vec3 light0_dir = vec3(0.0,1.0,0.0);\n"
    "	vec4 light0_col = vec4(0.7,0.9,1.0,1.0);\n"
    "	\n"
    "	vec3 light1_dir = vec3(0.0,0.0,1.0);\n"
    "	vec4 light1_col = vec4(1.0,1.0,0.7,1.0);\n"
    "	\n"
    "	gl_FragColor = \n"
    "	0.5*glazed_shader(mat_col, light0_col, light0_dir)+\n"
    "	0.5*glazed_shader(mat_col, light1_col, light1_dir);\n"
    "	\n"
    "	gl_FragColor.a = 1.0;\n"
    "}\n";
    
    
    const string ScalarFieldRenderer::vss =
    "	attribute float scalar;\n"
    "	varying vec3 _normal;\n"
    "	varying float s;\n"
    "	\n"
    "	void main(void)\n"
    "	{\n"
    "		gl_Position =  ftransform();\n"
    "		_normal = normalize(gl_NormalMatrix * gl_Normal);\n"
    "		s=scalar;\n"
    "	}\n";
    
    const string ScalarFieldRenderer::fss = 	
    "	varying vec3 _normal;\n"
    "	varying float s;\n"
    "	uniform float scalar_max;\n"
    "   uniform float gamma;\n"
    "	const vec3 light_dir = vec3(0,0,1);\n"
    "	\n"
    "	void main()\n"
    "	{\n"
    "       vec3 normal = normalize(_normal);\n"
    "		float dot_ln = max(0.0,dot(light_dir, normal));\n"
    "		\n"
    "		float s_norm = s/scalar_max;\n"
    "		float stripe_signal = 100.0 * s_norm;\n"
    "		vec4 stripe_col = abs(stripe_signal) < 3.14 ? vec4(0,0,0,0) : vec4(.1,.1,.1,0);\n"
    "		\n"
    "		gl_FragColor = s_norm * vec4(-1,0,1,0);\n"
    "       gl_FragColor *= dot_ln;\n"
    "       gl_FragColor.r = pow(gl_FragColor.r, 1.0/gamma);\n"
    "       gl_FragColor.g = pow(gl_FragColor.g, 1.0/gamma);\n"
    "       gl_FragColor.b = pow(gl_FragColor.b, 1.0/gamma);\n"
    "		gl_FragColor += stripe_col * smoothstep(0.8,1.0,cos(stripe_signal));\n"
    "	}\n";
    
    ScalarFieldRenderer::ScalarFieldRenderer(const Manifold& m, 
                                             bool smooth, 
                                             VertexAttributeVector<double>& field, 
                                             double max_val, float gamma): SimpleShaderRenderer(vss, fss)
    {
        
        GLint old_prog;
        glGetIntegerv(GL_CURRENT_PROGRAM, &old_prog);
        glUseProgram(prog);
        
        GLuint scalar_attrib = glGetAttribLocation(prog, "scalar");
        glUniform1fARB(glGetUniformLocationARB(prog, "scalar_max"), max_val);
        
        //    static float& gamma = CreateCVar("display.scalar_field_renderer.gamma",2.2f);
        glUniform1fARB(glGetUniformLocationARB(prog, "gamma"), gamma);
        glNewList(display_list,GL_COMPILE);
        
        for(FaceIDIterator f = m.faces_begin(); f != m.faces_end(); ++f){      
            if(!smooth) 
                glNormal3dv(normal(m, *f).get());
            if(no_edges(m, *f)== 3) 
                glBegin(GL_TRIANGLES);
            else 
                glBegin(GL_POLYGON);
            
            
            for(Walker w = m.walker(*f); !w.full_circle(); w = w.circulate_face_ccw()){
                Vec3d n(normal(m, w.vertex()));
                if(smooth) 
                    glNormal3dv(n.get());
                glVertexAttrib1d(scalar_attrib, field[w.vertex()]);
                glVertex3dv(m.pos(w.vertex()).get());
            }
            glEnd();
        }
        glEndList();	
        glUseProgram(old_prog);
        
    }

    const string PeriodicScalarFieldRenderer::vss =
    "	attribute float scalar;\n"
    "	varying vec3 _normal;\n"
    "	varying float s;\n"
    "	\n"
    "	void main(void)\n"
    "	{\n"
    "		gl_Position =  ftransform();\n"
    "		_normal = normalize(gl_NormalMatrix * gl_Normal);\n"
    "		s=scalar;\n"
    "	}\n";
    
    const string PeriodicScalarFieldRenderer::fss = 	
    "	varying vec3 _normal;\n"
    "	varying float s;\n"
    "   uniform float gamma;\n"
    "	const vec3 light_dir = vec3(0,0,1);\n"
    "	\n"
    "	void main()\n"
    "	{\n"
    "       vec3 normal = normalize(_normal);\n"
    "		float dot_ln = max(0.0,dot(light_dir, normal));\n"
    "		\n"
    "		gl_FragColor = sin(s) * vec4(-1,0,1,0);\n"
    "       gl_FragColor *= dot_ln;\n"
    "       gl_FragColor.r = pow(gl_FragColor.r, 1.0/gamma);\n"
    "       gl_FragColor.g = pow(gl_FragColor.g, 1.0/gamma);\n"
    "       gl_FragColor.b = pow(gl_FragColor.b, 1.0/gamma);\n"
    "	}\n";
    
    PeriodicScalarFieldRenderer::PeriodicScalarFieldRenderer(const Manifold& m, 
                                             bool smooth, 
                                             VertexAttributeVector<double>& field, 
                                             float gamma): SimpleShaderRenderer(vss, fss)
    {
        
        GLint old_prog;
        glGetIntegerv(GL_CURRENT_PROGRAM, &old_prog);
        glUseProgram(prog);
        
        GLuint scalar_attrib = glGetAttribLocation(prog, "scalar");
        
        //    static float& gamma = CreateCVar("display.scalar_field_renderer.gamma",2.2f);
        glUniform1fARB(glGetUniformLocationARB(prog, "gamma"), gamma);
        glNewList(display_list,GL_COMPILE);
        
        for(FaceIDIterator f = m.faces_begin(); f != m.faces_end(); ++f){      
            if(!smooth) 
                glNormal3dv(normal(m, *f).get());
            if(no_edges(m, *f)== 3) 
                glBegin(GL_TRIANGLES);
            else 
                glBegin(GL_POLYGON);
            
            
            for(Walker w = m.walker(*f); !w.full_circle(); w = w.circulate_face_ccw()){
                Vec3d n(normal(m, w.vertex()));
                if(smooth) 
                    glNormal3dv(n.get());
                glVertexAttrib1d(scalar_attrib, field[w.vertex()]);
                glVertex3dv(m.pos(w.vertex()).get());
            }
            glEnd();
        }
        glEndList();	
        glUseProgram(old_prog);
        
    }

    const string AmbientOcclusionRenderer::vss =
    "	attribute float scalar;\n"
    "	varying vec3 _normal;\n"
    "	varying float s;\n"
    "	\n"
    "	void main(void)\n"
    "	{\n"
    "		gl_Position =  ftransform();\n"
    "		_normal = normalize(gl_NormalMatrix * gl_Normal);\n"
    "		s=scalar;\n"
    "	}\n";
    
    const string AmbientOcclusionRenderer::fss = 	
    "	varying vec3 _normal;\n"
    "	varying float s;\n"
    "	uniform float scalar_max;\n"
    "	const vec3 light_dir = vec3(0,0,1);\n"
    "	\n"
    "	void main()\n"
    "	{\n"
    "   vec3 normal = normalize(_normal);\n"
    "		float dot_ln = max(0.0,dot(light_dir, normal));\n"
    "		\n"
    "		float s_norm = min(1.0,s/scalar_max+1.0);\n"
    "		\n"
    "		gl_FragColor = s_norm * vec4(1.0);\n"
    "       gl_FragColor *= dot_ln;\n"
    "       gl_FragColor.r = pow(gl_FragColor.r, 1.0);\n"
    "       gl_FragColor.g = pow(gl_FragColor.g, 1.0);\n"
    "       gl_FragColor.b = pow(gl_FragColor.b, 1.0);\n"
    "	}\n";
    
    AmbientOcclusionRenderer::AmbientOcclusionRenderer(const Manifold& m, bool smooth, VertexAttributeVector<double>& field, double max_val):
    SimpleShaderRenderer(vss,fss)
    {	
        GLint old_prog;
        glGetIntegerv(GL_CURRENT_PROGRAM, &old_prog);
        glUseProgram(prog);
        
        GLuint scalar_attrib = glGetAttribLocation(prog, "scalar");
        glUniform1fARB(glGetUniformLocationARB(prog, "scalar_max"), max_val);
        
        glNewList(display_list,GL_COMPILE);
        
        for(FaceIDIterator f = m.faces_begin(); f != m.faces_end(); ++f){
            
            if(!smooth) 
                glNormal3dv(normal(m, *f).get());
            if(no_edges(m, *f)== 3) 
                glBegin(GL_TRIANGLES);
            else 
                glBegin(GL_POLYGON);
            
            for(Walker w = m.walker(*f); !w.full_circle(); w = w.circulate_face_ccw())
            {
                Vec3d n(normal(m, w.vertex()));
                if(smooth) 
                    glNormal3dv(n.get());
                glVertexAttrib1d(scalar_attrib, field[w.vertex()]);
                glVertex3dv(m.pos(w.vertex()).get());
            }
            glEnd();
        }
        glEndList();	
        glUseProgram(old_prog);
        
    }
    
    
    LineFieldRenderer::LineFieldRenderer(const Manifold& m, bool smooth, VertexAttributeVector<Vec3d>& lines, float _r): 
    SimpleShaderRenderer(vss,fss), r(_r)
    {
        float noise_scale = 10.0f/r;
        float line_scale = 0.003f;
        
        GLint old_prog;
        glGetIntegerv(GL_CURRENT_PROGRAM, &old_prog);
        glUseProgram(prog);	
        glUniform1fARB(glGetUniformLocationARB(prog, "scale_line"),line_scale*noise_scale);
        glUniform1fARB(glGetUniformLocationARB(prog, "noise_scale"),noise_scale);
        glUniform1iARB(glGetUniformLocationARB(prog, "noise_tex"),0);
        GLuint direction = glGetAttribLocation(prog, "direction");	
        glNewList(display_list,GL_COMPILE);
        for(FaceIDIterator f = m.faces_begin(); f != m.faces_end(); ++f){
            if(!smooth) 
                glNormal3dv(normal(m, *f).get());
            if(no_edges(m, *f) == 3) 
                glBegin(GL_TRIANGLES);
            else 
                glBegin(GL_POLYGON);
            
            for(Walker w = m.walker(*f); !w.full_circle(); w = w.circulate_face_ccw()){
                Vec3d n(normal(m, w.vertex()));
                if(smooth) 
                    glNormal3dv(n.get());
                
                Vec3d d = lines[w.vertex()];
                d = normalize(d-n*dot(n,d));
                glVertexAttrib3dv(direction, d.get());
                glVertex3dv(m.pos(w.vertex()).get());
            }
            glEnd();
        }
        
        glBindTexture(GL_TEXTURE_3D, 0);
        glEndList();	
        glUseProgram(old_prog);
        
    }
    
    void LineFieldRenderer::draw()
    {
        GLint old_prog;
        glGetIntegerv(GL_CURRENT_PROGRAM, &old_prog);
        glUseProgram(prog);
        glBindTexture(GL_TEXTURE_3D, get_noise_texture_id());
        glCallList(display_list);
        glBindTexture(GL_TEXTURE_3D, 0);
        glUseProgram(old_prog);	
    }
    
    const string LineFieldRenderer::vss = 
    "attribute vec3 direction;\n"
    "varying vec3 _n;\n"
    "varying vec3 dir_obj;\n"
    "varying vec3 v_obj;\n"
    "\n"
    "void main(void)\n"
    "{\n"
    "	gl_Position = ftransform();\n"
    "   v_obj = gl_Vertex.xyz;\n"
    "	dir_obj = direction;\n"
    "	_n = normalize(gl_NormalMatrix * gl_Normal);\n"
    "}\n";
    
    const string LineFieldRenderer::fss =
    "uniform sampler3D noise_tex;\n"
    "uniform float scale_line;\n"
    "uniform float noise_scale;\n"
    "varying vec3 _n;\n"
    "varying vec3 dir_obj;\n"
    "varying vec3 v_obj;\n"
    "\n"
    "float tex(vec3 p) {return smoothstep(0.2,0.3,texture3D(noise_tex, p).x);}\n"
    "void main(void)\n"
    "{\n"
    "   vec3 n = normalize(_n);\n"
    "   float I = "
    "             tex(noise_scale*v_obj + 6.0*scale_line*dir_obj) + \n"
    "             tex(noise_scale*v_obj - 6.0*scale_line*dir_obj) + \n"
    "             tex(noise_scale*v_obj + 5.0*scale_line*dir_obj) + \n"
    "             tex(noise_scale*v_obj - 5.0*scale_line*dir_obj) + \n"
    "             tex(noise_scale*v_obj + 4.0*scale_line*dir_obj) + \n"
    "             tex(noise_scale*v_obj - 4.0*scale_line*dir_obj) + \n"
    "             tex(noise_scale*v_obj + 3.0*scale_line*dir_obj) + \n"
    "             tex(noise_scale*v_obj - 3.0*scale_line*dir_obj) + \n"
    "             tex(noise_scale*v_obj + 2.0*scale_line*dir_obj) + \n"
    "             tex(noise_scale*v_obj - 2.0*scale_line*dir_obj) + \n"
    "             tex(noise_scale*v_obj + 1.0*scale_line*dir_obj) + \n"
    "             tex(noise_scale*v_obj - 1.0*scale_line*dir_obj) + \n"
    "			  tex(noise_scale*v_obj); \n"
    "	\n"
    "   float diff = max(0.0,dot(n,vec3(0.0, 0.0, 1.0)));\n"
    "	gl_FragColor.rgb = vec3(diff*I*(1.0/13.0));\n"
    "	gl_FragColor.a = 1.0;\n"
    "}\n";
    
    const string DualVertexRenderer::vss = 
    "#version 120\n"
    "#extension GL_EXT_gpu_shader4 : enable\n"
    "varying vec4 diffuseIn;\n"
    "varying vec3 _normalIn;\n"
    "void main(void)\n"
    "{\n"
    "   diffuseIn = gl_Color;\n"
    "   _normalIn = normalize(gl_NormalMatrix*gl_Normal);\n"
    "   gl_Position =  ftransform();\n"
    "}\n";
    
    const string DualVertexRenderer::gss = 
    "#version 120\n"
    "#extension GL_EXT_gpu_shader4 : enable\n"
    "#extension GL_EXT_geometry_shader4 : enable\n"
    "\n"
    "varying in vec4 diffuseIn[3];\n"
    "varying in vec3 _normalIn[3];\n"
    "varying vec4 diffuse[3];\n"
    "varying float f;\n"
    "varying vec3 _normal;\n"
    "void main(void)\n"
    "{\n"
    "  diffuse[0] = diffuseIn[0];\n"
    "  diffuse[1] = diffuseIn[1];\n"
    "  diffuse[2] = diffuseIn[2];\n"
    "\n"
    "  f = diffuseIn[0].a;\n"
    "  gl_Position = gl_PositionIn[0];\n"
    "  _normal = _normalIn[0];\n"
    "  EmitVertex();\n"
    "	\n"
    "  f = diffuseIn[1].a;\n"
    "  gl_Position = gl_PositionIn[1];\n"
    "  _normal = _normalIn[1];\n"
    "  EmitVertex();\n"
    "\n"
    "  f = diffuseIn[2].a;\n"
    "  gl_Position = gl_PositionIn[2];\n"
    "  _normal = _normalIn[2];\n"
    "  EmitVertex();\n"
    "\n"
    "  EndPrimitive();\n"
    "}\n";
    
    const string DualVertexRenderer::fss =
    "#version 120\n"
    "#extension GL_EXT_gpu_shader4 : enable\n"
    "\n"
    "varying float f;\n"
    "varying vec4 diffuse[3];\n"
    "varying vec3 _normal;\n"
    "\n"
    "void main(void)\n"
    "{\n"
    "   vec3 normal = normalize(_normal);\n"
    "   float col_idx=0;\n"
    "   if(f>diffuse[0].g && f<diffuse[0].b)\n"
    "      col_idx = diffuse[0].r;\n"
    "   else if(f>diffuse[1].g && f<diffuse[1].b)\n"
    "      col_idx = diffuse[1].r;\n"
    "   else if(f>diffuse[2].g && f<diffuse[2].b)\n"
    "      col_idx = diffuse[2].r;\n"
    "   vec4 col = col_idx < .5 ? vec4(1,0,0,0) : vec4(0,0,1,0);\n"
    "\n"
    " 	gl_FragColor =col*dot(normal, vec3(0,0,1));\n"
    "}\n";
    
    
    
    
    DualVertexRenderer::DualVertexRenderer(const HMesh::Manifold& m, VertexAttributeVector<Vec4d>& field)
    {		
        // Create the program
        static GLuint prog = glCreateProgram();
        
        static bool was_here = false;
        if(!was_here)
        {
            was_here = true;
            // Create s	haders directly from file
            static GLuint vs = create_glsl_shader(GL_VERTEX_SHADER, vss);
            static GLuint gs = create_glsl_shader(GL_GEOMETRY_SHADER_EXT, gss);
            static GLuint fs = create_glsl_shader(GL_FRAGMENT_SHADER, fss);
            
            // Attach all shaders
            if(vs) glAttachShader(prog, vs);
            if(gs) glAttachShader(prog, gs);
            if(fs) glAttachShader(prog, fs);
            
            // Specify input and output for the geometry shader. Note that this must be
            // done before linking the program.
            glProgramParameteriEXT(prog,GL_GEOMETRY_INPUT_TYPE_EXT,GL_TRIANGLES);
            glProgramParameteriEXT(prog,GL_GEOMETRY_VERTICES_OUT_EXT,3);
            glProgramParameteriEXT(prog,GL_GEOMETRY_OUTPUT_TYPE_EXT,GL_TRIANGLE_STRIP);
            
            // Link the program object and print out the info log
            glLinkProgram(prog);
        }
        
        
        
        GLint old_prog;
        glGetIntegerv(GL_CURRENT_PROGRAM, &old_prog);
        
        glNewList(display_list,GL_COMPILE);
        glUseProgram(prog);
        for(FaceIDIterator f = m.faces_begin(); f != m.faces_end(); ++f){
            if(no_edges(m, *f) != 3) 
                continue;
            else 
                glBegin(GL_TRIANGLES);
            
            for(Walker w = m.walker(*f); !w.full_circle(); w = w.circulate_face_ccw()){
                Vec3d n(normal(m, w.vertex()));
                glNormal3dv(n.get());
                glColor4dv(field[w.vertex()].get());
                glVertex3dv(m.pos(w.vertex()).get());
            }
            glEnd();
        }
        glUseProgram(old_prog);
        glEndList();	
        
    }
    
}


