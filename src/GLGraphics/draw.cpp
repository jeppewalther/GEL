/* ----------------------------------------------------------------------- *
 * This file is part of GEL, http://www.imm.dtu.dk/GEL
 * Copyright (C) the authors and DTU Informatics
 * For license and list of authors, see ../../doc/intro.pdf
 * ----------------------------------------------------------------------- */

#include "../GL/glew.h"

#include "../CGLA/Mat4x4f.h"
#include "../CGLA/Vec3d.h"
#include "../HMesh/Manifold.h"

#include "draw.h"
#include "SinglePassWireframeRenderer.h"
#include "IDBufferWireFrameRenderer.h"
#include "SOIL.h"

namespace GLGraphics
{
    using namespace CGLA;
    using namespace HMesh;
    using namespace std;
    
    namespace
    {
        void set_material(const Geometry::Material& material)
        {
            if(material.has_texture && material.tex_id >=0)
            {
                glEnable(GL_TEXTURE_2D);
                glBindTexture(GL_TEXTURE_2D, material.tex_id);
                glTexEnvi(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_MODULATE);
            }
            else
                glDisable(GL_TEXTURE_2D);
            
            glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT, material.ambient);
            glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, material.diffuse);
            glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, material.specular);
            glMaterialf(GL_FRONT_AND_BACK, GL_SHININESS, material.shininess);
        }
    }
    
    void draw(const Geometry::IndexedFaceSet& geometry)
    {
        glBegin(GL_TRIANGLES);
        for(int i=0;i<geometry.no_faces();i++)
        {
            Vec3i g_face = geometry.face(i);
            Vec3f vert0 = geometry.vertex(g_face[0]);
            Vec3f vert1 = geometry.vertex(g_face[1]);
            Vec3f vert2 = geometry.vertex(g_face[2]);
            Vec3f norm = normalize(cross(vert1-vert0, vert2-vert0));
            glNormal3fv(norm.get());
            glVertex3fv(vert0.get());
            glVertex3fv(vert1.get());
            glVertex3fv(vert2.get());
        }
        glEnd();
    }
    
    void draw(const Geometry::TriMesh& tm, bool per_vertex_norms)
    {
        int old_mat_idx = -1;
        glBegin(GL_TRIANGLES);
        for(int i=0;i<tm.geometry.no_faces();i++)
        {
            int new_mat_idx = i<static_cast<int>(tm.mat_idx.size()) ? tm.mat_idx[i] : -1;
            if(new_mat_idx != old_mat_idx)
            {
                glEnd();
                set_material(tm.materials[tm.mat_idx[i]]);
                glBegin(GL_TRIANGLES);
                old_mat_idx = new_mat_idx;
            }
            Vec3i n_face = tm.normals.face(i);
            Vec3i g_face = tm.geometry.face(i);
            Vec3i t_face = tm.texcoords.face(i);
            
            if(!per_vertex_norms)
            {
                Vec3f vert0 = tm.geometry.vertex(g_face[0]);
                Vec3f vert1 = tm.geometry.vertex(g_face[1]);
                Vec3f vert2 = tm.geometry.vertex(g_face[2]);
                Vec3f norm = normalize(cross(vert1-vert0, vert2-vert0));
                glNormal3fv(norm.get());
            }
            for(int j=0;j<3;j++)
            {
                if(per_vertex_norms && n_face != Geometry::NULL_FACE)
                {
                    Vec3f norm = tm.normals.vertex(n_face[j]);
                    glNormal3fv(norm.get());
                }
                if(t_face != Geometry::NULL_FACE)
                {
                    Vec3f texc = tm.texcoords.vertex(t_face[j]);
                    glTexCoord2fv(texc.get());
                }
                Vec3f vert = tm.geometry.vertex(g_face[j]);
                glVertex3fv(vert.get());
            }
        }
        glEnd();
        glDisable(GL_TEXTURE_2D);
    }
    
    void load_textures(Geometry::TriMesh& tm)
    {
        for(unsigned int i=0;i<tm.materials.size(); ++i)
        {
            Geometry::Material& mat = tm.materials[i];
            if(mat.tex_name != "")
            {
                string name = mat.tex_path + mat.tex_name;
                mat.tex_id = SOIL_load_OGL_texture(name.data(), 0, 0, SOIL_FLAG_TEXTURE_REPEATS | SOIL_FLAG_INVERT_Y | SOIL_FLAG_POWER_OF_TWO);
            }
        }
    }

	/// Draw an object of type T which contains only triangles as wireframe. In practice T = Manifold or TriMesh.
	template<typename T>
    void draw_triangles_in_wireframe(T& m, bool per_vertex_norms, const CGLA::Vec3f& line_color)
	{
		static SinglePassWireframeRenderer swr;
		swr.enable(line_color);
		draw(m, per_vertex_norms);
		swr.disable();
	}
    
	template
    void draw_triangles_in_wireframe(HMesh::Manifold& m, bool per_vertex_norms, const CGLA::Vec3f& line_color);
    
	template
    void draw_triangles_in_wireframe(Geometry::TriMesh& m, bool per_vertex_norms, const CGLA::Vec3f& line_color);
    
    template<class T>
    void draw_wireframe_oldfashioned(const T& m, bool per_vertex_norms, const Vec3f& line_color)
    {
        // Store state that we change
        glPushAttrib(GL_POLYGON_BIT);
        GLboolean lights_on;
        glGetBooleanv(GL_LIGHTING, &lights_on);
        Vec4f current_color;
        glGetFloatv(GL_CURRENT_COLOR, &current_color[0]);
        
        // Draw filled
        glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
        draw(m, per_vertex_norms);
        
        // Draw lines
        glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
        glDisable(GL_LIGHTING);
        glEnable(GL_POLYGON_OFFSET_LINE);
        glPolygonOffset(0,-5);
        glColor3fv(line_color.get());
        draw(m, per_vertex_norms);
        
        // Put back old state
        glColor3fv(current_color.get());
        if(lights_on) glEnable(GL_LIGHTING);
        glPopAttrib();
    }
    template
    void draw_wireframe_oldfashioned(const HMesh::Manifold& m, bool per_vertex_norms, const Vec3f& line_color);
    
    template
    void draw_wireframe_oldfashioned(const Geometry::TriMesh& m, bool per_vertex_norms, const Vec3f& line_color);
    
    
    void draw(const Manifold& m, bool per_vertex_norms)
    {
        for(FaceIDIterator f = m.faces_begin(); f != m.faces_end(); ++f){
            if(!per_vertex_norms)
                glNormal3dv(normal(m, *f).get());
            if(no_edges(m, *f)== 3)
                glBegin(GL_TRIANGLES);
            else
                glBegin(GL_POLYGON);
            
            for(Walker w = m.walker(*f); !w.full_circle(); w = w.circulate_face_ccw()){
                Vec3d n = normal(m, w.vertex());
                if(per_vertex_norms)
                    glNormal3dv(n.get());
                glVertex3dv(m.pos(w.vertex()).get());
            }
            glEnd();
        }
    }
    
    
    
    void draw(const Geometry::AABox& box)
    {
        glBegin(GL_QUADS);
        Vec3f norm_neg[] = {Vec3f(0,0,-1), Vec3f(-1,0,0), Vec3f(0,-1,0)};
        Vec3f norm_pos[] = {Vec3f(0,0, 1), Vec3f( 1,0,0), Vec3f(0, 1,0)};
        for(int j=0;j<3;++j)
        {
            glNormal3fv(norm_neg[j].get());
            Vec3f p = box.get_pmin();
            glVertex3f(p[0], p[1], p[2]);
            p[(j+1)%3] = box.get_pmax()[(j+1)%3];
            glVertex3f(p[0], p[1], p[2]);
            p[j] = box.get_pmax()[j];
            glVertex3f(p[0], p[1], p[2]);
            p[(j+1)%3] = box.get_pmin()[(j+1)%3];
            glVertex3f(p[0], p[1], p[2]);
        }
        glEnd();
        glBegin(GL_QUADS);
        for(int j=0;j<3;++j)
        {
            glNormal3fv(norm_pos[j].get());
            Vec3f p = box.get_pmax();
            glVertex3f(p[0], p[1], p[2]);
            p[j] = box.get_pmin()[j];
            glVertex3f(p[0], p[1], p[2]);
            p[(j+1)%3] = box.get_pmin()[(j+1)%3];
            glVertex3f(p[0], p[1], p[2]);
            p[j] = box.get_pmax()[j];
            glVertex3f(p[0], p[1], p[2]);
        }
        glEnd();
    }
    
    void draw(const Geometry::OBox& box)
    {
        Mat4x4f m = identity_Mat4x4f();
        copy_matrix(box.get_rotation(), m);
        glPushMatrix();
        glMultMatrixf(m.get());
        draw(box.get_aabox());
        glPopMatrix();
    }
    
    /** Draw the tree. The first argument is the level counter, the second
     argument is the level at which to stop drawing. */
    template <class BoxType>
    void draw(const Geometry::BoundingINode<BoxType>& node, int level, int max_level)
    {
        if(level == max_level)
        {
            draw(node);
            return;
        }
        node->left->draw(level + 1, max_level);
        node->right->draw(level + 1, max_level);
    }
    
    template <class BoxType>
    void draw(const Geometry::BoundingLNode<BoxType>& node, int level, int max_level)
    {
#if USE_LEAF_BOXES
        draw(node); 
#endif
    }
    
    template <class BoxType>
    void draw(const Geometry::BoundingTree<BoxType>& tree, int max_level)
    {
        draw(*tree.root, 0, max_level);
    }
}
