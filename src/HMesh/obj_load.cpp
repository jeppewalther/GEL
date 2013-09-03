/* ----------------------------------------------------------------------- *
 * This file is part of GEL, http://www.imm.dtu.dk/GEL
 * Copyright (C) the authors and DTU Informatics
 * For license and list of authors, see ../../doc/intro.pdf
 * ----------------------------------------------------------------------- */

#include <fstream>
#include "obj_load.h"

#include "../Geometry/obj_load.h"
#include "../Geometry/TriMesh.h"

#include "Manifold.h"

using namespace std;
using namespace CGLA;

namespace HMesh
{
    using std::string;
    using Geometry::TriMesh;
    using Geometry::obj_load;
    
    bool obj_load(const string& filename, Manifold& m)
    {
        ifstream ifs(filename.data());
        
        if(ifs)
        {
            vector<Vec3f> vertices;
            vector<int> faces;
            vector<int> indices;
            while(ifs.good() && !ifs.eof())
            {
                string tok;
                ifs >> tok;
                if(tok == "v")
                {
                    float x,y,z;
                    ifs >> x >> y >> z;
                    vertices.push_back(Vec3f(x,y,z));
                    char line[1000];
                    ifs.getline(line, 998);
                }
                else if(tok == "f")
                {
                    char line[1000];
                    ifs.getline(line, 998);
                    char* pch = strtok(line, " \t");
                    int ctr = 0;
                    while(pch != 0)
                    {
                        int v;
                        sscanf(pch, "%d", &v);
                        indices.push_back(v-1);
                        pch = strtok(0, " \t");
                        ++ctr;
                    }
                    if(ctr)
                        faces.push_back(ctr);
                }
                else
                {
                    char line[1000];
                    ifs.getline(line, 998);
                }
            }
            m.clear();
            m.build(vertices.size(),
                    reinterpret_cast<float*>(&vertices[0]),
                    faces.size(),
                    &faces[0],
                    &indices[0]);
            
            return true;
        }
        return false;
        
    }
}