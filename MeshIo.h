//
// Created by 徐溶延 on 2019-07-16.
//

#ifndef RELIEF_GENERATION_MESHIO_H

#include <pmp/visualization/MeshViewer.h>
#include <pmp/SurfaceMeshIO.h>
#include <string>
#include <iostream>
#define RELIEF_GENERATION_MESHIO_H
using namespace pmp;
using namespace std;

class MeshIo {
public:
    int readMesh(SurfaceMesh&, const string&);

    int saveMesh(const SurfaceMesh &mesh, string path);

    int showMesh(const string&);
};


#endif //RELIEF_GENERATION_MESHIO_H
