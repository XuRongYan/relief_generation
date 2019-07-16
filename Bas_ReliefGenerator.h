//
// Created by 徐溶延 on 2019-07-16.
//

#ifndef RELIEF_GENERATION_BAS_RELIEFGENERATOR_H

#include "MeshIo.h"
#include "MeshAlgorithm.h"

#define RELIEF_GENERATION_BAS_RELIEFGENERATOR_H
using namespace std;
using namespace pmp;

class Bas_ReliefGenerator {
public:
    int generate(SurfaceMesh& mesh, const string &modelPath, const string &resultPath, Scalar);


private:
    int compression(SurfaceMesh &, Scalar);
};


#endif //RELIEF_GENERATION_BAS_RELIEFGENERATOR_H
