//
// Created by 徐溶延 on 2019-07-16.
//

#ifndef RELIEF_GENERATION_MESHALGORITHM_H

#include <pmp/algorithms/NormalCone.h>
#include <pmp/algorithms/SurfaceNormals.h>
#include <pmp/algorithms/SurfaceRemeshing.h>
#include <pmp/algorithms/SurfaceSimplification.h>
#include <pmp/algorithms/DifferentialGeometry.h>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <iostream>

#define RELIEF_GENERATION_MESHALGORITHM_H
using namespace pmp;
using namespace std;

typedef Eigen::SparseMatrix<Scalar> SparseMatrixType;
typedef SurfaceMesh::FaceAroundVertexCirculator FaceNeighbors;
typedef SurfaceMesh::VertexAroundVertexCirculator VertexNeighbors;

class MeshAlgorithm {
public:
    Scalar remeshing(SurfaceMesh &);

    int simplification(SurfaceMesh &, Scalar);

    Eigen::Vector3f get_gradient_Phi(const SurfaceMesh&, const Face&, const Vertex &);

    Eigen::MatrixXf generate_b(SurfaceMesh&);

    SparseMatrixType generate_Laplace(const SurfaceMesh &mesh);

    int solve_poison(SurfaceMesh &, const SparseMatrixType &, const Eigen::MatrixXf &);

    Eigen::Vector3f getGradientField(const SurfaceMesh &mesh, const Face &f);
};


#endif //RELIEF_GENERATION_MESHALGORITHM_H
