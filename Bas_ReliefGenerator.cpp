//
// Created by 徐溶延 on 2019-07-16.
//

#include "Bas_ReliefGenerator.h"

int Bas_ReliefGenerator::generate(SurfaceMesh& mesh, const string &modelPath, const string &resultPath, const Scalar beta) {
    MeshIo io;
    MeshAlgorithm algorithm;
    io.readMesh(mesh, modelPath);
    Scalar avg_scalar = algorithm.remeshing(mesh);
    algorithm.simplification(mesh, avg_scalar);
    SparseMatrixType Ls = algorithm.generate_Laplace(mesh);
    Eigen::MatrixXf b = algorithm.generate_b(mesh);
    algorithm.solve_poison(mesh, Ls, b);
    compression(mesh, beta);
    io.saveMesh(mesh, resultPath);
    io.showMesh(resultPath);
    return 0;
}

int Bas_ReliefGenerator::compression(SurfaceMesh &mesh, Scalar beta) {
    cout << "正在对z轴线性压缩" << endl;
    Scalar z_i, z_min = 0x3fffffff, z_max = -0x3fffffff;
    for (auto v : mesh.vertices()) {
        Point p = mesh.position(v);
        if (p[1] > z_max) z_max = p[1];
        if (p[1] < z_min) z_min = p[1];
    }

    cout << "z_max = " << z_max << ", z_min = " << z_min << endl;
    assert(z_min <= z_max);
    if (z_min < 0) {
        for (auto v : mesh.vertices()) {
            Point &p = mesh.position(v);
            p[1] -= z_min;
        }
        z_min = 0;
        z_max -= z_min;
    }

    for (auto v : mesh.vertices()) {
        Point &p = mesh.position(v);
        z_i = p[1];
        //执行线性压缩，原压缩公式似乎有问题，这里在beta之后加入了z_i
        p[1] = beta * z_i * (z_i - z_min) / (z_max - z_min);
        //cout << z_i - p[1] << endl;
    }
    cout << "压缩完成" << endl;
    return 0;
}
