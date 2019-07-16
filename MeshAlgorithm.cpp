//
// Created by 徐溶延 on 2019-07-16.
//

#include "MeshAlgorithm.h"

Scalar MeshAlgorithm::remeshing(SurfaceMesh &mesh) {
    cout << "正在进行remesh" << endl;
    SurfaceRemeshing remesh(mesh);
    Scalar avg_scalar = 0;
    for (auto e : mesh.edges()) {
        avg_scalar += mesh.edge_length(e);
    }
    avg_scalar /= mesh.n_edges();
    cout << "average edge length = " << avg_scalar << endl;
    remesh.uniform_remeshing(avg_scalar);
    cout << "remesh完成" << endl;
    cout << "vertices:" << mesh.n_vertices() << ", edges:" << mesh.n_edges() << ", faces:" << mesh.n_faces()
         << endl;
    return avg_scalar;
}

int MeshAlgorithm::simplification(SurfaceMesh &mesh, Scalar avg_scalar) {
    cout << "正在进行simplification" << endl;
    SurfaceSimplification simplification1(mesh);
    simplification1.initialize(5, avg_scalar, 10, 10, 0.001);
    int n = mesh.n_vertices() <= 10e3 ? mesh.n_vertices() : 5e3 + mesh.n_vertices() / 2;
    cout << "simplify之后的顶点数量为：" << n << endl;
    simplification1.simplify(n);
    cout << "simplify完成" << endl;
    cout << "vertices:" << mesh.n_vertices() << ", edges:" << mesh.n_edges() << ", faces:" << mesh.n_faces()
         << endl;
    return 0;
}

Eigen::Vector3f MeshAlgorithm::get_gradient_Phi(const SurfaceMesh &mesh, const Face &face, const Vertex &point) {
    vector<Point> points;
    int index = 0, i = 0;
    for (auto v : mesh.vertices(face)) {
        if (v == point) index = i;
        Point p = mesh.position(v);
        if (mesh.is_boundary(v)) p[1] = 0;
        points.push_back(p);
        i++;
    }
    assert(points.size() == 3);
    Eigen::Matrix3f m_left, m_right;
    m_right << 1, 0, -1,
            0, 1, -1,
            0, 0, 0;
    Point v1 = points[0] - points[2];
    Point v2 = points[1] - points[2];
    Normal n = SurfaceNormals::compute_face_normal(mesh, face);
    Eigen::Vector3f ev1, ev2, en;
    ev1 << v1[0], v1[1], v1[2];
    ev2 << v2[0], v2[1], v2[2];
    en << n[0], n[1], n[2];
    m_left << ev1.transpose(),
            ev2.transpose(),
            en.transpose();
    Eigen::Matrix3f res = m_left.inverse() * m_right;
    return res.col(index);
}

Eigen::MatrixXf MeshAlgorithm::generate_b(SurfaceMesh &mesh) {
    cout << "正在生成目标矩阵场" << endl;
    const int n = mesh.n_vertices();
    Eigen::VectorXf divw = Eigen::VectorXf::Zero(n);
    int idx = 0;
    //n*3
    for (auto v : mesh.vertices()) {
        Scalar res = 0;
        for (auto f : FaceNeighbors(&mesh, v)) {
            Point p[3];
            Scalar area = 0;
            if(mesh.is_boundary(f)) {
                int i = 0;
                for(auto vv: mesh.vertices(f)) {
                    p[i] = mesh.position(vv);
                    if(mesh.is_boundary(vv)) p[i][1] = 0;
                    i++;
                }
                area = triangle_area(p[0], p[1], p[2]);
            } else {
                area = triangle_area(mesh, f);
            }
            Eigen::Vector3f grad_ti = get_gradient_Phi(mesh, f, v);
            res += area * grad_ti.transpose() * getGradientField(mesh, f);
        }
        divw.row(idx++) << res;
    }
    cout << "目标矩阵场生成完毕" << endl;
    return divw;
}

int MeshAlgorithm::solve_poison(SurfaceMesh &mesh, const SparseMatrixType &Ls, const Eigen::MatrixXf &b) {
    Eigen::SparseMatrix<Scalar> ls_transpose = Ls.transpose();
    Eigen::SparseMatrix<Scalar> LsLs = ls_transpose * Ls;
    Eigen::SimplicialCholesky<Eigen::SparseMatrix<Scalar> > MatricsCholesky(LsLs);
    Eigen::VectorXf xyzRHS = ls_transpose * b;
    Eigen::VectorXf xyz = MatricsCholesky.solve(xyzRHS);
    for (auto v : mesh.vertices()) {
        Point &p = mesh.position(v);
        p[1] = xyz[v.idx()];
    }
    return 0;
}

Eigen::Vector3f MeshAlgorithm::getGradientField(const SurfaceMesh &mesh, const Face &f) {
    Vertex vi, vj, vk;
    int i = 0;
    for (auto v : mesh.vertices(f)) {
        if (i == 0) vi = v;
        else if (i == 1) vj = v;
        else vk = v;
        i++;
    }
    Point pi, pj, pk;
    pi = mesh.position(vi);
    pj = mesh.position(vj);
    pk = mesh.position(vk);
    if(mesh.is_boundary(vi)) pi[1] = 0;
    if(mesh.is_boundary(vj)) pj[1] = 0;
    if(mesh.is_boundary(vk)) pk[1] = 0;
    Eigen::Vector3f vpi, vpj, vpk;
    vpi << pi[0], pi[1], pi[2];
    vpj << pj[0], pj[1], pj[2];
    vpk << pk[0], pk[1], pk[2];
    Scalar vpij = vpj[1] - vpi[1];
    Scalar vpik = vpk[1] - vpi[1];
    Eigen::Vector3f v1 = vpij * get_gradient_Phi(mesh, f, vj);
    Eigen::Vector3f v2 = vpik * get_gradient_Phi(mesh, f, vk);
    Eigen::Vector3f w = v1 + v2;
    return w;

}

SparseMatrixType MeshAlgorithm::generate_Laplace(const SurfaceMesh &mesh) {
    cout << "正在生成拉普拉斯矩阵" << endl;
    int vn = mesh.n_vertices();
    int count0 = 0;
    vector<int> begin_N(vn);
    for (auto v : mesh.vertices()) {
        begin_N[v.idx()] = count0;
        count0 += mesh.valence(v) + 1;
    }
    typedef Eigen::Triplet<Scalar> T;
    vector<T> tripletList(count0);
    //n*n
    for (auto v : mesh.vertices()) {
        int i = 0;
        Scalar sum_weight = 0;
        for (auto nei_v : VertexNeighbors(&mesh, v)) {
            Edge e = mesh.find_edge(v, nei_v);
            Scalar weight = 0.5 * cotan_weight(mesh, e);
            tripletList[begin_N[v.idx()] + i + 1] = T(v.idx(), nei_v.idx(), -weight);
            i++;
            sum_weight += weight;
        }
        tripletList[begin_N[v.idx()]] = T(v.idx(), v.idx(), sum_weight);
    }
    Eigen::SparseMatrix<Scalar> Ls(vn, vn);
    Ls.setFromTriplets(tripletList.begin(), tripletList.end());
    cout << "拉普拉斯矩阵生成完毕" << endl;
    return Ls;
}
