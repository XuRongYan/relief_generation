#include <iostream>
#include <string>
#include <vector>
#include <pmp/visualization/MeshViewer.h>
#include <pmp/SurfaceMesh.h>
#include <pmp/SurfaceMeshIO.h>
#include <pmp/algorithms/NormalCone.h>
#include <pmp/algorithms/SurfaceNormals.h>
#include <pmp/algorithms/SurfaceRemeshing.h>
#include <pmp/algorithms/SurfaceSimplification.h>
#include <pmp/algorithms/DifferentialGeometry.h>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <GLFrame.h>

using namespace std;
using namespace pmp;


const string OBJ_DIR = "/Users/xry/Desktop/lessons/ComputerGraphics/obj/"; //obj文件父路径
const string OBJ_NAME = "half_horse.obj";                                  //obj文件名
const string OBJ_PLANE_NAME = "plane.obj";
const string OBJ_RESULT_PATH = "./";
const string OBJ_RESULT_NAME = "res.obj";
const string OBJ_RESULT_REMESH = "remesh.obj";
const string PROPERTY_VERTEX_NEIGHBOR = "vertex_around_vertex";
const string PROPERTY_FACE_NEIGHBOR = "face_around_vertex";
const Scalar BETA = 0.1;

typedef SurfaceMesh::FaceAroundVertexCirculator FaceNeighbors;
typedef SurfaceMesh::VertexAroundVertexCirculator VertexNeighbors;
typedef Eigen::SparseMatrix<Scalar> SparseMatrixType;

Scalar avg_gradient;

int readMesh(SurfaceMesh &, const string &);

int saveMesh(const SurfaceMesh &mesh, string path);

int showMesh(const string &);

int init(SurfaceMesh &);

int getFaceNeighborProperty(SurfaceMesh &);

Scalar remeshing(SurfaceMesh &);

int simplification(SurfaceMesh &, Scalar);

Eigen::Vector3f get_gradient_Phi(const SurfaceMesh &, const Face &, const Vertex &);

Eigen::MatrixXf generate_b(SurfaceMesh &);

SparseMatrixType generate_Laplace(const SurfaceMesh &mesh);

int solve_poison(SurfaceMesh &, const SparseMatrixType &, const Eigen::MatrixXf &);

int compression(SurfaceMesh &, Scalar);

int optimizeHighSlopes(SurfaceMesh &mesh);

Eigen::Vector3f getGradientField(const SurfaceMesh &mesh, const Face &f);


/**
 * 读取mesh
 * @param mesh
 * @return
 */
int readMesh(SurfaceMesh &mesh, const string &path) {
    cout << "正在读取mesh……, 路径：" << path << endl;
    bool success = mesh.read(path);
    if (success) {
        cout << "mesh读取成功!" << endl;
        cout << "vertices:" << mesh.n_vertices() << ", edges:" << mesh.n_edges() << ", faces:" << mesh.n_faces()
             << endl;
    } else {
        cout << "mesh读取失败！" << endl;
        exit(-1);
    }
    return 0;
}

/**
 * 保存mesh
 * @param mesh
 * @return
 */
int saveMesh(const SurfaceMesh &mesh, string path) {
    cout << "正在保存mesh……" << endl;
    bool success = mesh.write(path);
    cout << (success ? "保存成功" : "保存失败") << endl;
    return 0;
}

/**
 * 显示mesh
 * @param path
 * @return
 */
int showMesh(const string &path) {
    MeshViewer viewer(path.data(), 800, 450);
    viewer.load_mesh(path.data());
    viewer.draw(OBJ_DIR + OBJ_PLANE_NAME);
    return viewer.run();
}

/**
 * 获取顶点周围面的数量
 * @param mesh
 * @return
 */
int getFaceNeighborProperty(SurfaceMesh &mesh) {
    auto face_neighbor = mesh.add_vertex_property<size_t>(PROPERTY_FACE_NEIGHBOR);
    for (auto v : mesh.vertices()) {
        auto neighbors = FaceNeighbors(&mesh, v);
        face_neighbor[v] = neighbors.end() - neighbors.begin() + 1;
    }
    return 0;
}

int dragZ(SurfaceMesh& mesh, Scalar zpp) {
    for (auto v : mesh.vertices()) {
        Point &p = mesh.position(v);
        p[1] += zpp;
    }
    return 0;
}

/**
 * 检查顶点是否按id顺序遍历
 * @param mesh
 * @return
 */
int checkId(const SurfaceMesh &mesh) {
    int i = 0;
    for (auto v : mesh.vertices()) {
        if (v.idx() != i) {
            throw "disorderly vertices";
        }
    }
    return 0;
}

/**
 * 初始化
 * @param mesh
 * @return
 */
int init(SurfaceMesh &mesh) {

    getFaceNeighborProperty(mesh);
    try {
        checkId(mesh);
    } catch (const char *msg) {
        cout << msg << endl;
    }
    return 0;
}

/**
 * remesh，取网格的平均边长作为参数
 * @param mesh
 * @return
 */
Scalar remeshing(SurfaceMesh &mesh) {
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

/**
 * 简化网格，当三角形面片数量大于门限值时对网格进行简化
 * @param mesh
 * @param avg_scalar
 * @return
 */
int simplification(SurfaceMesh &mesh, const Scalar avg_scalar) {
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

/**
 * 获取梯度场
 * @param mesh
 * @param face
 * @param point
 * @return
 */
Eigen::Vector3f get_gradient_Phi(const SurfaceMesh &mesh, const Face &face, const Vertex &point) {
    vector<Point> points;
    int index = 0, i = 0;
    for (auto v : mesh.vertices(face)) {
        if (v == point) index = i;
        Point p = mesh.position(v);
        if(mesh.is_boundary(v)) p[1] = 0;
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

Eigen::Vector3f getGradientField(const SurfaceMesh &mesh, const Face &f) {
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

/**
 * 生成目标矩阵场2n*n
 * @param mesh
 * @return
 */
Eigen::MatrixXf generate_b(SurfaceMesh &mesh) {
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

/**
 * 生成拉普拉斯矩阵
 * @param mesh
 * @return
 */
SparseMatrixType generate_Laplace(const SurfaceMesh &mesh) {
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

/**
 * 最小二乘法求解方程
 * @param mesh
 * @param Ls
 * @param b
 * @return
 */
int solve_poison(SurfaceMesh &mesh, const SparseMatrixType &Ls, const Eigen::MatrixXf &b) {
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

/**
 * 对z轴线性压缩
 * @param mesh
 * @param beta
 * @return
 */
int compression(SurfaceMesh &mesh, Scalar beta) {
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

bool innerAngle(const SurfaceMesh &mesh, const Face &f) {
    vector<Vertex> vertices;
    int i = 0;
    for (auto v : mesh.vertices(f)) {
        vertices[i++] = v;
    }
    return true;
}


int main() {
    SurfaceMesh mesh, plane;
    readMesh(mesh, OBJ_DIR + OBJ_NAME);
    //dragZ(mesh, 2);
    Scalar avg_scalar = remeshing(mesh);
    simplification(mesh, avg_scalar);
    saveMesh(mesh, OBJ_RESULT_PATH + OBJ_RESULT_REMESH);
    init(mesh);
    Eigen::SparseMatrix<Scalar> Ls = generate_Laplace(mesh);
    Eigen::MatrixXf b = generate_b(mesh);
    solve_poison(mesh, Ls, b);
    compression(mesh, BETA);
    saveMesh(mesh, OBJ_RESULT_PATH + OBJ_RESULT_NAME);
    showMesh(OBJ_RESULT_PATH + OBJ_RESULT_NAME);
    return 0;
}