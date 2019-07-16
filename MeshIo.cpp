//
// Created by 徐溶延 on 2019-07-16.
//

#include "MeshIo.h"

/**
 * 读取mesh
 * @param mesh
 * @return
 */
int MeshIo::readMesh(SurfaceMesh &mesh, const string &path) {
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
int MeshIo::saveMesh(const SurfaceMesh &mesh, string path) {
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
int MeshIo::showMesh(const string &path) {
    MeshViewer viewer(path.data(), 800, 450);
    viewer.load_mesh(path.data());
    return viewer.run();
}
