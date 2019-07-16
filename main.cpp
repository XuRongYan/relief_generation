#include "Bas_ReliefGenerator.h"

using namespace std;
using namespace pmp;

const string OBJ_DIR = "/Users/xry/Desktop/lessons/ComputerGraphics/obj/"; //obj文件父路径
const string OBJ_NAME = "half_horse.obj";                                  //obj文件名
const string OBJ_RESULT_PATH = "./";
const string OBJ_RESULT_NAME = "res.obj";
const Scalar BETA = 0.1;

int main() {
    SurfaceMesh mesh;
    Bas_ReliefGenerator generator;
    generator.generate(mesh, OBJ_DIR + OBJ_NAME, OBJ_RESULT_PATH + OBJ_RESULT_NAME, BETA);
    return 0;
}