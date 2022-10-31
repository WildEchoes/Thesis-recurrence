#include <easy3d/core/surface_mesh.h>
#include <cstring>

using namespace easy3d;

class Delaunay {
public:
	// 构造函数
	Delaunay(std::string, int);

	// flip edge which can not be able to ...
	void Flip(SurfaceMesh::Edge);

	// update vertices
	void updateV(SurfaceMesh::Vertex);

	// face area
	float area(SurfaceMesh::Face);

	// excenter of a triangle
	vec2 excenter(SurfaceMesh::Face);
	vec3 excenter2(SurfaceMesh::Face);

	// Delaunay Triangulation
	void Triangulation();

	// save
	void save_mesh();

	// 析构函数
	~Delaunay();

private:
	SurfaceMesh *mesh_;  // 平面mesh
	std::string name_;   //文件名
	int n_;              // 循环次数
	//int n_v, n_e, n_f;  //vertices edges faces
};

float angle(float, float, float);