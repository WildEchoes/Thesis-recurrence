#include <easy3d/core/surface_mesh.h>
#include <cstring>

using namespace easy3d;

class Delaunay {
public:
	// ���캯��
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

	// ��������
	~Delaunay();

private:
	SurfaceMesh *mesh_;  // ƽ��mesh
	std::string name_;   //�ļ���
	int n_;              // ѭ������
	//int n_v, n_e, n_f;  //vertices edges faces
};

float angle(float, float, float);