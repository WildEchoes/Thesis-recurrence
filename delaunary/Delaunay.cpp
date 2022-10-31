#include "Delaunay.h"
#include <iostream>
#include <easy3d/fileio/surface_mesh_io.h>
#include <cmath>
#include <Eigen/Dense>

using namespace std;
using namespace Eigen;

#define pi 3.14159265358979323846

Delaunay::Delaunay(string name, int n)
{
	const string file_name = "DATA/" + name + ".obj";
	mesh_ = SurfaceMeshIO::load(file_name);
	if (mesh_)
		cout << "load OK! " << endl;
	name_ = name;
	n_ = n;
}

void Delaunay::Flip(SurfaceMesh::Edge e)
{
	// e is border continue
	if (mesh_->is_border(e) || !mesh_->is_flip_ok(e))
		return;

	// compute a1 a2, they are the two opposite corners of the edge(e)
	float len_e = mesh_->edge_length(e);

	auto h10 = mesh_->halfedge(e, 0);  //c
	auto h11 = mesh_->next(h10);       //a
	auto h12 = mesh_->next(h11);       //b

	float a1 = angle(mesh_->edge_length(mesh_->edge(h11)), mesh_->edge_length(mesh_->edge(h12)), len_e);

	auto h20 = mesh_->halfedge(e, 1);   //c
	auto h21 = mesh_->next(h20);		//a
	auto h22 = mesh_->next(h21);		//b

	float a2 = angle(mesh_->edge_length(mesh_->edge(h21)), mesh_->edge_length(mesh_->edge(h22)), len_e);

	// a = a1 + a2, if a > pi  --> flip edge
	float a = a1 + a2;
	if (a > pi)
		mesh_->flip(e);
}

void Delaunay::updateV(SurfaceMesh::Vertex v)
{
	// v is border -> continue
	if (mesh_->is_border(v))
		return;

	// update vertex position
	float total_area = 0;
	vec3 C = { 0.0,0.0,0.0 };
	for (auto f : mesh_->faces(v))
	{
		float f_area = area(f);
		total_area += f_area;

		vec3 cj = excenter2(f);
		C = C + f_area * cj;
	}
	C = C / total_area;

	// position 只支持三维坐标 只需将z设为0
	mesh_->position(v) = vec3(C[0], C[1], 0.0);
}

Delaunay::~Delaunay()
{
	delete mesh_;
}

// compute face area
float Delaunay::area(SurfaceMesh::Face f)
{
	float x[3], y[3];
	int i = 0;
	for (auto v : mesh_->vertices(f))
	{
		vec2 p = mesh_->position(v);
		x[i] = p.x;
		y[i++] = p.y;
	}
	float area = 0.5 *abs((x[0] * y[1] + x[1] * y[2] + x[2] * y[0] - x[0] * y[2] - x[1] * y[0] - x[2] * y[1]));
	return area;
}

vec2 Delaunay::excenter(SurfaceMesh::Face f)
{
	vec3 v[3];
	int i = 0;
	for (auto vv : mesh_->vertices(f))
		v[i++] = mesh_->position(vv);
	float x1, x2, x3, y1, y2, y3;
	x1 = v[0].x;   y1 = v[0].y;
	x2 = v[1].x;   y2 = v[1].y;
	x3 = v[2].x;   y3 = v[2].y;

	Matrix3f X, Y, xy;
	X << x1*x1 + y1 * y1, y1, 1,
		x2*x2 + y2 * y2, y2, 1,
		x3*x3 + y3 * y3, y3, 1;
	xy << x1, y1, 1,
		x2, y2, 1,
		x3, y3, 1;
	Y << x1, x1*x1 + y1 * y1, 1,
		x2, x2*x2 + y2 * y2, 1,
		x3, x3*x3 + y3 * y3, 1;

	float h_x, h_y, h_xy;
	h_x = X.determinant();
	h_y = Y.determinant();
	h_xy = xy.determinant();

	float x = 0.5*h_x / h_xy;
	float y = 0.5*h_y / h_xy;

	return vec2(x, y);
}

vec3 Delaunay::excenter2(SurfaceMesh::Face f)
{
	vec3 v[3];
	int i = 0;
	for (auto vv : mesh_->vertices(f))
		v[i++] = mesh_->position(vv);

	float x1, y1, x2, y2, x3, y3;
	x1 = v[0][0];   y1 = v[0][1];
	x2 = v[1][0];   y2 = v[1][1];
	x3 = v[2][0];   y3 = v[2][1];

	float a1, b1, c1, a2, b2, c2;
	a1 = 2 * (x2 - x1);	  a2 = 2 * (x3 - x2);	c1 = x2 * x2 + y2 * y2 - x1 * x1 - y1 * y1;
	b1 = 2 * (y2 - y1);	  b2 = 2 * (y3 - y2);	c2 = x3 * x3 + y3 * y3 - x2 * x2 - y2 * y2;

	float x = (b2 * c1 - b1 * c2) / (a1 * b2 - a2 * b1);
	float y = (a1 * c2 - a2 * c1) / (a1 * b2 - a2 * b1);

	return vec3(x, y, 0.0);
}

void Delaunay::Triangulation()
{
	// n为循环次数
	int i = 0;
	while (i++ < n_) 
	{
		for (SurfaceMesh::Edge e : mesh_->edges())
			Flip(e);
		for (SurfaceMesh::Vertex v : mesh_->vertices())
			updateV(v);
	}
	save_mesh();
}

void Delaunay::save_mesh()
{
	// Write the mesh to a new file.
	const std::string save_file_name = "./Result/" + name_ + "-" + to_string(n_) + ".obj";
	if (SurfaceMeshIO::save(save_file_name, mesh_))
		std::cout << "mesh saved to \'" << save_file_name << "\'" << std::endl;
	else
		std::cerr << "failed create the new file" << std::endl;
}

// return arccos(c)
float angle(float a, float b, float c)
{
	float cos_c = (a*a + b * b - c * c) / (2 * a*b);
	return (float)acos(cos_c);
}