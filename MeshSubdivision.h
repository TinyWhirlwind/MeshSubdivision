#ifndef MESH_SUBDIVSION_H
#define MESH_SUBDIVSION_H

#include <vector>
#include <set>
#include "common/point3.h"
#include "GeoData/MeshData.h"
#include "boost/smart_ptr.hpp"
using boost::shared_ptr;
typedef std::pair<int, int> MeshEdge;
namespace sn3DCore
{
	class sn3DMeshData;
}

class MeshSubdivision
{
public:
	MeshSubdivision();
	MeshSubdivision(boost::shared_ptr<sn3DMeshData> mesh, float length = 0.3f);
	~MeshSubdivision();

	void setMesh(boost::shared_ptr<sn3DMeshData> mesh);

	boost::shared_ptr<sn3DMeshData> getMesh();

	void setLength(float length = 0.3f);

	void localMidPointSubdivision();

	void localButterflySubdivision();

	void globalSubdivision();

	void clacWeightPoint(int vertexDegree, Vertex* curVertex, MeshEdge curEdge, Point3f& wPoint);
public:
	std::map<MeshEdge, Face> all_halfEdges;
private:
	Face* getFFp(Face face, const MeshEdge edge);//找到与当前面片某一边相邻的面片

private:
	boost::shared_ptr<sn3DMeshData> m_mesh;
	float m_edgeLength = 0.f;
};

#endif
