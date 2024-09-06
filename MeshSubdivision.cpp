#include "MeshSubdivision.h"
#include "GeoData/MeshTopology.h"
#include "FileFormat/sn3DIOManager.h"

class MeshSubdivision::PImpl
{
public:
    PImpl(boost::shared_ptr<sn3DMeshData> mesh)
    {
        m_mesh = mesh;
    }
    ~PImpl() {}
private:
    boost::shared_ptr<sn3DMeshData> m_mesh;
public:
    void MidSubdivision();
    void CentroidSubdivision();
    void LoopSubdivision();
    void LoopSubdivision2();
    void ButterflySubdivision();
	void LocalSubdivision();
};

void MeshSubdivision::MidSubdivision()
{
    impl_->MidSubdivision();
}

void MeshSubdivision::CentroidSubdivision()
{
    impl_->CentroidSubdivision();
}

void MeshSubdivision::LoopSubdivision()
{
    impl_->LoopSubdivision();
}

void MeshSubdivision::LoopSubdivision2()
{
    impl_->LoopSubdivision2();
}


void MeshSubdivision::ButterflySubdivision()
{
    impl_->ButterflySubdivision();
}

void MeshSubdivision::LocalSubdivision()
{
	impl_->LocalSubdivision();
}

void MeshSubdivision::PImpl::MidSubdivision()
{
	typedef std::pair<int, int> Edge;
	//边排序
	auto edgeId = [](int a, int b)
		{
			if (a > b)
				return std::make_pair(a, b);
			else
				return  std::make_pair(b, a);
		};

	boost::shared_ptr<sn3DMeshData> newMesh(new sn3DMeshData);
	std::map<Edge, Point3f> edgeToNewPt;//边-新点
	std::map<Edge, int> edgeToNewIndex;//边-新点Index
	for (auto itor :m_mesh->m_vertices)
	{
		newMesh->AddVertex(itor.P());
	}
	for (auto itor : m_mesh->m_faces)
	{
		for (int i = 0; i < 3; i++)
		{
			Edge e = edgeId(itor.V(i)->m_index, itor.V((i + 1) % 3)->m_index);
			Point3f p = (itor.V(i)->P() + itor.V((i + 1) % 3)->P()) * 0.5f;
			edgeToNewPt[e] = p;
		}
	}
	int startId = m_mesh->n_vertices();
	for (auto& itor : edgeToNewPt)
	{
		newMesh->AddVertex(itor.second);
		edgeToNewIndex[itor.first] = startId;
		startId++;
	}
	for (auto itor : m_mesh->m_faces)
	{
		int p0 = itor.V(0)->m_index;
		int p1 = itor.V(1)->m_index;
		int p2 = itor.V(2)->m_index;
		Edge e0 = edgeId(p0, p1);
		Edge e1 = edgeId(p1, p2);
		Edge e2 = edgeId(p2, p0);
		int new_p0 = edgeToNewIndex[e0];
		int new_p1 = edgeToNewIndex[e1];
		int new_p2 = edgeToNewIndex[e2];
		newMesh->AddFace(p0, new_p0, new_p2);
		newMesh->AddFace(new_p0, new_p1, new_p2);
		newMesh->AddFace(new_p0, p1, new_p1);
		newMesh->AddFace(new_p2, new_p1, p2);
	}
	m_mesh = newMesh;
}

void MeshSubdivision::PImpl::CentroidSubdivision()
{
	boost::shared_ptr<sn3DMeshData> newMesh(new sn3DMeshData);
	std::map<int, Point3f> faceToNewPt;//面-新点
	std::map<int, int> faceToNewIndex;//面-新点Index
	for (auto& itor : m_mesh->m_vertices)
	{
		newMesh->AddVertex(itor.P());
	}
	for (auto& itor : m_mesh->m_faces)
	{
		float avgEdgeLength = 0.f;
		for (int i = 0; i < 3; i++)
		{
			avgEdgeLength += (itor.V(i)->P() - itor.V((i + 1) % 3)->P()).Norm();
		}
		avgEdgeLength /= 3;
		if (avgEdgeLength < 0.45f)continue;
		Point3f p = (itor.V(0)->P() + itor.V(1)->P() + itor.V(2)->P()) / 3;
		faceToNewPt[itor.m_index] = p;
	}
	int startId = m_mesh->n_vertices();
	for (auto& itor : faceToNewPt)
	{
		newMesh->AddVertex(itor.second);
		faceToNewIndex[itor.first] = startId;
		startId++;
	}
	for (auto& itor : m_mesh->m_faces)
	{
		int p0 = itor.V(0)->m_index;
		int p1 = itor.V(1)->m_index;
		int p2 = itor.V(2)->m_index;
		int new_center = faceToNewIndex[itor.m_index];
		newMesh->AddFace(p0, new_center, p2);
		newMesh->AddFace(p0, p1, new_center);
		newMesh->AddFace(p1, p2, new_center);
	}
	m_mesh = newMesh;
}

void MeshSubdivision::PImpl::LoopSubdivision()
{
	std::vector<std::array<FT, 3> > points;
	std::vector<std::vector<std::size_t>> triangles;
	int vn = m_mesh->m_vertices.size();
	int fn = m_mesh->m_faces.size();
	for (int i = 0; i < vn; i++)
	{
		points.push_back(CGAL::make_array<FT>(m_mesh->V(i)->P().x, m_mesh->V(i)->P().y, m_mesh->V(i)->P().z));
	}
	for (int i = 0; i < fn; i++)
	{
		std::vector<std::size_t> tri;
		tri.push_back(m_mesh->F(i)->m_v[0]); tri.push_back(m_mesh->F(i)->m_v[1]); tri.push_back(m_mesh->F(i)->m_v[2]);
		triangles.push_back(tri);
	}
	Polyhedron polyhedron;
	std::vector<std::array<FT, 3>> outpoints;
	std::vector<std::vector<std::size_t>> outtriangles;
	CGAL::Polygon_mesh_processing::orient_polygon_soup(points, triangles);
	CGAL::Polygon_mesh_processing::repair_polygon_soup(points, triangles, CGAL::parameters::geom_traits(Array_traits()));
	CGAL::Polygon_mesh_processing::polygon_soup_to_polygon_mesh(points, triangles, polyhedron);
	CGAL::Subdivision_method_3::Loop_subdivision(polyhedron, 1);
	m_mesh->Clear();
	Point3f curPoint;
	for (int i = 0; i < outpoints.size(); i++)
	{
		int j = 0;
		for (auto& itor : outpoints[i])
		{
			curPoint[j] = CGAL::to_double(itor);
			j++;
		}
		m_mesh->AddVertex(curPoint);
	}
	for (int i = 0; i < outtriangles.size(); i++)
	{
		std::vector<int> curFace;
		for (auto itor : outtriangles[i])
		{
			curFace.push_back(static_cast<int>(itor));
		}
		m_mesh->AddFace(curFace);
	}
	m_mesh->DirtyTopology();
	m_mesh->DirtyBoundary();
	m_mesh->UpdateNormal();
	m_mesh->UpdateFaceTopology();
}

void MeshSubdivision::PImpl::LoopSubdivision2()
{
	PolygonMesh pMesh;
	int vn = m_mesh->m_vertices.size();
	int fn = m_mesh->m_faces.size();
	std::map<int, PolygonMesh::Vertex_index> vertex_descriptor;
	for (int i = 0; i < vn; i++)
	{
		PolygonMesh::Vertex_index ptor = pMesh.add_vertex(Point(m_mesh->V(i)->P().x, m_mesh->V(i)->P().y, m_mesh->V(i)->P().z));
		vertex_descriptor[i] = ptor;
	}
	for (int i = 0; i < fn; i++)
	{
		PolygonMesh::Vertex_index v0ptor = vertex_descriptor[m_mesh->F(i)->m_v[0]];
		PolygonMesh::Vertex_index v1ptor = vertex_descriptor[m_mesh->F(i)->m_v[1]];
		PolygonMesh::Vertex_index v2ptor = vertex_descriptor[m_mesh->F(i)->m_v[2]];
		pMesh.add_face(v0ptor, v1ptor, v2ptor);
	}

	// 创建 PQQ_stencil_3 对象
	typedef CGAL::PQQ_stencil_3<PolygonMesh, Vertex_pmap> PQQ_stencil;
	Vertex_pmap tva = get(CGAL::vertex_point, pMesh);
	CGAL::PQQ_stencil_3<PolygonMesh, Vertex_pmap> mask(&pMesh, tva);
	CGAL::Subdivision_method_3::PTQ(pMesh, mask, params::number_of_iterations(1));
	//CGAL::Subdivision_method_3::PTQ(pMesh, WLoop_mask_3<PolygonMesh>(pMesh), params::number_of_iterations(1));

	std::vector<Point3f> outpoints;
	std::vector<std::vector<int>> outtriangles;
	for (auto itor = pMesh.vertices_begin(); itor != pMesh.vertices_end(); ++itor)
	{
		double x = CGAL::to_double(pMesh.point(*itor).x());
		double y = CGAL::to_double(pMesh.point(*itor).y());
		double z = CGAL::to_double(pMesh.point(*itor).z());
		outpoints.push_back(Point3f(x, y, z));
	}
	for (auto itor = pMesh.faces_begin(); itor != pMesh.faces_end(); ++itor)
	{
		std::vector<int> face;
		auto hec = pMesh.halfedges_around_face(pMesh.halfedge(*itor)).begin();
		do {
			face.push_back(pMesh.target(*hec).idx());
		} while (++hec != pMesh.halfedges_around_face(pMesh.halfedge(*itor)).end());
		outtriangles.push_back(face);
	}

	m_mesh->Clear();
	Point3f curPoint;
	for (int i = 0; i < outpoints.size(); i++)
	{
		int j = 0;
		m_mesh->AddVertex(curPoint);
	}
	for (int i = 0; i < outtriangles.size(); i++)
	{
		std::vector<int> curFace;
		for (auto itor : outtriangles[i])
		{
			curFace.push_back(static_cast<int>(itor));
		}
		m_mesh->AddFace(curFace);
	}

	m_mesh->DirtyTopology();
	m_mesh->DirtyBoundary();
	m_mesh->UpdateNormal();
	m_mesh->UpdateFaceTopology();
}

void MeshSubdivision::PImpl::ButterflySubdivision()
{
    
}

void MeshSubdivision::PImpl::LocalSubdivision()
{
	std::vector<Face*> faces;
	typedef std::pair<int, int> Edge;
	//边排序
	auto edgeId = [](int a, int b)
		{
			if (a > b)
				return std::make_pair(a, b);
			else
				return  std::make_pair(b, a);
		};


	std::map<Edge, Point3f> newPointList;//边对应的新增点
	//mesh不需要细分的面片复制到newMesh 点全部复制 并计算新加点
	boost::shared_ptr<sn3DMeshData> newMesh(new sn3DMeshData);
	for (auto& itor : m_mesh->m_vertices)
	{
		newMesh->AddVertex(itor.P());
	}

	for (auto& itor : m_mesh->m_faces)
	{
		//边长过短不细分
		float e0Length = (itor.V(0)->P() - itor.V(1)->P()).Norm();
		float e1Length = (itor.V(1)->P() - itor.V(2)->P()).Norm();
		float e2Length = (itor.V(2)->P() - itor.V(0)->P()).Norm();
		if (itor.V(0)->IsS() || itor.V(1)->IsS() || itor.V(2)->IsS())
		{
			if (e0Length >= 0.45f || e1Length >= 0.45f || e2Length >= 0.45f)
			{
				Point3f addP0 = (itor.V(0)->P() + itor.V(1)->P()) * 0.5f;
				Point3f addP1 = (itor.V(1)->P() + itor.V(2)->P()) * 0.5f;
				Point3f addP2 = (itor.V(2)->P() + itor.V(0)->P()) * 0.5f;
				Edge e0 = edgeId(itor.V(0)->m_index, itor.V(1)->m_index);
				Edge e1 = edgeId(itor.V(1)->m_index, itor.V(2)->m_index);
				Edge e2 = edgeId(itor.V(2)->m_index, itor.V(0)->m_index);
				newPointList[e0] = addP0;
				newPointList[e1] = addP1;
				newPointList[e2] = addP2;
				continue;
			}
			continue;
		}
	}
	for (auto& itor : m_mesh->m_faces)
	{
		if (itor.V(0)->IsS() || itor.V(1)->IsS() || itor.V(2)->IsS())
		{
			itor.V(0)->SetS();
			itor.V(1)->SetS();
			itor.V(2)->SetS();
		}
		newMesh->AddFace(itor.m_v[0], itor.m_v[1], itor.m_v[2]);
	}
	for (auto& itor : m_mesh->m_faces)
	{
		if (itor.V(0)->IsS() || itor.V(1)->IsS() || itor.V(2)->IsS())
		{
			faces.push_back(&itor);
			continue;
		}
		newMesh->AddFace(itor.m_v[0], itor.m_v[1], itor.m_v[2]);
	}


	std::map<Edge, int> newPointIndex;//边对应的新增点的索引
	// 添加点
	int starId = newMesh->n_vertices();
	for (auto& pt : newPointList)
	{
		newMesh->AddVertex(pt.second);
		newPointIndex[pt.first] = starId;
		starId++;
	}
	
	// 添加面片
	std::vector<Edge> addPtEdgeList;//一个面片三条边含有新加点的边集
	for (auto itor : faces)
	{
		//统计面片新增点的个数
		addPtEdgeList.clear();
		int addPtEdgeNum = 0;
		for (int i = 0; i < 3; i++)
		{
			Edge e = edgeId(itor->V(i)->m_index, itor->V((i + 1) % 3)->m_index);
			if (newPointIndex.count(e))
			{
				addPtEdgeNum++;
				addPtEdgeList.push_back(e);
			}
		}

		int addPtIndex = -1;
		int linkPtIndex = -1;
		int nextPtIndex = -1;
		int otherPtIndex = -1;
		int secondAddPtIndex = -1;
		if (addPtEdgeNum == 0)
		{
			int p0 = itor->V(0)->m_index;
			int p1 = itor->V(1)->m_index;
			int p2 = itor->V(2)->m_index;
			newMesh->AddFace(p0, p1, p2);
		}
		else if (addPtEdgeNum == 1)
		{
			addPtIndex = newPointIndex[addPtEdgeList[0]];
			for (int i = 0; i < 3; i++)
			{
				if (itor->V(i)->m_index != addPtEdgeList[0].first && itor->V(i)->m_index != addPtEdgeList[0].second)
				{
					linkPtIndex = itor->V(i)->m_index;
					nextPtIndex = itor->V((i + 1) % 3)->m_index;
					otherPtIndex = itor->V((i + 2) % 3)->m_index;
					break;
				}
			}

			newMesh->AddFace(addPtIndex, otherPtIndex, linkPtIndex);
			newMesh->AddFace(addPtIndex, linkPtIndex, nextPtIndex);
		}
		else if (addPtEdgeNum == 2)
		{
			addPtIndex = newPointIndex[addPtEdgeList[0]];
			for (int i = 0; i < 3; i++)
			{
				if (itor->V(i)->m_index != addPtEdgeList[0].first && itor->V(i)->m_index != addPtEdgeList[0].second)
				{
					linkPtIndex = itor->V(i)->m_index;
					nextPtIndex = itor->V((i + 1) % 3)->m_index;
					otherPtIndex = itor->V((i + 2) % 3)->m_index;
					break;
				}
			}
			Edge nextEdge = edgeId(otherPtIndex, linkPtIndex);
			if (!newPointIndex.count(nextEdge))
			{
				secondAddPtIndex = newPointIndex[edgeId(linkPtIndex, nextPtIndex)];
				newMesh->AddFace(addPtIndex, otherPtIndex, linkPtIndex);
				newMesh->AddFace(addPtIndex, linkPtIndex, secondAddPtIndex);
				newMesh->AddFace(addPtIndex, secondAddPtIndex, nextPtIndex);
			}
			else
			{
				secondAddPtIndex = newPointIndex[edgeId(linkPtIndex, otherPtIndex)];
				newMesh->AddFace(addPtIndex, otherPtIndex, secondAddPtIndex);
				newMesh->AddFace(addPtIndex, secondAddPtIndex, linkPtIndex);
				newMesh->AddFace(addPtIndex, linkPtIndex, nextPtIndex);
			}
		}
		else if (addPtEdgeNum == 3)
		{
			int p0 = itor->V(0)->m_index;
			int p1 = itor->V(1)->m_index;
			int p2 = itor->V(2)->m_index;
			Edge e0 = edgeId(p0, p1);
			Edge e1 = edgeId(p1, p2);
			Edge e2 = edgeId(p2, p0);
			int new_p0 = newPointIndex[e0];
			int new_p1 = newPointIndex[e1];
			int new_p2 = newPointIndex[e2];

			newMesh->AddFace(p0, new_p0, new_p2);
			newMesh->AddFace(new_p0, new_p1, new_p2);
			newMesh->AddFace(new_p0, p1, new_p1);
			newMesh->AddFace(new_p2, new_p1, p2);
		}
		else
		{
			continue;
		}
	}

	newMesh->UpdateVertexTopology();
	newMesh->UpdateNormal();
	newMesh->UpdateEdge();
	newMesh->UpdateFaceTopology();
	newMesh->UpdateFaceTexCoord();
	m_mesh = newMesh;
	sn3DIOManager::instance().Write("D:\\data\\TempData\\initmesh.obj", newMesh.get());
}

