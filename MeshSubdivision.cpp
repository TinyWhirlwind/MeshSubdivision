#include "MeshSubdivision.h"
#include "GeoData/MeshTopology.h"
#include "FileFormat/sn3DIOManager.h"

#include <unordered_map>
#include <QDebug>

MeshSubdivision::MeshSubdivision()
{

}

MeshSubdivision::MeshSubdivision(boost::shared_ptr<sn3DMeshData> mesh, float length)
{
	m_mesh = mesh;
	m_edgeLength = length;
}

MeshSubdivision::~MeshSubdivision() {}

void MeshSubdivision::setMesh(boost::shared_ptr<sn3DMeshData> mesh)
{
	m_mesh = mesh;
}

boost::shared_ptr<sn3DMeshData> MeshSubdivision::getMesh()
{
	return m_mesh;
}

void MeshSubdivision::setLength(float length)
{
	m_edgeLength = length;
}

void MeshSubdivision::localMidPointSubdivision()
{
	//sn3DIOManager::instance().Write("D:\\data\\Subdivision\\initMesh.obj", mesh_.get());
	int faceNums = m_mesh->n_faces();
	for (int i = 0; i < faceNums; i++)
	{
		Face* curFace = m_mesh->F(i);
		if (curFace->V(0)->IsS() && curFace->V(1)->IsS() && curFace->V(2)->IsS())
		{
			curFace->SetS();
			float avgEdgeLength = 0;
			for (int i = 0; i < 3; i++)
			{
				Point3f v0 = curFace->V(i)->P();
				Point3f v1 = curFace->V((i + 1) % 3)->P();
				avgEdgeLength += (v0 - v1).Norm();
			}
			avgEdgeLength /= 3;
			if (avgEdgeLength <= 0.3f)
			{
				curFace->ClearS();
			}
		}
	}

	auto EdgeId = [](int a,int b) -> MeshEdge
	{
		if (a < b)
		{
			return std::make_pair(a, b);
		}
		else
		{
			return std::make_pair(b, a);
		}
	};

	std::map<MeshEdge, Face> ringFaceList;
	std::map<MeshEdge, Point3f> recordAddPts;
	for (auto& itor : m_mesh->m_faces)
	{
		if (itor.IsS())
		{
			itor.SetS();
			for (int i = 0; i < 3; i++)
			{
				Vertex* v1 = itor.V(i);
				Vertex* v2 = itor.V((i + 1) % 3);

				MeshEdge half_edge = std::make_pair(v1->m_index, v2->m_index);
				ringFaceList[half_edge] = itor;//有重复边

				MeshEdge edge = EdgeId(v1->m_index, v2->m_index);
				recordAddPts[edge] = (v1->P() + v2->P()) * 0.5f;//无重复边
			}
		}
	}
	if (ringFaceList.size() == 0)return;

	int num = m_mesh->n_vertices();
	//boost::shared_ptr<sn3DMeshData> aMesh(new sn3DMeshData);
	std::map<MeshEdge, int> addPtOnEdge;
	for (auto itor : recordAddPts)
	{
		m_mesh->AddVertex(itor.second);
		//aMesh->AddVertex(itor.second);
		addPtOnEdge[itor.first] = num;
		num++;
	}
	//sn3DIOManager::instance().Write("D:\\data\\Subdivision\\innerMesh.obj", aMesh.get());

	//细分边界面片 在ringFaceList中不重复为边界边
	//boost::shared_ptr<sn3DMeshData> bMesh(new sn3DMeshData);
	std::unordered_map<int, std::vector<MeshEdge>> boundaryPts;//记录每个边界面片上添加的新增点
	for (auto itor : ringFaceList)
	{
		MeshEdge swap_edge = std::make_pair(itor.first.second, itor.first.first);
		if (ringFaceList.find(swap_edge) == ringFaceList.end())
		{
			int v0 = itor.first.first;
			int v1 = itor.first.second;
			Face curFace = itor.second;
			auto ringFace = getFFp(curFace, itor.first);
			if (!ringFace || v0 < 0 || v1 < 0)
				continue;
			ringFace->SetM();
			boundaryPts[ringFace->m_index].push_back(EdgeId(v0, v1));//找到与边界边相邻的三角面片
			/*bMesh->AddVertex(ringFace->V(0)->P());
			bMesh->AddVertex(ringFace->V(1)->P());
			bMesh->AddVertex(ringFace->V(2)->P());
			qDebug() << "Edge_v0" << v0;
			qDebug() << "Edge_v0" << v1;
			qDebug() << "face_v0" << ringFace->V(0)->m_index;
			qDebug() << "face_v1" << ringFace->V(1)->m_index;
			qDebug() << "face_v2" << ringFace->V(2)->m_index;*/
		}
	}
	//sn3DIOManager::instance().Write("D:\\data\\Subdivision\\boundary_mesh.obj", bMesh.get());

	//链接细分面片
	boost::shared_ptr<sn3DMeshData> newMesh(new sn3DMeshData);
	for (auto itor : m_mesh->m_vertices)
	{
		newMesh->AddVertex(itor.P());
	}

	std::vector<int> faceList;
	for (auto& itor : m_mesh->m_faces)
	{
		int v0 = itor.V(0)->m_index;
		int v1 = itor.V(1)->m_index;
		int v2 = itor.V(2)->m_index;
		if (itor.IsS())
		{
			MeshEdge edge0 = EdgeId(v0, v1);
			MeshEdge edge1 = EdgeId(v1, v2);
			MeshEdge edge2 = EdgeId(v2, v0);
			int addV0 = addPtOnEdge[edge0];
			int addV1 = addPtOnEdge[edge1];
			int addV2 = addPtOnEdge[edge2];
			newMesh->AddFace(v0, addV0, addV2);
			newMesh->AddFace(addV0, v1, addV1);
			newMesh->AddFace(addV2, addV1, v2);
			newMesh->AddFace(addV0, addV1, addV2);
		}
		else if(itor.IsM())
		{
			int faceID = itor.m_index;
			faceList.push_back(faceID);

			int addNums = boundaryPts[faceID].size();
			if (addNums == 3)
			{
				int addV0 = addPtOnEdge[EdgeId(v0, v1)];
				int addV1 = addPtOnEdge[EdgeId(v1, v2)];
				int addV2 = addPtOnEdge[EdgeId(v2, v0)];
				
				newMesh->AddFace(v1, addV1, addV0);
				newMesh->AddFace(addV1, addV2, addV0);
				newMesh->AddFace(addV1, v2, addV2);
				newMesh->AddFace(addV0, addV2, v0);
			}
			else if(addNums == 2)
			{
				MeshEdge e0 = boundaryPts[faceID][0];
				MeshEdge e1 = boundaryPts[faceID][1];
				int addV0 = addPtOnEdge[EdgeId(e0.first, e0.second)];
				int addV1 = addPtOnEdge[EdgeId(e1.first, e1.second)];
				float linkLength0 = 0;
				float linkLength1 = 0;
				for (int i = 0; i < 3; i++)
				{
					if (itor.m_v[i] != e0.first && itor.m_v[i] != e0.second)
					{
						//与e0相对的点
						Point3f topPt = m_mesh->V(addV0)->P();
						linkLength0 = (itor.V(i)->P() - topPt).Norm();
						//将边的方向与三角形边顺序保持一直
						e0.first = itor.m_v[(i + 1) % 3];
						e0.second = itor.m_v[(i + 2) % 3];
					}
					if (itor.m_v[i] != e1.first && itor.m_v[i] != e1.second)
					{
						//与e0相对的点
						Point3f topPt = m_mesh->V(addV1)->P();
						linkLength1 = (itor.V(i)->P() - topPt).Norm();
						//将边的方向与三角形边顺序保持一直
						e1.first = itor.m_v[(i + 1) % 3];
						e1.second = itor.m_v[(i + 2) % 3];
					}
				}
				if (linkLength0 <= linkLength1)
				{
					//先连接e0后连接e1
					if (e1.second == e0.first)
					{
						newMesh->AddFace(e0.first, addV0, addV1);
						newMesh->AddFace(addV0, e1.first, addV1);
						newMesh->AddFace(addV0, e0.second, e1.first);
					}
					if (e0.second == e1.first)
					{
						newMesh->AddFace(e0.first, addV0, e1.second);
						newMesh->AddFace(addV0, e0.second, addV1);
						newMesh->AddFace(addV0, addV1, e1.second);
					}
				}
				else
				{
					//先连接e1后连接e0
					if (e0.second == e1.first)
					{
						newMesh->AddFace(e1.first, addV1, addV0);
						newMesh->AddFace(addV1, e0.first, addV0);
						newMesh->AddFace(addV1, e1.second, e0.first);
					}
					if (e1.second == e0.first)
					{
						newMesh->AddFace(e1.first, addV1, e0.second);
						newMesh->AddFace(addV1, e1.second, addV0);
						newMesh->AddFace(addV1, addV0, e0.second);
					}
				}
			}
			else if (addNums == 1)
			{
				MeshEdge edge = boundaryPts[faceID].front();
				int addV = addPtOnEdge[edge];
				if ((v0 == edge.first && v2 == edge.second) ||
					(v2 == edge.first && v0 == edge.second))
				{
					newMesh->AddFace(addV, v0, v1);
					newMesh->AddFace(addV, v1, v2);
				}
				else if (v2 == edge.first && v1 == edge.second ||
					(v1 == edge.first && v2 == edge.second))
				{
					newMesh->AddFace(addV, v2, v0);
					newMesh->AddFace(addV, v0, v1);
				}
				else if (v1 == edge.first && v0 == edge.second ||
					(v0 == edge.first && v1 == edge.second))
				{
					newMesh->AddFace(addV, v2, v0);
					newMesh->AddFace(addV, v1, v2);
				}
				else
				{
					qDebug() << "zpp: subdivision error!!!";
					return;
				}
			}
			else
			{
				qDebug() << "zpp: subdivision error!!!";
				return;
			}
		}
		else
		{
			newMesh->AddFace(v0, v1, v2);
		}
	}

	newMesh->UpdateFaceTopology();
	newMesh->UpdateFaceIdentity();
	newMesh->DirtyTopology();
	newMesh->UpdateNormal();

	m_mesh = newMesh;
	sn3DIOManager::instance().Write("D:\\data\\Subdivision\\subdivision_mesh.obj", m_mesh.get());
}

void MeshSubdivision::localButterflySubdivision()
{
	int faceNums = m_mesh->n_faces();
	for (int i = 0; i < faceNums; i++)
	{
		Face* curFace = m_mesh->F(i);
		if (curFace->V(0)->IsS() && curFace->V(1)->IsS() && curFace->V(2)->IsS())
		{
			curFace->SetS();
			float avgEdgeLength = 0;
			for (int i = 0; i < 3; i++)
			{
				Point3f v0 = curFace->V(i)->P();
				Point3f v1 = curFace->V((i + 1) % 3)->P();
				avgEdgeLength += (v0 - v1).Norm();
			}
			avgEdgeLength /= 3;
			if (avgEdgeLength <= 0.3f)
			{
				curFace->ClearS();
			}
		}
	}

	auto EdgeId = [](int a, int b) -> MeshEdge
		{
			if (a < b)
			{
				return std::make_pair(a, b);
			}
			else
			{
				return std::make_pair(b, a);
			}
		};

	std::map<MeshEdge, Face> subdvision_half_edges;//细分所有有向边 有重复边
	std::map<MeshEdge, int> subdvision_whole_edges;//细分所有无向边 无重复边
	std::map<MeshEdge, Point3f> edgePoints;//细分所有边上的插值点
	all_halfEdges.clear();//记录所有半边
	for (auto itor : m_mesh->m_faces)
	{
		for (int i = 0; i < 3; i++)
		{
			Vertex* from_vertex = m_mesh->V(i);
			Vertex* to_vertex = m_mesh->V((i + 1) % 3);
			MeshEdge half_edge = std::make_pair(from_vertex->m_index, to_vertex->m_index);
			MeshEdge whole_edge = EdgeId(from_vertex->m_index, to_vertex->m_index);
			if (itor.IsS())
			{
				subdvision_half_edges[half_edge] = itor;
				subdvision_whole_edges[whole_edge]++;
			}
			all_halfEdges[half_edge] = itor;
		}
	}
	for (auto itor : subdvision_whole_edges)
	{
		for (int i = 0; i < 3; i++)
		{
			Vertex* from_vertex = m_mesh->V(itor.first.first);
			Vertex* to_vertex = m_mesh->V(itor.first.second);
			MeshEdge whole_edge = itor.first;
			Point3f newInterpolationP = Point3f{ 0.f,0.f,0.f };
			if (edgePoints.find(whole_edge) != edgePoints.end())
				continue;

			VFIterator from_viter(from_vertex);
			VFIterator to_viter(to_vertex);
			if (from_vertex->IsB() && to_vertex->IsB())//边界边上插入
			{
				//边界边上插值 四点法
				Vertex* from_bb_vertex;
				Vertex* to_bb_vertex;
				for (; !from_viter.End(); ++from_viter)
				{
					Vertex* v1 = from_viter.f->V((from_viter.z + 1) % 3);
					Vertex* v2 = from_viter.f->V((from_viter.z + 2) % 3);
					if (v1->IsB() && v1->m_index != to_vertex->m_index)
					{
						from_bb_vertex = v1;
					}
					if (v2->IsB() && v2->m_index != to_vertex->m_index)
					{
						from_bb_vertex = v2;
					}
				}
				for (; !to_viter.End(); ++to_viter)
				{
					Vertex* v1 = to_viter.f->V((to_viter.z + 1) % 3);
					Vertex* v2 = to_viter.f->V((to_viter.z + 2) % 3);
					if (v1->IsB() && v1->m_index != from_vertex->m_index)
					{
						to_bb_vertex = v1;
					}
					if (v2->IsB() && v2->m_index != from_vertex->m_index)
					{
						to_bb_vertex = v2;
					}
				}
				newInterpolationP = (to_vertex->P() + from_vertex->P()) * 9 / 16 - (to_bb_vertex->P() + from_bb_vertex->P()) * 1 / 16;
			}
			else
			{
				//求左右顶点的度
				std::set<int> from_ringVertexs;
				std::set<int> to_ringVertexs;
				for (; !from_viter.End(); ++from_viter)
				{
					Vertex* v1 = from_viter.f->V((from_viter.z + 1) % 3);
					Vertex* v2 = from_viter.f->V((from_viter.z + 2) % 3);
					from_ringVertexs.insert(v1->m_index);
					from_ringVertexs.insert(v2->m_index);
				}
				for (; !to_viter.End(); ++to_viter)
				{
					Vertex* v1 = to_viter.f->V((to_viter.z + 1) % 3);
					Vertex* v2 = to_viter.f->V((to_viter.z + 2) % 3);
					to_ringVertexs.insert(v1->m_index);
					to_ringVertexs.insert(v2->m_index);
				}
				int from_degree = from_ringVertexs.size();
				int to_degree = to_ringVertexs.size();
				if (from_degree == 0 || to_degree == 0)
					continue;
				if (from_degree == 6 && to_degree == 6)
				{
					Point3f to_newInterpolationP;
					Point3f from_newInterpolationP;
					clacWeightPoint(6, to_vertex, whole_edge, to_newInterpolationP);
					clacWeightPoint(6, from_vertex, whole_edge, from_newInterpolationP);
					newInterpolationP = to_newInterpolationP + from_newInterpolationP;
				}
				else if (from_degree == 6 && to_degree != 6)
				{
					clacWeightPoint(to_degree, to_vertex, whole_edge, newInterpolationP);
					
				}
				else if (from_degree != 6 && to_degree == 6)
				{
					clacWeightPoint(from_degree, from_vertex, whole_edge, newInterpolationP);
				}
				else
				{
					Point3f to_newInterpolationP;
					Point3f from_newInterpolationP;
					clacWeightPoint(to_degree, to_vertex, whole_edge, to_newInterpolationP);
					clacWeightPoint(from_degree, from_vertex, whole_edge, from_newInterpolationP);
					newInterpolationP = (to_newInterpolationP + from_newInterpolationP).Normalize();
				}
			}
			edgePoints[whole_edge] = newInterpolationP;
		}
	}

	std::map<MeshEdge, int> interpolationPts_index;
	int num = m_mesh->n_vertices();
	for (auto itor : edgePoints)
	{
		m_mesh->AddVertex(itor.second);
		interpolationPts_index[itor.first] = num;
		num++;
	}

	std::unordered_map<int, std::vector<MeshEdge>> boundaryPts;//记录每个边界面片上添加的新增点
	for (auto itor : subdvision_half_edges)
	{
		MeshEdge swap_edge = std::make_pair(itor.first.second, itor.first.first);
		if (subdvision_half_edges.find(swap_edge) == subdvision_half_edges.end())
		{
			int v0 = itor.first.first;
			int v1 = itor.first.second;
			Face curFace = itor.second;
			auto ringFace = getFFp(curFace, itor.first);
			if (!ringFace || v0 < 0 || v1 < 0)
				continue;
			ringFace->SetM();
			boundaryPts[ringFace->m_index].push_back(EdgeId(v0, v1));//找到与边界边相邻的三角面片
		}
	}


	//链接细分面片
	boost::shared_ptr<sn3DMeshData> newMesh(new sn3DMeshData);
	for (auto itor : m_mesh->m_vertices)
	{
		newMesh->AddVertex(itor.P());
	}

	std::vector<int> faceList;
	for (auto& itor : m_mesh->m_faces)
	{
		int v0 = itor.V(0)->m_index;
		int v1 = itor.V(1)->m_index;
		int v2 = itor.V(2)->m_index;
		if (itor.IsS())
		{
			MeshEdge edge0 = EdgeId(v0, v1);
			MeshEdge edge1 = EdgeId(v1, v2);
			MeshEdge edge2 = EdgeId(v2, v0);
			int addV0 = interpolationPts_index[edge0];
			int addV1 = interpolationPts_index[edge1];
			int addV2 = interpolationPts_index[edge2];
			newMesh->AddFace(v0, addV0, addV2);
			newMesh->AddFace(addV0, v1, addV1);
			newMesh->AddFace(addV2, addV1, v2);
			newMesh->AddFace(addV0, addV1, addV2);
		}
		else if (itor.IsM())
		{
			int faceID = itor.m_index;
			faceList.push_back(faceID);

			int addNums = boundaryPts[faceID].size();
			if (addNums == 3)
			{
				int addV0 = interpolationPts_index[EdgeId(v0, v1)];
				int addV1 = interpolationPts_index[EdgeId(v1, v2)];
				int addV2 = interpolationPts_index[EdgeId(v2, v0)];

				newMesh->AddFace(v1, addV1, addV0);
				newMesh->AddFace(addV1, addV2, addV0);
				newMesh->AddFace(addV1, v2, addV2);
				newMesh->AddFace(addV0, addV2, v0);
			}
			else if (addNums == 2)
			{
				MeshEdge e0 = boundaryPts[faceID][0];
				MeshEdge e1 = boundaryPts[faceID][1];
				int addV0 = interpolationPts_index[EdgeId(e0.first, e0.second)];
				int addV1 = interpolationPts_index[EdgeId(e1.first, e1.second)];
				float linkLength0 = 0;
				float linkLength1 = 0;
				for (int i = 0; i < 3; i++)
				{
					if (itor.m_v[i] != e0.first && itor.m_v[i] != e0.second)
					{
						//与e0相对的点
						Point3f topPt = m_mesh->V(addV0)->P();
						linkLength0 = (itor.V(i)->P() - topPt).Norm();
						//将边的方向与三角形边顺序保持一直
						e0.first = itor.m_v[(i + 1) % 3];
						e0.second = itor.m_v[(i + 2) % 3];
					}
					if (itor.m_v[i] != e1.first && itor.m_v[i] != e1.second)
					{
						//与e0相对的点
						Point3f topPt = m_mesh->V(addV1)->P();
						linkLength1 = (itor.V(i)->P() - topPt).Norm();
						//将边的方向与三角形边顺序保持一直
						e1.first = itor.m_v[(i + 1) % 3];
						e1.second = itor.m_v[(i + 2) % 3];
					}
				}
				if (linkLength0 <= linkLength1)
				{
					//先连接e0后连接e1
					if (e1.second == e0.first)
					{
						newMesh->AddFace(e0.first, addV0, addV1);
						newMesh->AddFace(addV0, e1.first, addV1);
						newMesh->AddFace(addV0, e0.second, e1.first);
					}
					if (e0.second == e1.first)
					{
						newMesh->AddFace(e0.first, addV0, e1.second);
						newMesh->AddFace(addV0, e0.second, addV1);
						newMesh->AddFace(addV0, addV1, e1.second);
					}
				}
				else
				{
					//先连接e1后连接e0
					if (e0.second == e1.first)
					{
						newMesh->AddFace(e1.first, addV1, addV0);
						newMesh->AddFace(addV1, e0.first, addV0);
						newMesh->AddFace(addV1, e1.second, e0.first);
					}
					if (e1.second == e0.first)
					{
						newMesh->AddFace(e1.first, addV1, e0.second);
						newMesh->AddFace(addV1, e1.second, addV0);
						newMesh->AddFace(addV1, addV0, e0.second);
					}
				}
			}
			else if (addNums == 1)
			{
				MeshEdge edge = boundaryPts[faceID].front();
				int addV = interpolationPts_index[edge];
				if ((v0 == edge.first && v2 == edge.second) ||
					(v2 == edge.first && v0 == edge.second))
				{
					newMesh->AddFace(addV, v0, v1);
					newMesh->AddFace(addV, v1, v2);
				}
				else if (v2 == edge.first && v1 == edge.second ||
					(v1 == edge.first && v2 == edge.second))
				{
					newMesh->AddFace(addV, v2, v0);
					newMesh->AddFace(addV, v0, v1);
				}
				else if (v1 == edge.first && v0 == edge.second ||
					(v0 == edge.first && v1 == edge.second))
				{
					newMesh->AddFace(addV, v2, v0);
					newMesh->AddFace(addV, v1, v2);
				}
				else
				{
					qDebug() << "zpp: subdivision error!!!";
					return;
				}
			}
			else
			{
				qDebug() << "zpp: subdivision error!!!";
				return;
			}
		}
		else
		{
			newMesh->AddFace(v0, v1, v2);
		}
	}

	newMesh->UpdateFaceTopology();
	newMesh->UpdateFaceIdentity();
	newMesh->DirtyTopology();
	newMesh->UpdateNormal();

	m_mesh = newMesh;
}

void MeshSubdivision::clacWeightPoint(int vertexDegree, Vertex* curVertex, MeshEdge curEdge, Point3f& wPoint)
{
	if (vertexDegree < 3)
		return;
	Vertex* toVertex;
	if (curVertex->m_index == curEdge.first)
	{
		toVertex = m_mesh->V(curEdge.second);
	}
	else
	{
		toVertex = m_mesh->V(curEdge.first);
	}
	VFIterator viter(curVertex);
	VFIterator to_viter(toVertex);
	if (vertexDegree == 3)
	{
		Vertex* ee_vertex1;
		Vertex* ee_vertex2;
		MeshEdge bb_edge_CW = curEdge;
		MeshEdge bb_edge_CCW = std::make_pair(curEdge.second, curEdge.first);
		Face bb_face0 = all_halfEdges[bb_edge_CW];
		Face bb_face1 = all_halfEdges[bb_edge_CCW];
		for (int i = 0; i < 3; i++)
		{
			if (bb_face0.m_v[i] != curVertex->m_index && bb_face0.m_v[i] != toVertex->m_index)
			{
				ee_vertex1 = bb_face0.V(i);
			}
			if (bb_face1.m_v[i] != curVertex->m_index && bb_face1.m_v[i] != toVertex->m_index)
			{
				ee_vertex2 = bb_face1.V(i);
			}
		}
		wPoint = toVertex->P() * 5 / 12 - (ee_vertex1->P() + ee_vertex2->P()) * 1 / 12;
	}
	else if (vertexDegree == 4)
	{
		Vertex* ee_vertex1;
		Vertex* ee_vertex2;
		Vertex* ee_vertex3;
		MeshEdge bb_edge_CW = std::make_pair(curVertex->m_index, toVertex->m_index);
		MeshEdge bb_edge_CCW = std::make_pair(toVertex->m_index, curVertex->m_index);
		Face bb_face0 = all_halfEdges[bb_edge_CW];
		Face bb_face1 = all_halfEdges[bb_edge_CCW];
		for (int i = 0; i < 3; i++)
		{
			if (bb_face0.m_v[i] != toVertex->m_index && bb_face0.m_v[i] != curVertex->m_index)
			{
				ee_vertex1 = bb_face0.V(i);
			}
			if (bb_face1.m_v[i] != toVertex->m_index && bb_face1.m_v[i] != curVertex->m_index)
			{
				ee_vertex3 = bb_face1.V(i);
			}
		}

		MeshEdge ee_Edge0 = std::make_pair(ee_vertex1->m_index, toVertex->m_index);
		Face ee_face = all_halfEdges[ee_Edge0];
		if (ee_face.m_index == bb_face0.m_index)
		{
			std::swap(ee_Edge0.first, ee_Edge0.second);
		}
		for (int i = 0; i < 3; i++)
		{
			if (ee_face.m_v[i] != ee_Edge0.first && ee_face.m_v[i] != ee_Edge0.second)
			{
				ee_vertex2 = ee_face.V(i);
			}
		}
		wPoint = curVertex->P() * 3 / 8 - ee_vertex2->P() * 1 / 8;
	}
	else if (vertexDegree >= 5)
	{
		if (vertexDegree == 6)
		{
			MeshEdge bb_edge_CW = std::make_pair(curVertex->m_index, toVertex->m_index);
			MeshEdge bb_edge_CCW = std::make_pair(toVertex->m_index, curVertex->m_index);
			Face bb_face0 = all_halfEdges[bb_edge_CW];
			Face bb_face1 = all_halfEdges[bb_edge_CCW];
			Vertex* bb_vertex0;
			Vertex* bb_vertex1;
			for (int i = 0; i < 3; i++)
			{
				if (bb_face0.m_v[i] != curEdge.first && bb_face0.m_v[i] != curEdge.second)
				{
					bb_vertex0 = bb_face0.V(i);
				}
				if (bb_face1.m_v[i] != curEdge.first && bb_face1.m_v[i] != curEdge.second)
				{
					bb_vertex1 = bb_face1.V(i);
				}
			}
			MeshEdge cc_edge_CW0 = std::make_pair(curVertex->m_index, bb_vertex0->m_index);
			MeshEdge cc_edge_CW1 = std::make_pair(curVertex->m_index, bb_vertex1->m_index);
			Face cc_face0 = all_halfEdges[cc_edge_CW0];
			Face cc_face1 = all_halfEdges[cc_edge_CW1];
			if (cc_face0.m_index == bb_face0.m_index)
			{
				std::swap(cc_edge_CW0.first, cc_edge_CW0.second);
				cc_face0 = all_halfEdges[cc_edge_CW0];
			}
			if (cc_face1.m_index == bb_face1.m_index)
			{
				std::swap(cc_edge_CW1.first, cc_edge_CW1.second);
				cc_face1 = all_halfEdges[cc_edge_CW1];
			}
			Vertex* cc_vertex0;
			Vertex* cc_vertex1;
			for (int i = 0; i < 3; i++)
			{
				if (cc_face0.m_v[i] != cc_edge_CW0.first && cc_face0.m_v[i] != cc_edge_CW0.second)
				{
					cc_vertex0 = cc_face0.V(i);
				}
				if (cc_face1.m_v[i] != cc_edge_CW1.first && cc_face1.m_v[i] != cc_edge_CW1.second)
				{
					cc_vertex1 = cc_face1.V(i);
				}
			}
			//这里b点一半 方便左右两边求和
			wPoint = curVertex->P() * 1 / 2 + (bb_vertex0->P() + bb_vertex1->P()) * 1 / 16 - (cc_vertex0->P() + cc_vertex1->P()) * 1 / 16;
		}
		else
		{
			for (auto& itor : m_mesh->m_vertices)
			{
				itor.ClearS();
			}
			std::vector<int> ee_indexs;
			for (; !viter.End(); ++viter)
			{
				Vertex* v1 = viter.f->V((viter.z + 1) % 3);
				Vertex* v2 = viter.f->V((viter.z + 2) % 3);
				if (!v1->IsS())
				{
					ee_indexs.push_back(v1->m_index);
				}
				if (!v2->IsS())
				{
					ee_indexs.push_back(v2->m_index);
				}
				v1->SetS();
				v2->SetS();
			}
			int startIndex = -1;
			int ringSize = ee_indexs.size();
			for (int i = 0; i < ee_indexs.size(); i++)
			{
				if (ee_indexs[i] == toVertex->m_index)
				{
					startIndex = i;
					break;
				}
			}
			std::vector<int> ee_points;
			for (int i = startIndex; i < startIndex + ringSize; i++)
			{
				ee_points.push_back(ee_indexs[i % ringSize]);
			}
			for (int i = 0; i < ee_points.size(); i++)
			{
				float ei_w = (1 / 4 + cos(2 * i * M_PI / ringSize) + 1 / 2 * cos(4 * i * M_PI)) / ringSize;
				wPoint += m_mesh->V(ee_points[i])->P() * ei_w;
			}
		}
		
	}
}

Face* MeshSubdivision::getFFp(Face face, const MeshEdge edge)
{
	Face* ringFace = nullptr;
	Vertex* from_vertex = m_mesh->V(edge.first);
	Vertex* to_vertex = m_mesh->V(edge.second);
	VFIterator vfi(from_vertex);
	for (; !vfi.End(); ++vfi)
	{
		Vertex* v1 = vfi.f->V((vfi.z + 1) % 3);
		Vertex* v2 = vfi.f->V((vfi.z + 2) % 3);
		if (v2->m_index == to_vertex->m_index)
		{
			ringFace = vfi.f;
			break;
		}
	}
	return ringFace;
}

void MeshSubdivision::globalSubdivision()
{
	//更新边
}
