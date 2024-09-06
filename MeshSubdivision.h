#ifndef _MeshSubdivision_
#define _MeshSubdivision_
#include "GeoData/MeshData.h"
#include "boost/smart_ptr.hpp"
#include <CGAL/subdivision_method_3.h>
#include <CGAL/Subdivision_method_3/subdivision_masks_3.h>
#include <CGAL/Subdivision_method_3/subdivision_hosts_3.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/Surface_mesh.h>

#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Polygon_mesh_processing/repair_polygon_soup.h>
#include <CGAL/Polygon_mesh_processing/orient_polygon_soup.h>
#include <CGAL/Polygon_mesh_processing/polygon_soup_to_polygon_mesh.h>
#include <CGAL/Polygon_mesh_processing/polygon_mesh_to_polygon_soup.h>

typedef CGAL::Exact_predicates_exact_constructions_kernel  K;
typedef K::FT  FT;
typedef std::array<FT, 3>                                       Custom_point;
typedef CGAL::Polyhedron_3<K> Polyhedron;
typedef K::Point_3 Point;
typedef CGAL::Surface_mesh<Point> PolygonMesh;
typedef typename boost::property_map<PolygonMesh, CGAL::vertex_point_t>::type Vertex_pmap;
namespace params = CGAL::parameters;
namespace sn3DCore
{
    class sn3DMeshData;
}
struct Array_traits
{
    struct Equal_3
    {
        bool operator()(const Custom_point& p, const Custom_point& q) const {
            return (p == q);
        }
    };
    struct Less_xyz_3
    {
        bool operator()(const Custom_point& p, const Custom_point& q) const {
            return std::lexicographical_compare(p.begin(), p.end(), q.begin(), q.end());
        }
    };
    Equal_3 equal_3_object() const { return Equal_3(); }
    Less_xyz_3 less_xyz_3_object() const { return Less_xyz_3(); }
};
class MeshSubdivision
{
public:
    MeshSubdivision(boost::shared_ptr<sn3DMeshData> mesh);
    ~MeshSubdivision();

public:
    //不改变模型轮廓细分
    void MidSubdivision();
    void CentroidSubdivision();

    //改变模型轮廓细分
    void LoopSubdivision();
    void LoopSubdivision2();
    void ButterflySubdivision();
    
    //局部细分
    void LocalSubdivision();
private:
    class PImpl;
    boost::shared_ptr<PImpl> impl_;
};

template <class Poly>
class WLoop_mask_3 {
    typedef Poly                                         PolygonMesh;
    typedef typename boost::graph_traits<PolygonMesh>::vertex_descriptor   vertex_descriptor;
    typedef typename boost::graph_traits<PolygonMesh>::halfedge_descriptor halfedge_descriptor;
    typedef typename boost::property_traits<Vertex_pmap>::value_type Point;
    typedef typename boost::property_traits<Vertex_pmap>::reference Point_ref;
    PolygonMesh& pmesh;
    Vertex_pmap vpm;
public:
    WLoop_mask_3(PolygonMesh& pmesh)
        : pmesh(pmesh), vpm(get(CGAL::vertex_point, pmesh))
    {}
    void edge_node(halfedge_descriptor hd, Point& pt) {
        Point_ref p1 = get(vpm, target(hd, pmesh));
        Point_ref p2 = get(vpm, target(opposite(hd, pmesh), pmesh));
        Point_ref f1 = get(vpm, target(next(hd, pmesh), pmesh));
        Point_ref f2 = get(vpm, target(next(opposite(hd, pmesh), pmesh), pmesh));
        pt = Point((3 * (p1[0] + p2[0]) + f1[0] + f2[0]) / 8,
            (3 * (p1[1] + p2[1]) + f1[1] + f2[1]) / 8,
            (3 * (p1[2] + p2[2]) + f1[2] + f2[2]) / 8);
    }
    void vertex_node(vertex_descriptor vd, Point& pt) {
        double R[] = { 0.0, 0.0, 0.0 };
        Point_ref S = get(vpm, vd);
        std::size_t n = 0;
        for (halfedge_descriptor hd : halfedges_around_target(vd, pmesh)) {
            ++n;
            Point_ref p = get(vpm, target(opposite(hd, pmesh), pmesh));
            R[0] += p[0];         R[1] += p[1];         R[2] += p[2];
        }
        if (n == 6) {
            pt = Point((10 * S[0] + R[0]) / 16, (10 * S[1] + R[1]) / 16, (10 * S[2] + R[2]) / 16);
        }
        else if (n == 3) {
            double B = (5.0 / 8.0 - std::sqrt(3 + 2 * std::cos(6.283 / n)) / 64.0) / n;
            double A = 1 - n * B;
            pt = Point((A * S[0] + B * R[0]), (A * S[1] + B * R[1]), (A * S[2] + B * R[2]));
        }
        else {
            double B = 3.0 / 8.0 / n;
            double A = 1 - n * B;
            pt = Point((A * S[0] + B * R[0]), (A * S[1] + B * R[1]), (A * S[2] + B * R[2]));
        }
    }
    void border_node(halfedge_descriptor hd, Point& ept, Point& vpt) {
        Point_ref ep1 = get(vpm, target(hd, pmesh));
        Point_ref ep2 = get(vpm, target(opposite(hd, pmesh), pmesh));
        ept = Point((ep1[0] + ep2[0]) / 2, (ep1[1] + ep2[1]) / 2, (ep1[2] + ep2[2]) / 2);
        Halfedge_around_target_circulator<Poly> vcir(hd, pmesh);
        Point_ref vp1 = get(vpm, target(opposite(*vcir, pmesh), pmesh));
        Point_ref vp0 = get(vpm, target(*vcir, pmesh));
        --vcir;
        Point_ref vp_1 = get(vpm, target(opposite(*vcir, pmesh), pmesh));
        vpt = Point((vp_1[0] + 6 * vp0[0] + vp1[0]) / 8,
            (vp_1[1] + 6 * vp0[1] + vp1[1]) / 8,
            (vp_1[2] + 6 * vp0[2] + vp1[2]) / 8);
    }
};
#endif