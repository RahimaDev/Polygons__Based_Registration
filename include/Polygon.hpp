#include"Estimation.hpp"
#include "LiteGeom/LgDistance3.hpp"

#include "LiteGeom/LgPoint3.hpp"
#include "LiteGeom/LgPoint2.hpp"
#include <CGAL/General_polygon_with_holes_2.h>
#include <pcl/surface/concave_hull.h>


#include <iostream>
#include <fstream>


#include <CGAL/Boolean_set_operations_2.h>
#include <CGAL/intersections.h>
#include <boost/config.hpp>
#include <boost/version.hpp>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>


#include <CGAL/Polygon_2.h>
#include <CGAL/Polygon_with_holes_2.h>

#include <CGAL/algorithm.h>
#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/Alpha_shape_2.h>

#include <CGAL/Constrained_Delaunay_triangulation_2.h>
#include <CGAL/Constrained_triangulation_plus_2.h>





typedef CGAL::Exact_predicates_exact_constructions_kernel K;
typedef K::FT FT;
typedef K::Point_2  Point;
typedef K::Segment_2  Segment;
typedef K::Vector_2 Vector;

typedef K::Triangle_2 Triangle;
typedef CGAL::Polygon_2<K>                           Polygon;
typedef CGAL::Polygon_with_holes_2<K>                Polygon_with_holes;
typedef CGAL::Alpha_shape_vertex_base_2<K> Vb;
typedef CGAL::Alpha_shape_face_base_2<K>  Fb;
typedef CGAL::Triangulation_data_structure_2<Vb,Fb> Tds;
typedef CGAL::Delaunay_triangulation_2<K,Tds> Triangulation_2;

typedef CGAL::Alpha_shape_2<Triangulation_2>  Alpha_shape_2;
typedef Alpha_shape_2::Alpha_shape_edges_iterator Alpha_shape_edges_iterator;
typedef Alpha_shape_2::Alpha_shape_vertices_iterator Alpha_shape_vertices_iterator;
typedef Alpha_shape_2::Face_handle Face_handle;
typedef Alpha_shape_2::Edge Edge;
typedef Alpha_shape_2::Vertex_handle Vertex_handle;
typedef Polygon_with_holes::Hole_const_iterator  Hole_const_iterator;
typedef Polygon_with_holes::Hole_iterator   Hole_iterator;


/////////// project a point on a plane
Lg::Point3f  projection(Lg::Point3f  p, Lg::Point3f O,Lg::Point3f N)
{

    Lg::Point3f  proj;

    Lg::Point3f V=p-O;
    double dist=V*N;
    proj=p-(dist*N);

    return proj;

}
//////////////////

Lg::Point3f best_m(Lg::Point3f n)

{
    Lg::Point3f x(1,0,0), y(0,1,0), z(0,0,1);
    Lg::Point3f nx=n^x, ny=n^y;
    float nnx = nx.Norm2(), nny=ny.Norm2();
    if(nnx>nny) return nx;
    else return ny;


}
////////////////// Alpha shape
std::list<Polygon_with_holes> Alpha_shape(  pcl::ModelCoefficients::Ptr  pl, pcl::PointCloud<pcl::PointXYZ>::Ptr cloud,float alpha)
{        
    Lg::Point3f n(pl->values[0], pl->values[1],pl->values[2]);
    //n.Normalize();
    Lg::Point3f m =  best_m(n);
    m.Normalize();
    Lg::Point3f k = n^m;
    k.Normalize();

    vector<Lg::Point2f> points_2d(cloud->points.size());
    std::vector<K::Point_2> points_cgal(points_2d.size());
    // TODO CALCULER O
    // Lg::Point3f O = -pl->values[3]/n.Norm()*n;
    Lg::Point3f O=( - pl->values[3]/((pl->values[0]*pl->values[0])+(pl->values[1]*pl->values[1])+(pl->values[2]*pl->values[2])))*n;

    for(int i = 0; i < cloud->points.size(); ++i)
    {
        Lg::Point3f P(cloud->points[i].x,cloud->points[i].y,cloud->points[i].z);

        Lg::Point3f  P0=P-O;
        //float xp=P0*m, yp=P0*k;
        float xp=P0*m, yp=P0*k;
        points_2d[i] = Lg::Point2f(xp, yp); // CgalPoint(xp, yp);
        points_cgal[i] = Point(xp,yp);
    }
    // ALPHASHAPE CGAL INTEGRATION
    Alpha_shape_2 A(points_cgal.begin(), points_cgal.end(),
                    FT(alpha),
                    Alpha_shape_2::REGULARIZED);
    std::cout << "##### Join triangles #####" << std::endl;

    std::list<Polygon> triangles;
    for(typename Alpha_shape_2::Finite_faces_iterator fit = A.finite_faces_begin();
        fit != A.finite_faces_end();
        ++fit){

        if(A.classify(fit) == Alpha_shape_2::INTERIOR){
            Triangle triangle = A.triangle(fit);
            Polygon poly;
            poly.push_back(triangle.vertex(0));
            poly.push_back(triangle.vertex(1));
            poly.push_back(triangle.vertex(2));
            triangles.push_back(poly);
        }
    }

    std::list<Polygon_with_holes> res;
    CGAL::join(triangles.begin(), triangles.end(), std::back_inserter (res));
    for(typename std::list<Polygon_with_holes>::iterator it_ring = res.begin();
        it_ring != res.end(); it_ring++){
        //CGAL::draw(*it_ring);
        Polygon outer = it_ring->outer_boundary();
        //std::cout << "=== Outer ===" << std::endl;

        for(typename Polygon::Vertex_const_iterator it_vertex = outer.vertices_begin();
            it_vertex != outer.vertices_end(); it_vertex++){
            //std::cout << it_vertex->x() << "\t" << it_vertex->y() << std::endl;
        }
        for(typename Polygon_with_holes::Hole_const_iterator it_hole = it_ring->holes_begin();
            it_hole != it_ring->holes_end(); it_hole++)
        {
            Polygon inner = *it_hole;
            // std::cout << "=== Inner ===" << std::endl;
            for(typename Polygon::Vertex_const_iterator it_vertex = inner.vertices_begin();
                it_vertex != inner.vertices_end(); it_vertex++){
                //std::cout << it_vertex->x() << "\t" << it_vertex->y() << std::endl;
            }}

    }
    return res;
}

//// 3D model of 2D polygon

poly_holes Polygon3D_creation(Polygon_with_holes& p,pcl::ModelCoefficients::Ptr  pl)
{poly_holes Res;

    Lg::Point3f n(pl->values[0], pl->values[1],pl->values[2]);
    Lg::Point3f m = best_m(n);
    m.Normalize();
    Lg::Point3f k = n^m;
    k.Normalize();
    Lg::Point3f O=( - pl->values[3]/((pl->values[0]*pl->values[0])+(pl->values[1]*pl->values[1])+(pl->values[2]*pl->values[2])))*n;
    std::vector<Point> VEC,VEC1;
    std::vector<std::vector<Point>>VEC2;
    Polygon outer = p.outer_boundary();
    ///////
    for(typename Polygon::Vertex_const_iterator it_vertex = outer.vertices_begin();
        it_vertex != outer.vertices_end(); it_vertex++){

        VEC1.push_back(Point( it_vertex->x(),it_vertex->y() ));
    }
    for(typename Polygon_with_holes::Hole_const_iterator it_hole = p.holes_begin();
        it_hole != p.holes_end(); it_hole++)
    {
        Polygon inner = *it_hole;
        std::vector<Point>V;
        for(typename Polygon::Vertex_const_iterator it_vertex = inner.vertices_begin();
            it_vertex != inner.vertices_end(); it_vertex++){

            V.push_back(Point( it_vertex->x(),it_vertex->y() ));
        }
        VEC2.push_back(V);}



    pcl::PointCloud<pcl::PointXYZ>::Ptr polygon_cloud(new pcl::PointCloud<pcl::PointXYZ>);
    pcl::PointCloud<pcl::PointXYZ>::Ptr out(new pcl::PointCloud<pcl::PointXYZ>);
    std::vector<pcl::PointCloud<pcl::PointXYZ>::Ptr> Hol;


    out->width    = VEC1.size();
    out->height   = 1;
    out->is_dense = false;
    for(int i=0; i<VEC1.size(); i++)
    {Lg::Point3f A;
        A =O +CGAL::to_double( VEC1[i].x())*m +CGAL::to_double( VEC1[i].y())*k;

        pcl::PointXYZ B=pcl::PointXYZ(A.x(),A.y(), A.z());
        out->points.push_back(B);

    }

    for(int i=0; i<VEC2.size(); i++)
    {std::vector<Point> V2=VEC2[i];
        pcl::PointCloud<pcl::PointXYZ>::Ptr In(new pcl::PointCloud<pcl::PointXYZ>);
        In->width    = V2.size();
        In->height   = 1;
        In->is_dense = false;

        for(int j=0; j<V2.size(); j++)
        {Lg::Point3f A;
            A =O +CGAL::to_double( V2[j].x())*m +CGAL::to_double( V2[j].y())*k;

            pcl::PointXYZ B=pcl::PointXYZ(A.x(),A.y(), A.z());
            In->points.push_back(B);

        }
        Hol.push_back(In);

    }

    Res.outer=out;
    Res.inners=Hol;
    return Res;
}   



////// 3D polygon projection on a plane

Polygon_with_holes  polygon_creation(  pcl::ModelCoefficients::Ptr pl, poly_holes& pol)
{   

    Polygon outer_p ;
    Lg::Point3f n(pl->values[0], pl->values[1], pl->values[2]);

    //n.Normalize();
    Lg::Point3f m = best_m(n);
    m.Normalize();
    Lg::Point3f k = n^m;
    k.Normalize();
    Lg::Point3f O=(-pl->values[3]/(n*n))*n;


    vector<Lg::Point2f> points_2d(pol.outer->points.size());
    std::vector<K::Point_2> points_cgal(points_2d.size());
    // TODO CALCULER O
    //Lg::Point3f O=(-pl->values[3]/((pl->values[0]*pl->values[0])+(pl->values[1]*pl->values[1])+(pl->values[2]*pl->values[2])))*n;



    for(int i = 0; i < pol.outer->points.size(); ++i)
    {//pcl::PointXYZ pt1=project(cloud->points[i],  pl);
        Lg::Point3f P=Lg::Point3f(pol.outer->points[i].x,pol.outer->points[i].y,pol.outer->points[i].z);

        Lg::Point3f P0=P-O;
        float xp=P0*m, yp=P0*k;
        outer_p.push_back(K::Point_2(xp, yp));

        


    }
    Polygon_with_holes poly;
    poly.outer_boundary()=outer_p ;

    for(int j=0; j<pol.inners.size(); j++)
    {Polygon Hol;
        for(int l=0; l<pol.inners[j]->points.size(); l++)
        {
            Lg::Point3f P=Lg::Point3f(pol.inners[j]->points[l].x,pol.inners[j]->points[l].y,pol.inners[j]->points[l].z);
            Lg::Point3f P0 = Lg::Point3f(P-O);
            float xp=P0*m, yp=P0*k;

            Hol.push_back(K::Point_2(xp, yp));


        }
        poly.add_hole( Hol ) ;
    }
    return poly;
}
