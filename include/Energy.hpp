#include"Transformation.hpp"

#include <boost/geometry/geometry.hpp>
#include <boost/geometry/geometries/point_xy.hpp>
#include <boost/geometry/geometries/polygon.hpp>
#include <boost/geometry/geometries/adapted/c_array.hpp>
#include <boost/geometry/geometries/multi_polygon.hpp>
BOOST_GEOMETRY_REGISTER_C_ARRAY_CS(cs::cartesian)
using namespace boost::geometry;
typedef model::d2::point_xy<double> point_2d;
typedef model::polygon<point_2d> polygon_2d;
typedef model::box<point_2d> box_2d;

///// transformation of a plane

pcl::ModelCoefficients::Ptr Transform_plane(pcl::ModelCoefficients::Ptr& pl, Eigen::Isometry3d& T)
{
    pcl::ModelCoefficients::Ptr pl1(new pcl::ModelCoefficients);
    pcl::ModelCoefficients::Ptr pl2(new pcl::ModelCoefficients);
    Lg::Point3f n(pl->values[0], pl->values[1], pl->values[2]);


    Eigen::Vector3d V=Eigen::Vector3d(n.x(),n.y(),n.z());

    //V1=(trans.inverse()).transpose()*V;
    Eigen::Vector3d V1=T.linear()*V;
    Lg::Point3f O=( - pl->values[3]/(n*n))*n;

    Eigen::Vector4d V2=Eigen::Vector4d (O.x(), O.y(),O.z(),1);
    Eigen::Vector4d V3=T.matrix()*V2;
    Lg::Point3f O1=Lg::Point3f(V3.x(),V3.y(),V3.z());


    Lg::Point3f n1=Lg::Point3f(V1[0],V1[1],V1[2]);
    float nor=n1.Norm();
    n1.Normalize();

    float d1=-(n1*O1);



    pl1->values.resize (4);
    pl1->values[0]=n1.x();
    pl1->values[1]=n1.y();
    pl1->values[2]=n1.z();
    pl1->values[3]=d1;
    pl2=Normalize(pl1);
    return pl2;
}
/////////////////3D point cloud  transformation
pcl::PointCloud<pcl::PointXYZ>::Ptr Tr_cloud(pcl::PointCloud<pcl::PointXYZ>::Ptr cloud, Eigen::Isometry3d& T)
{
    pcl::PointCloud<pcl::PointXYZ>::Ptr cloud1(new pcl::PointCloud<pcl::PointXYZ>);
    pcl::PointCloud<pcl::PointXYZ>::Ptr cloud2(new pcl::PointCloud<pcl::PointXYZ>);
    Eigen::Vector4f centroid;
    pcl::compute3DCentroid(*cloud,centroid);
    pcl::PointXYZ C=pcl::PointXYZ(centroid[0],centroid[1],centroid[2]);
    for(int i=0; i<cloud->points.size(); i++)
    {
        pcl::PointXYZ P;
        P.x=cloud->points[i].x-C.x;
        P.y=cloud->points[i].y-C.y;
        P.z=cloud->points[i].z-C.z;
        cloud1->points.push_back(P);
    }
    pcl::transformPointCloud (*cloud1, *cloud1, T.matrix()
                              );

    for(int i=0; i<cloud1->points.size(); i++)
    {
        pcl::PointXYZ P;
        P.x=cloud1->points[i].x+C.x;
        P.y=cloud1->points[i].y+C.y;
        P.z=cloud1->points[i].z+C.z;
        cloud2->points.push_back(P);

    }
    return cloud2;
}

/////////////////3D polygon  transformation

poly_holes TR_hol(poly_holes& poly_h1,Eigen::Isometry3d& T)
{
    poly_holes trans_holes ;
    pcl::PointCloud<pcl::PointXYZ>::Ptr cloud1(new pcl::PointCloud<pcl::PointXYZ>);

    cloud1=Tr_cloud(poly_h1.outer,T);

    trans_holes.outer=cloud1;
    std::vector<pcl::PointCloud<pcl::PointXYZ>::Ptr> inners_clouds;
    for(int s=0; s<poly_h1.inners.size(); s++)
    {
        pcl::PointCloud<pcl::PointXYZ>::Ptr cloud3(new pcl::PointCloud<pcl::PointXYZ>);

        cloud3=Tr_cloud(poly_h1.inners[s],T);


        inners_clouds.push_back(cloud3);
    }


    trans_holes.inners=inners_clouds;
    return trans_holes;
}

///// quality criterion 


double quality_criterion (std::vector<std::pair<Cluster, Cluster>>&H,Eigen::Isometry3d T, double th)
{
    double summ=0;

    int start_s = clock();
    for(int q=0; q<H.size(); q++)
    {
        Cluster C1=H[q].first;
        Cluster C2=H[q].second;



        for(int i=0;i<C1.Planes.size(); i++)
        {


            for(int  j=0;j<C2.Planes.size(); j++){



                std::vector<poly_holes> poly_h1=C1.Planes[i].poly;
                std::vector<poly_holes> poly_h2=C2.Planes[j].poly;
                for(int h=0; h<poly_h1.size(); h++)
                {

                    poly_holes Pol=TR_hol( poly_h1[h], T);


                    for(int f=0; f<poly_h2.size(); f++)
                    {

                        pcl::ModelCoefficients::Ptr pl(new pcl::ModelCoefficients);
                        pl=Transform_plane(C1.Planes[i].Cof,T);
                        pcl::ModelCoefficients::Ptr besect=Besector_plane(pl,C2.Planes[j].Cof);
                        pcl::ModelCoefficients::Ptr  besect_N =Normalize(besect);
                        Eigen::Vector4f centroid1;
                        pcl::compute3DCentroid(*Pol.outer,centroid1);
                        pcl::PointXYZ pt1=pcl::PointXYZ(centroid1[0],centroid1[1],centroid1[2]);
                        Eigen::Vector4f centroid2;
                        pcl::compute3DCentroid(*poly_h2[f].outer,centroid2);
                        pcl::PointXYZ pt2=pcl::PointXYZ(centroid2[0],centroid2[1],centroid2[2]);
                        double D1=point_plane_distance(pt1,besect_N);
                        double D2=point_plane_distance(pt2,besect_N);
                        int s1=Pol.inners.size();
                        int s2=  poly_h2[f].inners.size();
                        float dist=D1+D2;
                        double totalArea=0.0;
                        double d1=max(0.0,((th*th)-(dist*dist)));



                        if(d1>0){

                            Polygon_with_holes  p=polygon_creation(besect_N,Pol);




                            Polygon_with_holes p1=polygon_creation(besect_N,poly_h2[f]);

                            //creat boost polygons from cgal polygons
                            polygon_2d poly1;
                            polygon_2d poly2;
                            std::vector<point_2d>  points;
                            Polygon q=p.outer_boundary();
                            for(int k=0; k<q.size();k++)
                            {double x=CGAL::to_double(q[k].x());

                                double y=CGAL::to_double(q[k].y());
                                points.push_back(point_2d(x,y));
                            }
                            boost::geometry::assign_points(poly1, points);
                            correct(poly1);
                            Polygon q1=p1.outer_boundary();
                            std::vector<point_2d> points1;
                            for(int k=0; k<q1.size();k++)
                            {

                                double x=CGAL::to_double(q1[k].x());

                                double y=CGAL::to_double(q1[k].y());
                                points1.push_back(point_2d(x,y));


                            }
                            boost::geometry::assign_points(poly2, points1);
                            correct(poly2);

                            poly1.inners().resize(s1);
                            ///////////////////////

                            for(typename Polygon_with_holes::Hole_const_iterator it_hole = p.holes_begin();
                                it_hole != p.holes_end(); it_hole++)
                            {
                                Polygon inner = *it_hole;
                                std::vector<point_2d> PP;
                                for(typename Polygon::Vertex_const_iterator it_vertex = inner.vertices_begin();
                                    it_vertex != inner.vertices_end(); it_vertex++){
                                    PP.push_back(point_2d( CGAL::to_double(it_vertex->x()),CGAL::to_double(it_vertex->y())) );

                                }


                                model::ring<point_2d> inn ;
                                assign_points(inn, PP);
                                poly1.inners().push_back(inn);
                            }
                            correct(poly1);
                            /////////
                            poly2.inners().resize(s2);
                            ///////////////////////

                            for(typename Polygon_with_holes::Hole_const_iterator it_hole = p1.holes_begin();
                                it_hole != p1.holes_end(); it_hole++)
                            {
                                Polygon inner = *it_hole;
                                std::vector<point_2d> PP;
                                for(typename Polygon::Vertex_const_iterator it_vertex = inner.vertices_begin();
                                    it_vertex != inner.vertices_end(); it_vertex++){
                                    PP.push_back(point_2d( CGAL::to_double(it_vertex->x()),CGAL::to_double(it_vertex->y())) );

                                }


                                model::ring<point_2d>inn ;
                                assign_points(inn, PP);
                                poly2.inners().push_back(inn);
                            }
                            correct(poly2);


                            ////////////////

                            typedef std::vector<polygon_2d> polygon_list;
                            polygon_list v;
                            intersection(poly1, poly2, v);
                            for (polygon_list::const_iterator it = v.begin(); it != v.end(); ++it)
                            {
                                totalArea+=boost::geometry::area(*it);
                            }}


                        summ+=totalArea*(d1/(th*th));

                    }}}
        }}

    int stop_s = clock();
    cout << "hypo TIME:  " << (stop_s-start_s)/double(CLOCKS_PER_SEC) << endl;
    return summ;
}     


