#include "Energy.hpp"

int main (int argc, char** argv)
{
    cout << "Usage: " << argv[0] << " file.ply file2.ply   d_thr" << endl;
    if(argc<3) return 1;

    int m = 1;
    string file="", file2="";
    if (argc > m)
        file = argv[m++];

    if (argc > m)
        file2 = argv[m++];

    

    double inlier_threshold1 = 0.04;
    double inlier_threshold2 = 0.02;
    int nb_iter=1000;
    int min_plane_size1 = 1000;
    int min_plane_size2 = 100;
    double alpha=0.05;
    double d_thr;
    
    if (argc > m) d_thr = atof(argv[m++]);
    pcl::PointCloud<pcl::PointXYZ>::Ptr cloud1(new pcl::PointCloud<pcl::PointXYZ>);
    pcl::PointCloud<pcl::PointXYZ>::Ptr cloud2(new pcl::PointCloud<pcl::PointXYZ>);
    pcl::PointCloud<pcl::PointXYZ>::Ptr cloud3(new pcl::PointCloud<pcl::PointXYZ>);

    if (pcl::io::loadPLYFile<pcl::PointXYZ>(file, *cloud1) == -1) //* load the file
    {
        //PCL_ERROR("Couldn't read file " + file + "\n");
        return (-1);
    }

    if (pcl::io::loadPLYFile<pcl::PointXYZ>(file2, *cloud2) == -1) //* load the file
    {
        //PCL_ERROR("Couldn't read file " + file + "\n");
        return (-1);
    }
    


    Eigen::Isometry3d TR= Eigen::Isometry3d::Identity() ;
    double max_Ene=0;

    /////////////////


    srand (time(NULL));
    int start_s = clock();
    vector<vector<poly_holes>>P_H1;
    vector<vector<poly_holes>>P_H2;
    vector<Plane>Planes1;
    vector<Plane>Planes2;
    vector<Cluster> clu1,clu2;
    int S1=0,S2=0;
    Planes_data result1,result2;
    do{
        pcl::PointCloud<pcl::PointXYZ>::Ptr cloudCopy(new pcl::PointCloud<pcl::PointXYZ>);
        pcl::PointCloud<pcl::PointXYZ>::Ptr cloud2Copy(new pcl::PointCloud<pcl::PointXYZ>);

        pcl::copyPointCloud<pcl::PointXYZ,pcl::PointXYZ>(*cloud1, * cloudCopy);
        pcl::copyPointCloud<pcl::PointXYZ,pcl::PointXYZ>(*cloud2, * cloud2Copy);





        if(S1<3)
        {
            result1 = plane_detection(cloudCopy,inlier_threshold1, nb_iter,min_plane_size1,alpha);
            for(int i=0;i<result1.Pl.size(); i++)

            {vector<poly_holes> p_h;

                list<Polygon_with_holes> res=Alpha_shape(  result1.Pl[i], result1.In_cloud[i],alpha);

                vector<pcl::PointCloud<pcl::PointXYZ>::Ptr>Coll;
                int k=0;
                for(typename list<Polygon_with_holes>::iterator it_ring = res.begin();
                    it_ring != res.end(); it_ring++){

                    Polygon_with_holes PP=(*it_ring);

                    pcl::PointCloud<pcl::PointXYZ>::Ptr cll(new pcl::PointCloud<pcl::PointXYZ>);

                    poly_holes Res=Polygon3D_creation(PP , result1.Pl[i]);
                    p_h.push_back(Res);
                }
                P_H1.push_back(p_h);
            }


            for(int i=0;i<result1.Pl.size(); i++)
            {
                vector<poly_holes> poly_h1=P_H1[i];

                Plane pl;
                Eigen::Vector3d N(result1.Pl[i]->values[0],result1.Pl[i]->values[1],result1.Pl[i]->values[2]);
                pl.normal=N;
                pl.d=result1.Pl[i]->values[3];
                pl.Cloud=result1.In_cloud[i];
                pl.size=result1.In_cloud[i]->points.size();
                pl.poly=poly_h1;
                pl.Cof=result1.Pl[i];
                Planes1.push_back(pl);

            }
            //// sorting planes according to inliers number
            for(int i=0;i<Planes1.size() - 1;i++)
            {
                for(int j=i+1;j<Planes1.size();j++)
                {
                    if(Planes1[j].size> Planes1[i].size)
                    {
                        Plane temp=Planes1[j];
                        Planes1[j]= Planes1[i];
                        Planes1[i]=temp;
                    }
                }
            }

            clu1=Cluster_generation(Planes1);
            S1=clu1.size();

        }

        if(S2<3)
        {     result2 = plane_detection(cloud2Copy,inlier_threshold2, nb_iter,min_plane_size2,alpha);
            for(int i=0;i<result2.Pl.size(); i++)

            {vector<poly_holes> p_h;

                list<Polygon_with_holes> res =Alpha_shape(  result2.Pl[i], result2.In_cloud[i],alpha);

                vector<pcl::PointCloud<pcl::PointXYZ>::Ptr>Coll;
                int k=0;
                for(typename list<Polygon_with_holes>::iterator it_ring = res.begin();
                    it_ring != res.end(); it_ring++){

                    Polygon_with_holes PP=(*it_ring);

                    pcl::PointCloud<pcl::PointXYZ>::Ptr cll(new pcl::PointCloud<pcl::PointXYZ>);

                    poly_holes Res=Polygon3D_creation(PP , result2.Pl[i]);
                    p_h.push_back(Res);
                }
                P_H2.push_back(p_h);
            }

            //////////////////////////////////////////////////////////////


            //////////

            for(int i=0;i<result2.Pl.size(); i++)
            {
                vector<poly_holes> poly_h2=P_H2[i];

                Plane pl;
                Eigen::Vector3d N(result2.Pl[i]->values[0],result2.Pl[i]->values[1],result2.Pl[i]->values[2]);
                pl.normal=N;
                pl.d=result2.Pl[i]->values[3];
                pl.Cloud=result2.In_cloud[i];
                pl.size=result2.In_cloud[i]->points.size();
                pl.poly=poly_h2;
                pl.Cof=result2.Pl[i];
                Planes2.push_back(pl);
            }
            //// sorting planes according to inliers number

            for(int i=0;i<Planes2.size() - 1;i++)
            {
                for(int j=i+1;j<Planes2.size();j++)
                {
                    if(Planes2[j].size> Planes2[i].size)
                    {
                        Plane temp=Planes2[j];
                        Planes2[j]= Planes2[i];
                        Planes2[i]=temp;
                    }
                }
            }

            //// Planes clustering



            clu2=Cluster_generation(Planes2);
            S2=clu2.size();
        }
        if(S1<3)
        {std::cout<<"give me inlier_threshold for scan1"<<std::endl;
            cin >> inlier_threshold1;
            std::cout<<"give me min_plane_size for scan1"<<std::endl;
            cin >> min_plane_size1;
        }
        else
            if(S2<3)
            {std::cout<<"give me inlier_threshold for scan2"<<std::endl;
                cin >> inlier_threshold2;
                std::cout<<"give me min_plane_size for scan2"<<std::endl;
                cin >> min_plane_size2;
            }

    }
    while((S1<3)||(S2<3));
    std::cout<<"size1="<<result1.Pl.size()<<std::endl;
    std::cout<<"size2="<<result2.Pl.size()<<std::endl;
    Cluster H1,H2,V1_x,V2_x,V1_y,V2_y;

    direction_clusters(clu1);
    direction_clusters(clu2);
    for(int i=0; i<clu1.size(); i++)
    {
        if(clu1[i].label=="HO")
            H1=clu1[i];
        else if(clu1[i].label=="V_X")
            V1_x=clu1[i];
        else
            V1_y=clu1[i];
    }
    for(int i=0; i<clu2.size(); i++)
    {
        if(clu2[i].label=="HO")
            H2=clu2[i];
        else if(clu2[i].label=="V_X")
            V2_x=clu2[i];
        else
            V2_y=clu2[i];
    }


    vector<Cluster> hh1=separ_clus(H1);
    vector<Cluster> hh2=separ_clus(H2);


    Eigen::Vector3d v1,v2,h1,h2,Z(0,0,1),X(1,0,0),Y(0,1,0);
    pair<Eigen::Vector3d,Eigen::Vector3d> vec1= best_vec(hh1,hh2,Z);
    v1=vec1.first;
    v2=vec1.second;
    vector<Cluster> vv1_1=separ_clus(V1_x);
    vector<Cluster> vv2_1=separ_clus(V2_x);
    pair<Eigen::Vector3d,Eigen::Vector3d> vec2= best_vec(vv1_1,vv2_1,X);

    vector<Cluster> vv1_2=separ_clus(V1_y);
    vector<Cluster> vv2_2=separ_clus(V2_y);

    pair<Eigen::Vector3d,Eigen::Vector3d> vec3= best_vec(vv1_2,vv2_2,Y);
    if(Angle(vec2.first,vec2.second)< Angle(vec3.first,vec3.second))
    {
        h1=vec2.first;
        h2=vec2.second;
    }
    else
    {
        h1=vec3.first;
        h2=vec3.second;
    }

    Eigen::Matrix3d R= Eigen::Matrix3d::Identity();
    R=Rotation(v1,v2,h1,h2);
    std::vector<std::pair<Cluster, Cluster>>H;
    H.push_back(std::make_pair(V1_x,V2_x));
    H.push_back(std::make_pair(V1_y,V2_y));
    H.push_back(std::make_pair(H1,H2));
    /////////////////////// RANSAC /////////////////////////////////

    for(int i=0; i<20000; i++)
    {
        int a=rand() % H1.Planes.size();
        int b=rand() % H2.Planes.size();
        int c=rand() % V1_y.Planes.size();
        int d=rand() % V2_y.Planes.size();
        int e=rand() % V1_x.Planes.size();
        int f=rand() % V2_x.Planes.size();
        Plane P1,P2,P3,P4,P5,P6;
        P1=H1.Planes[a]; P2=H2.Planes[b];P3=V1_y.Planes[c];P4=V2_y.Planes[d];P5=V1_x.Planes[e];P6=V2_x.Planes[f];

        vector<pair<Plane, Plane>>hyp;

        hyp.push_back(make_pair(P5,P6));
        hyp.push_back(make_pair(P3,P4));
        hyp.push_back(make_pair(P1,P2));
        Eigen::Vector3d v1,v2,h1,h2;



        Eigen::Vector3d t=translation(hyp, R);
        Eigen::Isometry3d M = Eigen::Isometry3d::Identity();
        M.linear()=R;
        M.translation()=t;
        double Ene=quality_criterion(H,M,d_thr);

        if(Ene>max_Ene)
        {
            max_Ene=Ene;
            TR.matrix()=M.matrix();
        }
    }

    std::cout<<"max_Ene="<<max_Ene<<std::endl;
    std::cout<<"transformation"<<std::endl;
    std::cout<< TR.matrix()<<std::endl;
    cloud3=Tr_cloud(cloud1,TR);

    

    pcl::io::savePLYFile("TR_indoor_scan.ply",*cloud3);
    int stop_s = clock();
    cout << "EXECUTION TIME:  " << (double)(stop_s-start_s)/(CLOCKS_PER_SEC) << endl;
    
    

    
    angles(TR.linear());


    //////////////////////////////////////////////////////////
    return 0;

}

