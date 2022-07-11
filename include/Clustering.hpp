#include"Polygon.hpp"
/////

struct Cluster
{
    std::vector<Plane> Planes;
    Eigen::Vector3d clus_normal;
    string label="GENERAL";
};


Eigen::Vector3d wei_centroid(Cluster C)
{double nor_x=0;
    double nor_y=0;
    double nor_z=0;
    double s_w=0;
    for(int i=0; i<C.Planes.size(); i++)
    { double w=C.Planes[i].size/C.Planes[0].size;
        nor_x+=w*C.Planes[i].normal.x();
        nor_y+=w*C.Planes[i].normal.y();
        nor_z+=w*C.Planes[i].normal.z();
        s_w+=w;
    }
    Eigen::Vector3d w_normal;
    w_normal.x()=nor_x/s_w;
    w_normal.y()=nor_y/s_w;
    w_normal.z()=nor_z/s_w;
    w_normal.normalize();
    return w_normal;
}

//////////////////////////////////////

int k_min(std::vector<Eigen::Vector3d>& vec, Eigen::Vector3d N)
{
    int l=0;
    double min=10.0;
    for(int i=0; i<vec.size(); i++)
    {
        double val=1-std::abs(vec[i].dot(N));
        if(val<min)
        {min=val;
            l=i;
        }
    }
    return l;
}

/////////////////// Clusters generation

std::vector<Cluster> Cluster_generation(std::vector<Plane>& pl)
{std::vector<Cluster>res;
    Cluster C1,C2,C3;

    C1.Planes.push_back(pl[0]);
    pl.erase(pl.begin());

    int M=-1;
    int m=0;
    do{
        if((Angle(C1.Planes[0].normal, pl[m].normal)>85)&&(Angle(C1.Planes[0].normal, pl[m].normal)<95))
        {

            M=m;

        }
        m++;

    }
    while((M==-1)&&(m<pl.size()));
    if(M!=-1)
    {
        C2.Planes.push_back(pl[M]);
        pl.erase(pl.begin()+M);
    }
    int F=-1;
    int f=0;
    do{

        if((Angle(C1.Planes[0].normal, pl[f].normal)>85)&&(Angle(C2.Planes[0].normal, pl[f].normal)>85)&&(Angle(C1.Planes[0].normal, pl[f].normal)<95)&&(Angle(C2.Planes[0].normal, pl[f].normal)<95))
        {
            F=f;

        }
        f++;
    }
    while((F==-1)&&(f<pl.size()));
    if(F!=-1)
    {
        C3.Planes.push_back(pl[F]);
        pl.erase(pl.begin()+F);}
    bool d;
    //do{
    Eigen::Vector3d N1= wei_centroid(C1);
    Eigen::Vector3d N2= wei_centroid(C2);
    Eigen::Vector3d N3= wei_centroid(C3);
    std::vector<Eigen::Vector3d>vec;
    vec.push_back(N1);
    vec.push_back(N2);

    vec.push_back(N3);

    for(int i=0; i<pl.size(); i++)
    {
        if(pl[i].processed==false)
        {
            pl[i].processed=true;
            int k=k_min(vec,pl[i].normal);
            double val=1-std::abs(vec[k].dot(pl[i].normal));
            //if(val<0.003)
            if(val<0.09)
            {
                if(k==0)
                    C1.Planes.push_back(pl[i]);



                else if(k==1)
                    C2.Planes.push_back(pl[i]);



                else if(k==2)
                    C3.Planes.push_back(pl[i]);


            }
            else
                std::cout<<"we reject the plane"<<std::endl;
        }}



    if(C1.Planes.size()>0)
        res.push_back(C1);
    if(C2.Planes.size()>0)
        res.push_back(C2);
    if(C3.Planes.size()>0)
        res.push_back(C3);

    return res;
}

////////////// horizental cluster selection

void direction_clusters(vector<Cluster>&cl)
{
    Eigen::Vector3d X(1,0,0),Y(0,1,0), Z(0,0,1);
    vector<Eigen::Vector3d>V;
    V.push_back(X);
    V.push_back(Y);
    V.push_back(Z);
    for(int i=0; i<cl.size(); i++)
    {

        Eigen::Vector3d N=cl[i].Planes[0].normal;
        int id=k_min(V,N);
        if(id==0)
            cl[i].label="V_X";
        else if(id==1)
            cl[i].label="V_Y";
        else
            cl[i].label="HO";
    }
}
std::vector<Cluster> separ_clus(Cluster C)
{
    Cluster C1,C2;
    Eigen::Vector3d N=C.Planes[0].normal;
    for(int i=0; i<C.Planes.size(); i++)
    {
        if(Angle(N,C.Planes[i].normal)<20)
            C1.Planes.push_back(C.Planes[i]);
        else
            C2.Planes.push_back(C.Planes[i]);
    }
    std::vector<Cluster>res;
    if(C1.Planes.size()>0)
        res.push_back(C1);
    if(C2.Planes.size()>0)
        res.push_back(C2);
    return res;
}

///////////////////////////////////////////////////

pair<Eigen::Vector3d,Eigen::Vector3d> best_vec (vector<Cluster>&H1, vector<Cluster>&H2, Eigen::Vector3d V) 
{
    pair<Eigen::Vector3d,Eigen::Vector3d> vec;
    Eigen::Vector3d u1,u2;
    if((H1.size()>1)&&(H2.size()>1))
    {
        if((Angle(wei_centroid(H1[0]),V)<Angle(wei_centroid(H1[1]),V))&&(Angle(wei_centroid(H2[0]),V)<Angle(wei_centroid(H2[1]),V)))
        {

            u1=wei_centroid(H1[0]);
            u2=wei_centroid(H2[0]);

        }
        else if((Angle(wei_centroid(H1[0]),V)<Angle(wei_centroid(H1[1]),V))&&(Angle(wei_centroid(H2[1]),V)<Angle(wei_centroid(H2[0]),V)))
        {
            u1=wei_centroid(H1[0]);
            u2=wei_centroid(H2[1]);
        }
        else if((Angle(wei_centroid(H1[1]),V)<Angle(wei_centroid(H1[0]),V))&&(Angle(wei_centroid(H2[0]),V)<Angle(wei_centroid(H2[1]),V)))
        {
            u1=wei_centroid(H1[1]);
            u2=wei_centroid(H2[0]);
        }
        else if((Angle(wei_centroid(H1[1]),V)<Angle(wei_centroid(H1[0]),V))&&(Angle(wei_centroid(H2[1]),V)<Angle(wei_centroid(H2[0]),V)))
        {
            u1=wei_centroid(H1[1]);
            u2=wei_centroid(H2[1]);
        }
    }
    else if((H1.size()>1)&&(H2.size()==1))
    {



        if(Angle(wei_centroid(H1[0]),wei_centroid(H2[0]))<Angle(wei_centroid(H1[1]),wei_centroid(H2[0])))
        {

            u1=wei_centroid(H1[0]);
            u2=wei_centroid(H2[0]);

        }
        else
        {
            u1=wei_centroid(H1[1]);
            u2=wei_centroid(H2[0]);

        }}
    else if((H2.size()>1)&&(H1.size()==1))
    {



        if(Angle(wei_centroid(H2[0]),wei_centroid(H1[0]))<Angle(wei_centroid(H2[1]),wei_centroid(H1[0])))
        {

            u1=wei_centroid(H1[0]);
            u2=wei_centroid(H2[0]);

        }
        else
        {
            u1=wei_centroid(H1[0]);
            u2=wei_centroid(H2[1]);

        }}
    else
    { u1=wei_centroid(H1[0]);
        u2=wei_centroid(H2[0]);
    }
    vec=make_pair(u1,u2);

    return vec;
}
