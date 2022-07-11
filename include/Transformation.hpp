#include"Clustering.hpp"



double angles(Eigen::Matrix3d R)
{

    double thetaX ;
    double thetaY ;
    double thetaZ;
    if ( R(0,2) < +1)
    {
        if ( R(0,2) > -1)
        {
            thetaY = asin (R(0,2) ) ;
            thetaX = atan2 (-R(1,2) , R(2,2) ) ;
            thetaZ = atan2 (-R(0,1) , R(0,0) ) ;
        }
        else
        {
            thetaY = -M_PI / 2 ;
            thetaX = -atan2( R(1,0) , R(1,1) ) ;
            thetaZ = 0;
        }
    }
    else
    {
        thetaY=+M_PI / 2 ;
        thetaX=atan2 ( R(1,0) , R(1,1) ) ;
        thetaZ=0;
    }

    double M=(thetaX+thetaY+thetaZ)/3;
    std::cout<<"mean="<<M<<std::endl;
    return M;
}

Eigen::Vector3d V_proj(Eigen::Vector3d& U1, Eigen::Vector3d& V1)
{
    double  v2_ls = U1.norm()*U1.norm();
    //std::cout<<"the values="<<v2_ls<<std::endl;
    Eigen::Vector3d Vec=V1-(U1*(V1.dot(U1)/v2_ls));
    Eigen::Vector3d Vec1=Vec.normalized();
    return Vec1;
}


Eigen::Matrix3d Rotation(Eigen::Vector3d& u0,Eigen::Vector3d& u1,Eigen::Vector3d& v0,Eigen::Vector3d& v1)
{
    Eigen::Vector3d v00=V_proj(u0,v0);


    Eigen::Vector3d s=u0.cross(v00);
    Eigen::Vector3d s0=s.normalized();
    Eigen::Vector3d v11=V_proj(u1,v1);

    Eigen::Vector3d q=u1.cross(v11);

    Eigen::Vector3d q0=q.normalized();
    Eigen::Matrix3d M1;
    M1.col(0)=u0;
    M1.col(1)=v00;
    M1.col(2)=s0;
    Eigen::Matrix3d M2;
    M2.col(0)=u1;
    M2.col(1)=v11;
    M2.col(2)=q0;
    Eigen::Matrix3d R=M2*M1.inverse();
    return R;
}






Eigen::Vector3d translation(vector<pair<Plane, Plane>>&vec, Eigen::Matrix3d& R)
{
    Eigen::Vector3d N1=R*vec[0].first.normal;
    Eigen::Vector3d N2=R*vec[1].first.normal;
    Eigen::Vector3d N3=R*vec[2].first.normal;
    //std::cout<<vec.size()<<std::endl;
    Eigen::Matrix3d M1;
    M1.row(0)=N1;

    M1.row(1)=N2;
    M1.row(2)=N3;

    Eigen::Vector3d D;
    D.x()=vec[0].first.d-vec[0].second.d;
    D.y()=vec[1].first.d-vec[1].second.d;
    D.z()=vec[2].first.d-vec[2].second.d;
    Eigen::Vector3d D1=D.normalized();

    Eigen::Vector3d trans=M1.inverse()*D;
    //cout<<trans<<endl;
    return trans;
}
