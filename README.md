# planar polygons based registration framework by Rahima Djahel(ENPC), Bruno Vallet(IGN) and  Pascal Monasse(ENPC)

Required dependencies: PCL, Eigen, CGAL,Boost

Poly_registration: an efficient algorithm to register an indoor and outdoor scan.

Test:



./Poly_registration ../data/indoor_scan.ply ../data/STR_Int.ply 30

where:

indoor_scan: the indoor scan
STR_Int: indoor points detected from exterior scans 
30 : distance threshold (can be adapted by the users)

to visualize the result:

cloudcompare.CloudCompare  TR_indoor_scan.ply  ../data/outdoor_scan.ply

Remark :for poly_registration

As the size of the 3D point cloud representing the indoor points seen from the outdoor is very small compared to the size of the indoor scan, we preferred to use two different thresholds of inliers (to detect planar polygons) as well as two different values for the minimum size of a planar region.

that's why we have fixed these two parameters in Poly_registration.cpp in order to be compatible with our data size. 



If you use our algorithm in any of your publications or projects, please cite our paper:

Djahel, Rahima, Bruno Vallet, and Pascal Monasse. "Towards Efficient Indoor/outdoor Registration Using Planar Polygons." ISPRS Annals of the Photogrammetry, Remote Sensing and Spatial Information Sciences 2 (2021): 51-58.


If you have any questions, you can send an email to :

rahima.djahel@enpc.fr

rdjahel@gmail.com
