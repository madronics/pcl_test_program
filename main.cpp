#include <pcl/console/parse.h>
#include <pcl/point_types.h>
#include <pcl/point_cloud.h>
#include <pcl/point_representation.h>

#include <pcl/io/pcd_io.h>
#include <pcl/io/ply_io.h>
#include <pcl/conversions.h>
#include <pcl/filters/uniform_sampling.h>
#include <pcl/features/normal_3d.h>
#include <pcl/features/fpfh.h>
#include <pcl/registration/correspondence_estimation.h>
#include <pcl/registration/correspondence_rejection_distance.h>
#include <pcl/registration/transformation_estimation_svd.h>

#include <QString>
#include <QFile>
#include <QVector>
#include <QMatrix4x4>
#include <QMatrix3x3>
#include <QVector3D>
#include <QTime>

#include "stl_reader.h"

#include "ap.h"
#include "dataanalysis.h"

#include "mathstatistics.h"

#include <gp_Trsf.hxx>
#include <gp_Ax1.hxx>
#include <gp_Ax3.hxx>
#include <gp_Pnt.hxx>
#include <gp_Vec.hxx>

#include <random>

#define PI 3.14159265359

using namespace pcl;
using namespace pcl::io;
using namespace pcl::console;
using namespace pcl::registration;
PointCloud<PointXYZ>::Ptr src, tgt, rot, temp;
gp_Ax3 src_ax, tgt_ax, src_pca, tgt_pca;

////////////////////////////////////////////////////////////////////////////////

struct point
{
    point(double x, double y, double z) : x(x), y(y), z(z) {}
    double x, y, z;
};

QVector<point> srcPoints, tgtPoints, rotPoints;

struct eulerAngles
{
    eulerAngles() : A(0), B(0), C(0) {}
    eulerAngles(double a, double b, double c) : A(a), B(b), C(c) {}
    double A, B, C;
};

////////////////////////////////////////////////////////////////////////////////
eulerAngles srcAngles, tgtAngles, rotAngles, outAngles;

QMatrix3x3 *srcMatrix = nullptr;
QMatrix3x3 *tgtMatrix = nullptr;
QMatrix3x3 *rotMatrix = nullptr;
QMatrix3x3 *outMatrix = nullptr;

QVector3D vecA(0, 0, 1);
QVector3D vecB(0, 1, 0);
QVector3D vecC(1, 0, 0);
///////////////////////////////////////////////////////////////////////////////////

void AddGausssianNoise (const PointCloud<PointXYZ>::ConstPtr &xyz_cloud, PointCloud<PointXYZ>::Ptr &xyz_cloud_filtered, float mean, float standard_deviation)
{
  print_highlight ("Adding Gaussian noise with mean 0.0 and standard deviation %f\n", standard_deviation);

  xyz_cloud_filtered->points.resize (xyz_cloud->size ());
  xyz_cloud_filtered->header = xyz_cloud->header;
  xyz_cloud_filtered->width = xyz_cloud->width;
  xyz_cloud_filtered->height = xyz_cloud->height;


  std::random_device rd;
  std::mt19937 rng(rd());
  std::normal_distribution<float> nd (mean, standard_deviation);

  for (std::size_t point_i = 0; point_i < xyz_cloud->size (); ++point_i)
  {
    (*xyz_cloud_filtered)[point_i].x = (*xyz_cloud)[point_i].x + nd (rng);
    (*xyz_cloud_filtered)[point_i].y = (*xyz_cloud)[point_i].y + nd (rng);
    (*xyz_cloud_filtered)[point_i].z = (*xyz_cloud)[point_i].z + nd (rng);
  }
}


///////////////////////////////////////////////////////////////////////////////////
QString FileToText(QString filename)
{
    QString data;
    QFile f(filename);
    if(f.open(QIODevice::ReadOnly))
    {
        data = f.readAll();
        f.close();
    }
    return data;
}

QString DeleteExcessSpace(QString s)
{
    int index = 0;
    s = s.trimmed();

    while (true)
    {
        index = s.indexOf("  ",index);
        if (index == -1) {
            break;
        }
        s.replace(index, 2, " ");
    }
    return s;
}

void Simplification(QVector<point> &points)
{
    qDebug() << "Simplification" << points.size();

    qSort(points.begin(),points.end(),[](point p1, point p2){
        if(p1.x > p2.x)
            return false;
        else if(p1.x < p2.x)
            return true;
        if(p1.y > p2.y)
            return false;
        else if(p1.y < p2.y)
            return true;
        if(p1.z > p2.z)
            return false;
        else if(p1.z < p2.z)
            return true;
        else
            return false;
    });

    QVector<point> pointsSimp;

    point lastPoint = *points.begin();
    pointsSimp.push_back(lastPoint);

    for (auto it = points.begin(); it != points.end();it++)
    {
        if(!(it->x == lastPoint.x && it->y == lastPoint.y && it->z == lastPoint.z))
        {
            pointsSimp.push_back(*it);
        }
        lastPoint = *it;
    }

    points = pointsSimp;
    qDebug() << "Simplification fin " << points.size();
}

void LoadFromCsvFile(QString filename, QVector<point> &points)
{
    QString data = FileToText(filename);
    QStringList lineList = data.split("\n");

    points.clear();

    for (auto l:lineList)
    {
        QStringList wordList = l.split(',');
        if(wordList.size() < 3)
            continue;
        double x = wordList.at(0).trimmed().toDouble();
        double y = wordList.at(1).trimmed().toDouble();
        double z = wordList.at(2).trimmed().toDouble();
        points.push_back(point(x, y, z));
    }
}

void LoadFromStlFile(QString filename, QVector<point> &points)
{
    points.clear();

    try {
        stl_reader::StlMesh <float, unsigned int> mesh (filename.toStdString());

        qDebug() << mesh.num_tris();

        for(size_t itri = 0; itri < mesh.num_tris(); ++itri) {
            for(size_t icorner = 0; icorner < 3; ++icorner) {
                const float* c = mesh.tri_corner_coords (itri, icorner);
                points.push_back(point(c[0], c[1], c[2]));
            }
        }
    }
    catch (std::exception& e) {
        qDebug() << e.what();
    }

    Simplification(points);
}

void RotatePoints(QVector<point> &sourcePoints, QVector<point> &rotatePoints, const float angle, QVector3D &rotVec)
{
    double x, y, z;
    double mean_x = 0;
    double mean_y = 0;
    double mean_z = 0;
    rotatePoints.clear();

    for (auto p: sourcePoints)
    {
        mean_x += p.x;
        mean_y += p.y;
        mean_z += p.z;
    }

    mean_x /= sourcePoints.length();
    mean_y /= sourcePoints.length();
    mean_z /= sourcePoints.length();

    qDebug() << "mean" << mean_x << mean_y << mean_z;

    gp_Trsf trsf;

    trsf.SetRotation(gp_Ax1(gp_Pnt(mean_x, mean_y, mean_z),gp_Vec(rotVec.x(),rotVec.y(),rotVec.z())), angle / 180 * PI);

    for (auto p: sourcePoints)
    {
        gp_Pnt gpPoint(p.x, p.y, p.z);

        gpPoint.Transform(trsf);

        x = gpPoint.X();
        y = gpPoint.Y();
        z = gpPoint.Z();

        rotatePoints.push_back(point(x, y, z));
    }

//    QMatrix4x4 trsf;

//    trsf.translate(mean_x, mean_y, mean_z);
//    trsf.rotate(angle, rotVec);
//    trsf.translate(-mean_x, -mean_y, -mean_z);

//    for (auto p: sourcePoints)
//    {
//        QMatrix4x4 mat(0, 0, 0, p.x,
//                       0, 0, 0, p.y,
//                       0, 0, 0, p.z,
//                       0, 0, 0, 1);
//        mat = trsf * mat;

//        x = mat(0, 3);
//        y = mat(1, 3);
//        z = mat(2, 3);

//        rotatePoints.push_back(point(x, y, z));
//    }
}

void RotatePoints(PointCloud<PointXYZ>::Ptr &sourcePoints, PointCloud<PointXYZ>::Ptr &rotatePoints, const float angle, QVector3D &rotVec)
{
    double x, y, z;
    double mean_x = 0;
    double mean_y = 0;
    double mean_z = 0;
    rotatePoints->clear();

    for (auto p: *sourcePoints)
    {
        mean_x += p.x;
        mean_y += p.y;
        mean_z += p.z;
    }

    mean_x /= sourcePoints->size();
    mean_y /= sourcePoints->size();
    mean_z /= sourcePoints->size();

    qDebug() << "mean" << mean_x << mean_y << mean_z;

    gp_Trsf trsf;

    trsf.SetRotation(gp_Ax1(gp_Pnt(mean_x, mean_y, mean_z),gp_Vec(rotVec.x(),rotVec.y(),rotVec.z())), angle / 180 * PI);

    for (auto p: *sourcePoints)
    {
        gp_Pnt gpPoint(p.x, p.y, p.z);

        gpPoint.Transform(trsf);

        x = gpPoint.X();
        y = gpPoint.Y();
        z = gpPoint.Z();

        rotatePoints->push_back(PointXYZ(x, y, z));
    }

//    QMatrix4x4 trsf;

//    trsf.translate(mean_x, mean_y, mean_z);
//    trsf.rotate(angle, rotVec);
//    trsf.translate(-mean_x, -mean_y, -mean_z);

//    for (auto p: *sourcePoints)
//    {
//        QMatrix4x4 mat(0, 0, 0, p.x,
//                       0, 0, 0, p.y,
//                       0, 0, 0, p.z,
//                       0, 0, 0, 1);
//        mat = trsf * mat;

//        x = mat(0, 3);
//        y = mat(1, 3);
//        z = mat(2, 3);

//        rotatePoints->push_back(PointXYZ(x, y, z));
//    }
}



void Points2CloudPtr(QVector<point> &points, PointCloud<PointXYZ>::Ptr &cloud)
{
    cloud->clear();

    for(auto p: points)
    {
        cloud->push_back(PointXYZ(p.x, p.y, p.z));
    }
}

void PcaTest(PointCloud<PointXYZ>::Ptr &points, QMatrix3x3 **matrix)
{
    alglib::real_2d_array ptInput;
    QVector<double> p;

    for (auto i: *points)
    {
        p.append(i.x);
        p.append(i.y);
        p.append(i.z);
    }

    QTime time;
    time.start();

    ptInput.setcontent(points->size(), 3, p.data());
    alglib::ae_int_t info;
    alglib::real_1d_array eigValues;
    alglib::real_2d_array eigVectors;

    alglib::pcabuildbasis(ptInput, points->size(), 3, info, eigValues, eigVectors);

    qDebug() << "Pca test time :" << time.elapsed()/1000.0;

    // now the vectors can be accessed as follows:

    double basis0_x = eigVectors[0][0];
    double basis0_y = eigVectors[1][0];
    double basis0_z = eigVectors[2][0];

    double basis1_x = eigVectors[0][1];
    double basis1_y = eigVectors[1][1];
    double basis1_z = eigVectors[2][1];

    double basis2_x = eigVectors[0][2];
    double basis2_y = eigVectors[1][2];
    double basis2_z = eigVectors[2][2];

    float vectors [] = { basis0_x , basis0_y , basis0_z,
                         basis1_x , basis1_y , basis1_z,
                         basis2_x , basis2_y , basis2_z};

    if(*matrix != nullptr)
        delete *matrix;

    *matrix = new QMatrix3x3(vectors);

    qDebug() << **matrix;
}

void PcaTest(PointCloud<PointXYZ>::Ptr &points, gp_Ax3 &ax_pca)
{
    double mean_x = 0;
    double mean_y = 0;
    double mean_z = 0;

    for (auto p: *points)
    {
        mean_x += p.x;
        mean_y += p.y;
        mean_z += p.z;
    }

    mean_x /= points->size();
    mean_y /= points->size();
    mean_z /= points->size();

    qDebug() << "mean" << mean_x << mean_y << mean_z;
    alglib::real_2d_array ptInput;
    QVector<double> p;

    for (auto i: *points)
    {
        p.append(i.x);
        p.append(i.y);
        p.append(i.z);
    }

    QTime time;
    time.start();

    ptInput.setcontent(points->size(), 3, p.data());
    alglib::ae_int_t info;
    alglib::real_1d_array eigValues;
    alglib::real_2d_array eigVectors;

    alglib::pcabuildbasis(ptInput, points->size(), 3, info, eigValues, eigVectors);

    qDebug() << "Pca test time :" << time.elapsed()/1000.0;

    // now the vectors can be accessed as follows:

    double basis0_x = eigVectors[0][0];
    double basis0_y = eigVectors[1][0];
    double basis0_z = eigVectors[2][0];

    double basis1_x = eigVectors[0][1];
    double basis1_y = eigVectors[1][1];
    double basis1_z = eigVectors[2][1];

    double basis2_x = eigVectors[0][2];
    double basis2_y = eigVectors[1][2];
    double basis2_z = eigVectors[2][2];

    ax_pca.SetLocation(gp_Pnt(mean_x, mean_y, mean_z));
    ax_pca.SetDirection(gp_Dir(basis2_x, basis2_y, basis2_z));
    ax_pca.SetXDirection(gp_Dir(basis0_x, basis0_y, basis0_z));
    ax_pca.SetYDirection(gp_Dir(basis1_x, basis1_y, basis1_z));
}

void vec2euler(QMatrix3x3* ptrMatrix, eulerAngles &e)
{
    try
    {
        auto &m = *ptrMatrix;
        if (m(0, 2) < +1)
        {
         if (m(0, 2) > -1)
         {
          e.B = asin(m(0, 2));
          e.C = atan2(-m(1, 2),m(2, 2));
          e.A = atan2(-m(0, 1),m(0, 0));
         }
         else // r02 = -1
         {
          // Not a unique solution: thetaZ - thetaX = atan2(r10,r11)
          e.B = -PI/2;
          e.C = -atan2(m(1, 0),m(1, 1));
          e.A = 0;
         }
        }
        else // r02 = +1
        {
         // Not a unique solution: thetaZ + thetaX = atan2(r10,r11)
         e.B = +PI/2;
         e.C = atan2(m(1, 0),m(1, 1));
         e.A = 0;
        }
    }

    catch(...)
    {

    }

    e.A *= 180/PI;
    e.B *= 180/PI;
    e.C *= 180/PI;

    qDebug() << QString("A: %1, B: %2, C: %3").arg(e.A).arg(e.B).arg(e.C);
}

void CalculateError(PointCloud<PointXYZ>::Ptr& src, PointCloud<PointXYZ>::Ptr& tgt)
{
    if(src->size() != tgt->size())
        return;

    QVector<double> distances;

    auto size = src->size();
    for(int i = 0; i < size; i++)
    {
        PointXYZ *srcPoint = src->data() + i;
        PointXYZ *tgtPoint = tgt->data() + i;
        double distance = sqrt(pow(srcPoint->x - tgtPoint->x, 2) + pow(srcPoint->y - tgtPoint->y, 2) + pow(srcPoint->z - tgtPoint->z, 2));
        distances.push_back(distance);
    }

    qDebug() << "***********************************************************";
    qDebug() << "Mean : " << MathStatistics::getMean(distances);
    qDebug() << "Max : " << MathStatistics::max(distances);
    qDebug() << "Median : " << MathStatistics::median(distances);
    qDebug() << "Standard Deviation : " << MathStatistics::getStdDev(distances);
    qDebug() << "***********************************************************";
}

void init()
{
    QString srcPath = "E:\\Codes\\Dataset\\face.csv";
    QString tgtPath = "E:\\Codes\\Dataset\\face_distorted.csv";

    QString stlPath = "E:\\Codes\\Dataset\\Gyroid.stl";

    QString bunnyPath = "E:\\Codes\\Dataset\\Bunny.ply";

    PLYReader reader;

    src.reset (new PointCloud<PointXYZ>);
    tgt.reset (new PointCloud<PointXYZ>);
    rot.reset (new PointCloud<PointXYZ>);
    temp.reset (new PointCloud<PointXYZ>);

    reader.read(bunnyPath.toStdString(), *src);

    float angle = 180;
    QVector3D vec(1, 1, 0);

    RotatePoints(src, tgt, angle, vec);

    *temp = *tgt;

    AddGausssianNoise(temp, tgt, 0.0, 0.004);

    savePLYFileBinary("src.ply", *src);
    savePLYFileBinary("tgt.ply", *tgt);

////    LoadFromCsvFile(srcPath, srcPoints);
////    LoadFromCsvFile(tgtPath, tgtPoints);
////    qDebug() << "Load From CSV OK";
////    qDebug() << "Source Lenght " << srcPoints.length();
////    qDebug() << "Target Lenght " << tgtPoints.length();
//    LoadFromStlFile(stlPath, srcPoints);
//    RotatePoints(srcPoints, rotPoints, angle, vec);
//    qDebug() << "Rotate OK";
//    Points2CloudPtr(srcPoints, src);
////    Points2CloudPtr(tgtPoints, tgt);
//    Points2CloudPtr(rotPoints, tgt);
//    qDebug() << "CloudPtr OK";

    qDebug() << "test src";
    PcaTest(src, &srcMatrix);
    qDebug() << "test tgt";
    PcaTest(tgt, &tgtMatrix);

    // Pca Test
    PcaTest(src, src_pca);
    PcaTest(tgt, tgt_pca);

    qDebug() << "Angle : " << src_pca.Angle(tgt_pca) * 180 / PI;

    vec2euler(srcMatrix, srcAngles);
    vec2euler(tgtMatrix, tgtAngles);

    qDebug() << "Source Targer Error";
    CalculateError(src, tgt);

    RotatePoints(src, rot, -srcAngles.B, vecB);
    *temp = *rot;
    RotatePoints(temp, rot, -srcAngles.C, vecC);
    //    *temp = *rot;
    //    RotatePoints(temp, rot, -srcAngles.C, vecC);
    PcaTest(rot, &rotMatrix);
    *src = *rot;

    RotatePoints(tgt, rot, -tgtAngles.B, vecB);
    *temp = *rot;
    RotatePoints(temp, rot, -tgtAngles.C, vecC);
//    *temp = *rot;
//    RotatePoints(temp, rot, -tgtAngles.C, vecC);
    PcaTest(rot, &rotMatrix);
    *tgt = *rot;

    RotatePoints(tgt, rot, srcAngles.C, vecC);
    *temp = *rot;
    RotatePoints(temp, rot, srcAngles.B, vecB);
//    *temp = *rot;
//    RotatePoints(temp, rot, srcAngles.C, vecC);
    PcaTest(rot, &rotMatrix);
    *temp = *rot;

    savePLYFileBinary("temp.ply", *temp);

    qDebug() << "Source Targer Error Çevirdikten Sonra";
    CalculateError(src, tgt);
}

////////////////////////////////////////////////////////////////////////////////
void
estimateKeypoints (const PointCloud<PointXYZ>::Ptr &src,
                   const PointCloud<PointXYZ>::Ptr &tgt,
                   PointCloud<PointXYZ> &keypoints_src,
                   PointCloud<PointXYZ> &keypoints_tgt)
{
  // Get an uniform grid of keypoints
  UniformSampling<PointXYZ> uniform;
  uniform.setRadiusSearch (1);  // 1m

  uniform.setInputCloud (src);
  uniform.filter (keypoints_src);

  uniform.setInputCloud (tgt);
  uniform.filter (keypoints_tgt);

  // For debugging purposes only: uncomment the lines below and use pcl_viewer to view the results, i.e.:
  // pcl_viewer source_pcd keypoints_src.pcd -ps 1 -ps 10
  savePCDFileBinary ("keypoints_src.pcd", keypoints_src);
  savePCDFileBinary ("keypoints_tgt.pcd", keypoints_tgt);
}

////////////////////////////////////////////////////////////////////////////////
void
estimateNormals (const PointCloud<PointXYZ>::Ptr &src,
                 const PointCloud<PointXYZ>::Ptr &tgt,
                 PointCloud<Normal> &normals_src,
                 PointCloud<Normal> &normals_tgt)
{
  NormalEstimation<PointXYZ, Normal> normal_est;
  normal_est.setInputCloud (src);
  normal_est.setRadiusSearch (0.5);  // 50cm
  normal_est.compute (normals_src);

  normal_est.setInputCloud (tgt);
  normal_est.compute (normals_tgt);

  // For debugging purposes only: uncomment the lines below and use pcl_viewer to view the results, i.e.:
  // pcl_viewer normals_src.pcd
  PointCloud<PointNormal> s, t;
  copyPointCloud (*src, s);
  copyPointCloud (normals_src, s);
  copyPointCloud (*tgt, t);
  copyPointCloud (normals_tgt, t);
  savePCDFileBinary ("normals_src.pcd", s);
  savePCDFileBinary ("normals_tgt.pcd", t);
}

////////////////////////////////////////////////////////////////////////////////
void
estimateFPFH (const PointCloud<PointXYZ>::Ptr &src,
              const PointCloud<PointXYZ>::Ptr &tgt,
              const PointCloud<Normal>::Ptr &normals_src,
              const PointCloud<Normal>::Ptr &normals_tgt,
              const PointCloud<PointXYZ>::Ptr &keypoints_src,
              const PointCloud<PointXYZ>::Ptr &keypoints_tgt,
              PointCloud<FPFHSignature33> &fpfhs_src,
              PointCloud<FPFHSignature33> &fpfhs_tgt)
{
  FPFHEstimation<PointXYZ, Normal, FPFHSignature33> fpfh_est;
  fpfh_est.setInputCloud (keypoints_src);
  fpfh_est.setInputNormals (normals_src);
  fpfh_est.setRadiusSearch (1); // 1m
  fpfh_est.setSearchSurface (src);
  fpfh_est.compute (fpfhs_src);

  fpfh_est.setInputCloud (keypoints_tgt);
  fpfh_est.setInputNormals (normals_tgt);
  fpfh_est.setSearchSurface (tgt);
  fpfh_est.compute (fpfhs_tgt);

  // For debugging purposes only: uncomment the lines below and use pcl_viewer to view the results, i.e.:
  // pcl_viewer fpfhs_src.pcd
  PCLPointCloud2 s, t, out;
  toPCLPointCloud2 (*keypoints_src, s); toPCLPointCloud2 (fpfhs_src, t); concatenateFields (s, t, out);
  savePCDFile ("fpfhs_src.pcd", out);
  toPCLPointCloud2 (*keypoints_tgt, s); toPCLPointCloud2 (fpfhs_tgt, t); concatenateFields (s, t, out);
  savePCDFile ("fpfhs_tgt.pcd", out);
}

////////////////////////////////////////////////////////////////////////////////
void
findCorrespondences (const PointCloud<FPFHSignature33>::Ptr &fpfhs_src,
                     const PointCloud<FPFHSignature33>::Ptr &fpfhs_tgt,
                     Correspondences &all_correspondences)
{
  CorrespondenceEstimation<FPFHSignature33, FPFHSignature33> est;
  est.setInputCloud (fpfhs_src);
  est.setInputTarget (fpfhs_tgt);
  est.determineReciprocalCorrespondences (all_correspondences);
}

////////////////////////////////////////////////////////////////////////////////
void
rejectBadCorrespondences (const CorrespondencesPtr &all_correspondences,
                          const PointCloud<PointXYZ>::Ptr &keypoints_src,
                          const PointCloud<PointXYZ>::Ptr &keypoints_tgt,
                          Correspondences &remaining_correspondences)
{
  CorrespondenceRejectorDistance rej;
  rej.setInputSource<PointXYZ> (keypoints_src);
  rej.setInputTarget<PointXYZ> (keypoints_tgt);
  rej.setMaximumDistance (1);    // 1m
  rej.setInputCorrespondences (all_correspondences);
  rej.getCorrespondences (remaining_correspondences);
}


////////////////////////////////////////////////////////////////////////////////
void
computeTransformation (const PointCloud<PointXYZ>::Ptr &src,
                       const PointCloud<PointXYZ>::Ptr &tgt,
                       Eigen::Matrix4f &transform)
{
  // Get an uniform grid of keypoints
  PointCloud<PointXYZ>::Ptr keypoints_src (new PointCloud<PointXYZ>),
                            keypoints_tgt (new PointCloud<PointXYZ>);

  estimateKeypoints (src, tgt, *keypoints_src, *keypoints_tgt);
  print_info ("Found %zu and %zu keypoints for the source and target datasets.\n", static_cast<std::size_t>(keypoints_src->size ()), static_cast<std::size_t>(keypoints_tgt->size ()));

  // Compute normals for all points keypoint
  PointCloud<Normal>::Ptr normals_src (new PointCloud<Normal>),
                          normals_tgt (new PointCloud<Normal>);
  estimateNormals (src, tgt, *normals_src, *normals_tgt);
  print_info ("Estimated %zu and %zu normals for the source and target datasets.\n", static_cast<std::size_t>(normals_src->size ()), static_cast<std::size_t>(normals_tgt->size ()));

  // Compute FPFH features at each keypoint
  PointCloud<FPFHSignature33>::Ptr fpfhs_src (new PointCloud<FPFHSignature33>),
                                   fpfhs_tgt (new PointCloud<FPFHSignature33>);
  estimateFPFH (src, tgt, normals_src, normals_tgt, keypoints_src, keypoints_tgt, *fpfhs_src, *fpfhs_tgt);

  // Copy the data and save it to disk
/*  PointCloud<PointNormal> s, t;
  copyPointCloud (*keypoints_src, s);
  copyPointCloud (normals_src, s);
  copyPointCloud (*keypoints_tgt, t);
  copyPointCloud (normals_tgt, t);*/

  // Find correspondences between keypoints in FPFH space
  CorrespondencesPtr all_correspondences (new Correspondences),
                     good_correspondences (new Correspondences);
  findCorrespondences (fpfhs_src, fpfhs_tgt, *all_correspondences);

  // Reject correspondences based on their XYZ distance
  rejectBadCorrespondences (all_correspondences, keypoints_src, keypoints_tgt, *good_correspondences);

  for (int i = 0; i < good_correspondences->size (); ++i)
    std::cerr << good_correspondences->at (i) << std::endl;
  // Obtain the best transformation between the two sets of keypoints given the remaining correspondences
  TransformationEstimationSVD<PointXYZ, PointXYZ> trans_est;
  trans_est.estimateRigidTransformation (*keypoints_src, *keypoints_tgt, *good_correspondences, transform);
}

/* ---[ */
int
main (int argc, char** argv)
{
  init();
  // Parse the command line arguments for .pcd files
//  std::vector<int> p_file_indices;
//  p_file_indices = parse_file_extension_argument (argc, argv, ".pcd");
//  if (p_file_indices.size () != 2)
//  {
//    print_error ("Need one input source PCD file and one input target PCD file to continue.\n");
//    print_error ("Example: %s source.pcd target.pcd\n", argv[0]);
//    return (-1);
//  }

//  // Load the files
//  print_info ("Loading %s as source and %s as target...\n", argv[p_file_indices[0]], argv[p_file_indices[1]]);
//  src.reset (new PointCloud<PointXYZ>);
//  tgt.reset (new PointCloud<PointXYZ>);
//  if (loadPCDFile (argv[p_file_indices[0]], *src) == -1 || loadPCDFile (argv[p_file_indices[1]], *tgt) == -1)
//  {
//    print_error ("Error reading the input files!\n");
//    return (-1);
//  }

  // Compute the best transformtion
  QTime time;
  time.start();
  Eigen::Matrix4f transform;
  computeTransformation (src, tgt, transform);

  Eigen::Matrix4f rotate;

  qDebug() << "Transformasyon hesaplama süresi : " << time.elapsed() / 1000.0;

  std::cerr << transform << std::endl;
  // Transform the data and write it to disk
  PointCloud<PointXYZ> output;
  PointCloud<PointXYZ>::Ptr out;
  out.reset(new PointCloud<PointXYZ>);

  transformPointCloud (*src, output, transform);
  transformPointCloud (*src, *out, transform);
  savePLYFileBinary("source_transformed_pca.ply", output);

  qDebug() << time.elapsed() / 1000.0;

  qDebug() << srcAngles.A << srcAngles.B << srcAngles.C;

  RotatePoints(out, rot, srcAngles.C, vecC);
  *temp = *rot;
  RotatePoints(temp, rot, srcAngles.B, vecB);
  *out = *rot;

  savePCDFileBinary ("source_transformed.pcd", output);
  savePLYFileBinary("source_transformed.ply", output);

  savePLYFileBinary("src_pca.ply", *src);
  savePLYFileBinary("tgt_pca.ply", *tgt);
}
