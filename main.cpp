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

#include <pcl/sample_consensus/ransac.h>
#include <pcl/sample_consensus/sac_model_registration.h>
#include <pcl/sample_consensus/sac_model_plane.h>
#include <pcl/sample_consensus/sac_model_sphere.h>
#include <pcl/sample_consensus/sac_model_circle3d.h>
#include <pcl/sample_consensus/sac_model_cone.h>
#include <pcl/sample_consensus/sac_model_cylinder.h>

#include <pcl/common/distances.h>
#include <pcl/search/kdtree.h>

#include <QString>
#include <QFile>
#include <QVector>
#include <QMatrix4x4>
#include <QMatrix3x3>
#include <QVector3D>
#include <QTime>
#include <QRandomGenerator>

#include "stl_reader.h"

#include "ap.h"
#include "dataanalysis.h"

#include "mathstatistics.h"
#include <gp_Pnt.hxx>
#include <gp_Trsf.hxx>
#include <gp_Ax1.hxx>
#include <gp_Ax3.hxx>
#include <gp_Pnt.hxx>
#include <gp_Vec.hxx>
#include <gp_Quaternion.hxx>

#include <random>

#define PI 3.14159265359

using namespace pcl;
using namespace pcl::io;
using namespace pcl::console;
using namespace pcl::registration;
////////////////////////////////////////////////////////////////////////////////

void TransformationPointCloud(PointCloud<PointXYZ>::Ptr &source, PointCloud<PointXYZ>::Ptr &target, gp_Trsf &transform);

struct point
{
    point(double x, double y, double z) : x(x), y(y), z(z) {}
    double x, y, z;
};

struct eulerAngles
{
    eulerAngles() : A(0), B(0), C(0) {}
    eulerAngles(double a, double b, double c) : A(a), B(b), C(c) {}
    double A, B, C;
};

////////////////////////////////////////////////////////////////////////////////
PointCloud<PointXYZ>::Ptr src, tgt, src_rot, tempCloud, src_king;
QVector<point> srcPoints, tgtPoints, rotPoints;

gp_Ax3 src_pca_ax, tgt_pca_ax;

///////////////////////////////////////////////////////////////////////////////////

void AddGausssianNoise (const PointCloud<PointXYZ>::ConstPtr &xyz_cloud, PointCloud<PointXYZ>::Ptr &xyz_cloud_filtered, float mean, float standard_deviation)
{
//  print_highlight ("Adding Gaussian noise with mean 0.0 and standard deviation %f\n", standard_deviation);

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

void RotatePoints(QVector<point> &sourcePoints, QVector<point> &rotatePoints, const float angle, gp_Vec &rotVec)
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

    gp_Trsf trsf;

    trsf.SetRotation(gp_Ax1(gp_Pnt(mean_x, mean_y, mean_z), rotVec), angle / 180 * PI);

    for (auto p: sourcePoints)
    {
        gp_Pnt gpPoint(p.x, p.y, p.z);

        gpPoint.Transform(trsf);

        x = gpPoint.X();
        y = gpPoint.Y();
        z = gpPoint.Z();

        rotatePoints.push_back(point(x, y, z));
    }
}

void RotatePoints(PointCloud<PointXYZ>::Ptr &sourcePoints, PointCloud<PointXYZ>::Ptr &rotatePoints, const float angle, const gp_Vec &rotVec)
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

    gp_Trsf trsf;

    trsf.SetRotation(gp_Ax1(gp_Pnt(mean_x, mean_y, mean_z), rotVec), angle / 180 * PI);

    TransformationPointCloud(sourcePoints, rotatePoints, trsf);

//    for (auto p: *sourcePoints)
//    {
//        gp_Pnt gpPoint(p.x, p.y, p.z);

//        gpPoint.Transform(trsf);

//        x = gpPoint.X();
//        y = gpPoint.Y();
//        z = gpPoint.Z();

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

//    qDebug() << "mean" << mean_x << mean_y << mean_z;
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

//    qDebug() << "Pca test time :" << time.elapsed()/1000.0;

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

//    ax_pca.SetLocation(gp_Pnt(mean_x, mean_y, mean_z));
//    ax_pca.SetXDirection(gp_Dir(basis1_x, basis1_y, basis1_z));
//    ax_pca.SetYDirection(gp_Dir(basis2_x, basis2_y, basis2_z));
//    ax_pca.SetDirection(gp_Dir(basis0_x, basis0_y, basis0_z));

    gp_Ax3 temp_ax(gp_Pnt(mean_x, mean_y, mean_z),gp_Dir(basis0_x, basis0_y, basis0_z),gp_Dir(basis1_x, basis1_y, basis1_z));
    ax_pca = temp_ax;
}

void PcaTest(PointCloud<PointXYZ>::Ptr &points, gp_Pnt &mean_pt ,QVector<gp_Dir> &dir)
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

//    qDebug() << "mean" << mean_x << mean_y << mean_z;
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

//    qDebug() << "Pca test time :" << time.elapsed()/1000.0;

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

    mean_pt.SetXYZ(gp_XYZ(mean_x, mean_y, mean_z));

    dir.clear();
    dir.push_back(gp_Dir(basis0_x, basis0_y, basis0_z));
    dir.push_back(gp_Dir(basis1_x, basis1_y, basis1_z));
    dir.push_back(gp_Dir(basis2_x, basis2_y, basis2_z));
}

void vec2euler(QMatrix3x3* ptrMatrix, eulerAngles &e)
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
    e.A *= 180/PI;
    e.B *= 180/PI;
    e.C *= 180/PI;
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

void TransformationPointCloud(PointCloud<PointXYZ>::Ptr &source, PointCloud<PointXYZ>::Ptr &target, gp_Trsf &transform)
{
    target->clear();

    for(auto p: *source)
    {
        gp_Pnt point(p.x, p.y, p.z);
        point.Transform(transform);

        target->push_back(PointXYZ(point.X(), point.Y(), point.Z()));
    }
}

void RansacPoints(PointCloud<PointXYZ>::Ptr &source, PointCloud<PointXYZ>::Ptr &target)
{
    PointCloud<PointXYZ>::Ptr temp(new PointCloud<PointXYZ>());

    QTime time;

    time.start();
    std::vector<int> inliers;
    bool random = false;

    pcl::SampleConsensusModelRegistration<pcl::PointXYZ>::Ptr model_r(new pcl::SampleConsensusModelRegistration<pcl::PointXYZ> (source, random));
    pcl::SampleConsensusModelSphere<pcl::PointXYZ>::Ptr model_s(new pcl::SampleConsensusModelSphere<pcl::PointXYZ> (source, random));
    pcl::SampleConsensusModelPlane<pcl::PointXYZ>::Ptr model_p(new pcl::SampleConsensusModelPlane<pcl::PointXYZ> (source, random));

    pcl::RandomSampleConsensus<pcl::PointXYZ> ransac (model_p);
    ransac.setDistanceThreshold (10);
    ransac.computeModel();
    ransac.getInliers(inliers);

    pcl::copyPointCloud (*source, inliers, *temp);

    pcl::SampleConsensusModelCircle3D<pcl::PointXYZ>::Ptr model_circ(new pcl::SampleConsensusModelCircle3D<pcl::PointXYZ> (temp, random));

    pcl::RandomSampleConsensus<pcl::PointXYZ> ransac2 (model_circ);
    ransac2.setDistanceThreshold (5);
    ransac2.computeModel();
    ransac2.getInliers(inliers);

    pcl::copyPointCloud (*temp, inliers, *target);

}

double KdtreeTest(PointCloud<PointXYZ>::Ptr &source, PointCloud<PointXYZ>::Ptr &target)
{
    pcl::search::KdTree<PointXYZ> kdtree;
    kdtree.setInputCloud (target);

    unsigned int no_of_neighbors = 1;

    std::vector<int> k_indices;
    k_indices.resize (no_of_neighbors);
    std::vector<float> k_distances;
    k_distances.resize (no_of_neighbors);
    QVector<double> distances;

    for (const auto &point : *source)
    {
        kdtree.nearestKSearch (point, no_of_neighbors, k_indices, k_distances);
//        qDebug() << "*******************************************";
//        qDebug() << k_indices;
//        qDebug() << k_distances;
        distances.push_back(k_distances.at(0));

//        qDebug() << gp_Pnt(point.x,point.y,point.z).Distance(gp_Pnt((*target)[k_indices.at(0)].x,(*target)[k_indices.at(0)].y,(*target)[k_indices.at(0)].z));
    }

    return MathStatistics::getMean(distances) + MathStatistics::getVariance(distances) * 3;
}

gp_Trsf GetBestTransform(PointCloud<PointXYZ>::Ptr &Source, PointCloud<PointXYZ>::Ptr &Target)
{
    PointCloud<PointXYZ>::Ptr temp(new PointCloud<PointXYZ>());
    PointCloud<PointXYZ>::Ptr source_rot(new PointCloud<PointXYZ>());
    PointCloud<PointXYZ>::Ptr ransac_cloud(new PointCloud<PointXYZ>());
    gp_Trsf transform, trsf_result;
    gp_Ax3 source_ax, target_ax;
    PcaTest(Source, source_ax);
    PcaTest(Target, target_ax);

    transform.SetDisplacement(source_ax, target_ax);
    TransformationPointCloud(Source, source_rot, transform);
    RansacPoints(source_rot, ransac_cloud);

//    double min_error = KdtreeTest(ransac_cloud, Target);
    double min_error = KdtreeTest(source_rot, Target);
    trsf_result = transform;

    for(int i = 0; i < 360; i += 90)
    {
        gp_Ax3 temp_ax = target_ax;
        gp_Trsf x_rot;
        x_rot.SetRotation(gp_Ax1(temp_ax.Location(), temp_ax.XDirection()), i * PI / 180);
        temp_ax.Transform(x_rot);
        for(int j = 0; j < 360; j += 90)
        {
            gp_Trsf y_rot;
            y_rot.SetRotation(gp_Ax1(temp_ax.Location(), temp_ax.YDirection()), j * PI / 180);

            gp_Trsf temp_trsf;
            temp_trsf.PreMultiply(x_rot);
            temp_trsf.PreMultiply(y_rot);

//            TransformationPointCloud(ransac_cloud, temp, temp_trsf);
            TransformationPointCloud(source_rot, temp, temp_trsf);

            double error = KdtreeTest(temp, Target);

            gp_Trsf t;
            t.PreMultiply(transform);
            t.PreMultiply(temp_trsf);

            if(error < min_error)
            {
                error = min_error;
                trsf_result = t;
            }
        }
    }

    for(int i = -90; i <= 90; i += 180)
//    for(int i = 0; i < 360; i += 90)
    {
        gp_Ax3 temp_ax = target_ax;
        gp_Trsf z_rot;
        z_rot.SetRotation(gp_Ax1(temp_ax.Location(), temp_ax.Direction()), i * PI / 180);
        temp_ax.Transform(z_rot);
        for(int j = 0; j < 360; j += 90)
        {
            gp_Trsf y_rot;
            y_rot.SetRotation(gp_Ax1(temp_ax.Location(), temp_ax.YDirection()), j * PI / 180);

            gp_Trsf temp_trsf;
            temp_trsf.PreMultiply(z_rot);
            temp_trsf.PreMultiply(y_rot);

//            TransformationPointCloud(ransac_cloud, temp, temp_trsf);
            TransformationPointCloud(source_rot, temp, temp_trsf);

            double error = KdtreeTest(temp, Target);

            gp_Trsf t;
            t.PreMultiply(transform);
            t.PreMultiply(temp_trsf);

            if(error < min_error)
            {
                error = min_error;
                trsf_result = t;
            }
        }
    }

    return trsf_result;
}

gp_Trsf Evren(PointCloud<PointXYZ>::Ptr &source, PointCloud<PointXYZ>::Ptr target)
{
    PointCloud<PointXYZ>::Ptr source_rot;
    source_rot.reset(new PointCloud<PointXYZ>());
    gp_Pnt source_mean_pt;
    QVector<gp_Dir> source_vec;
    QVector<gp_Ax3> source_ax_vec;
    gp_Ax3 target_ax;

    gp_Trsf result;
    double min_error = KdtreeTest(source, target);

    PcaTest(source, source_mean_pt, source_vec);
    PcaTest(target, target_ax);

    for (int i = 0; i < 3; ++i)
    {
        for (int j = 0; j < 3; ++j)
        {
            if(i == j)
                continue;
            gp_Ax3 ax(source_mean_pt,source_vec.at(i),source_vec.at(j));
            for (int k = 0; k < 8; ++k)
            {
                gp_Ax3 temp_ax = ax;
                int temp = k;
                if(temp % 2 == 1)
                    temp_ax.XReverse();
                if(temp % 2 == 1)
                    temp_ax.YReverse();
                temp >>= 1;
                if(temp % 2 == 1)
                    temp_ax.ZReverse();
                source_ax_vec.push_back(temp_ax);
            }
        }
    }


    for(auto source_ax: source_ax_vec)
    {
        gp_Trsf trsf;
        trsf.SetDisplacement(source_ax, target_ax);

//        gp_Ax3 temp_ax = source_ax;
//        temp_ax.Transform(trsf);

//        if(!temp_ax.IsCoplanar(target_ax, 0.01, 0.01))
//        {
//            qDebug() << "temp ax " << temp_ax.Location().X() << temp_ax.Location().Y() << temp_ax.Location().Z();
//            qDebug() << "target ax " << target_ax.Location().X() << target_ax.Location().Y() << target_ax.Location().Z();
//        }


        TransformationPointCloud(source, source_rot, trsf);
        double error = KdtreeTest(source_rot,target);

        if(error < min_error)
        {
            min_error = error;
            result = trsf;
        }

    }
    return result;
}

void TestProgram()
{
    QString bunnyPath = "E:\\Codes\\Dataset\\Bunny.ply";
    PLYReader reader;

    double mean_x = 0;
    double mean_y = 0;
    double mean_z = 0;

    int count_ok = 0, count_fault = 0;

    double test_ok_value = 1;

    PointCloud<PointXYZ>::Ptr source, target, source_result, temp;
    source.reset (new PointCloud<PointXYZ>);
    target.reset (new PointCloud<PointXYZ>);
    source_result.reset (new PointCloud<PointXYZ>);
    temp.reset (new PointCloud<PointXYZ>);

    gp_Trsf trsf_scale;
    trsf_scale.SetScaleFactor(1000);

    reader.read(bunnyPath.toStdString(), *temp);

    TransformationPointCloud(temp, source, trsf_scale);

    savePLYFileBinary("source.ply", *source);

    for (auto p: *source)
    {
        mean_x += p.x;
        mean_y += p.y;
        mean_z += p.z;
    }

    mean_x /= source->size();
    mean_y /= source->size();
    mean_z /= source->size();

    for (int i = 0; i < 10; ++i)
    {
        int x = QRandomGenerator::global()->bounded(-100, 100);
        int y = QRandomGenerator::global()->bounded(-100, 100);
        int z = QRandomGenerator::global()->bounded(-100, 100);

        qDebug() << "*******************************************";

        gp_Vec vec(x, y, z);
        vec.Normalize();

        gp_Ax1 ax(gp_Pnt(mean_x, mean_y, mean_z), vec);

        for(double i = -250; i <= 250; i += 12.5)
        {
            gp_Trsf trsf_rot;
            trsf_rot.SetRotation(ax, i);
            TransformationPointCloud(source, target, trsf_rot);
//            RotatePoints(source, target, i, vec);
//            RotatePoints(source, temp, i, vec);
//            AddGausssianNoise(temp, target, 0, 1.5);
//            gp_Trsf transform = GetBestTransform(source, target);
            gp_Trsf transform = Evren(source, target);
            TransformationPointCloud(source, source_result, transform);

            int error = KdtreeTest(source_result, target);
            if(error > test_ok_value)
            {
                qDebug() << "----------------------------------------------";
                qDebug() << QString("X = %1, Y = %2, Z = %3, Angle = %4, Error =").arg(x).arg(y).arg(z).arg(i) << error;
                savePLYFileBinary(QString("_%1_%2_%3_%4_source_result.ply").arg(x).arg(y).arg(z).arg(i).toStdString(), *source_result);
                savePLYFileBinary(QString("_%1_%2_%3_%4_target.ply").arg(x).arg(y).arg(z).arg(i).toStdString(), *target);

                auto q_rot = trsf_rot.GetRotation();
                auto q_trsf = transform.GetRotation();

                gp_Trsf trsf2rot;
                trsf2rot.PreMultiply(transform.Inverted());
                trsf2rot.PreMultiply(trsf_rot);

                auto q_2 = trsf2rot.GetRotation();

                double A_rot, B_rot, C_rot, A_trsf, B_trsf, C_trsf, A_2, B_2, C_2;

                q_rot.GetEulerAngles(gp_Intrinsic_ZYX, A_rot, B_rot, C_rot);
                q_trsf.GetEulerAngles(gp_Intrinsic_ZYX, A_trsf, B_trsf, C_trsf);
                q_2.GetEulerAngles(gp_Intrinsic_ZYX, A_2, B_2, C_2);

                A_rot *= 180 / PI;
                B_rot *= 180 / PI;
                C_rot *= 180 / PI;

                B_trsf *= 180 / PI;
                A_trsf *= 180 / PI;
                C_trsf *= 180 / PI;

                A_2 *= 180 / PI;
                B_2 *= 180 / PI;
                C_2 *= 180 / PI;

                qDebug() << "Rotation" << trsf_rot.TranslationPart().X() << trsf_rot.TranslationPart().Y() << trsf_rot.TranslationPart().Z() << A_rot << B_rot << C_rot;
                qDebug() << "Transform" << transform.TranslationPart().X() << transform.TranslationPart().Y() << transform.TranslationPart().Z() << A_trsf << B_trsf << C_trsf;
                qDebug() << "trsf2rot" << trsf2rot.TranslationPart().X() << trsf2rot.TranslationPart().Y() << trsf2rot.TranslationPart().Z() << A_2 << B_2 << C_2;

                TransformationPointCloud(source_result, temp, trsf2rot);

                qDebug() << "Error" << KdtreeTest(temp, target);

                qDebug() << "----------------------------------------------";
                count_fault++;
            }
            else
                count_ok++;
        }

        for(double i = -180; i <= 180; i += 30)
        {
            gp_Trsf trsf_rot;
            trsf_rot.SetRotation(ax, i);
            TransformationPointCloud(source, target, trsf_rot);
//            RotatePoints(source, target, i, vec);
//            RotatePoints(source, temp, i, vec);
//            AddGausssianNoise(temp, target, 0, 1.5);
//            gp_Trsf transform = GetBestTransform(source, target);
            gp_Trsf transform = Evren(source, target);
            TransformationPointCloud(source, source_result, transform);

            int error = KdtreeTest(source_result, target);
            if(error > test_ok_value)
            {
                qDebug() << "----------------------------------------------";
                qDebug() << QString("X = %1, Y = %2, Z = %3, Angle = %4, Error =").arg(x).arg(y).arg(z).arg(i) << error;
                savePLYFileBinary(QString("_%1_%2_%3_%4_source_result.ply").arg(x).arg(y).arg(z).arg(i).toStdString(), *source_result);
                savePLYFileBinary(QString("_%1_%2_%3_%4_target.ply").arg(x).arg(y).arg(z).arg(i).toStdString(), *target);

                auto q_rot = trsf_rot.GetRotation();
                auto q_trsf = transform.GetRotation();

                gp_Trsf trsf2rot;
                trsf2rot.PreMultiply(transform.Inverted());
                trsf2rot.PreMultiply(trsf_rot);

                auto q_2 = trsf2rot.GetRotation();

                double A_rot, B_rot, C_rot, A_trsf, B_trsf, C_trsf, A_2, B_2, C_2;

                q_rot.GetEulerAngles(gp_Intrinsic_ZYX, A_rot, B_rot, C_rot);
                q_trsf.GetEulerAngles(gp_Intrinsic_ZYX, A_trsf, B_trsf, C_trsf);
                q_2.GetEulerAngles(gp_Intrinsic_ZYX, A_2, B_2, C_2);

                A_rot *= 180 / PI;
                B_rot *= 180 / PI;
                C_rot *= 180 / PI;

                B_trsf *= 180 / PI;
                A_trsf *= 180 / PI;
                C_trsf *= 180 / PI;

                A_2 *= 180 / PI;
                B_2 *= 180 / PI;
                C_2 *= 180 / PI;

                qDebug() << "Rotation" << trsf_rot.TranslationPart().X() << trsf_rot.TranslationPart().Y() << trsf_rot.TranslationPart().Z() << A_rot << B_rot << C_rot;
                qDebug() << "Transform" << transform.TranslationPart().X() << transform.TranslationPart().Y() << transform.TranslationPart().Z() << A_trsf << B_trsf << C_trsf;
                qDebug() << "trsf2rot" << trsf2rot.TranslationPart().X() << trsf2rot.TranslationPart().Y() << trsf2rot.TranslationPart().Z() << A_2 << B_2 << C_2;

                TransformationPointCloud(source_result, temp, trsf2rot);

                qDebug() << "Error" << KdtreeTest(temp, target);

                qDebug() << "----------------------------------------------";
                count_fault++;
            }
            else
                count_ok++;
        }
    }

    qDebug() << "Başarılı" << count_ok;
    qDebug() << "Başarısız" << count_fault;



}

void init()
{
    QString bunnyPath = "E:\\Codes\\Dataset\\Bunny.ply";
    PLYReader reader;

    src.reset (new PointCloud<PointXYZ>);
    tgt.reset (new PointCloud<PointXYZ>);
    src_rot.reset (new PointCloud<PointXYZ>);
    tempCloud.reset (new PointCloud<PointXYZ>);
    src_king.reset (new PointCloud<PointXYZ>);

    gp_Trsf tscale;
    tscale.SetScaleFactor(1000);

    reader.read(bunnyPath.toStdString(), *tempCloud);

    TransformationPointCloud(tempCloud, src, tscale);

    float angle = 330;
    gp_Vec vec(1,-2,3);

//    float angle = 180;
//    gp_Vec vec(1,0,0);

    RotatePoints(src, tgt, angle, vec);

    int mean_x = 0;
    int mean_y = 0;
    int mean_z = 0;

    for (auto p: *src)
    {
        mean_x += p.x;
        mean_y += p.y;
        mean_z += p.z;
    }

    mean_x /= src->size();
    mean_y /= src->size();
    mean_z /= src->size();

    gp_Trsf t1;

    t1.SetRotation(gp_Ax1(gp_Pnt(mean_x,mean_y,mean_z), vec), angle / 180 * PI);

    *tempCloud = *tgt;

    AddGausssianNoise(tempCloud, tgt, 0.0, 2);

    savePLYFileBinary("src.ply", *src);
    savePLYFileBinary("tgt.ply", *tgt);

    // Pca Test
    PcaTest(src, src_pca_ax);
    PcaTest(tgt, tgt_pca_ax);

    gp_Trsf trsf;
    trsf.SetDisplacement(src_pca_ax, tgt_pca_ax);
    TransformationPointCloud(src, src_rot, trsf);

    savePLYFileBinary("src_rot.ply", *src_rot);

    gp_Trsf trsf_king = GetBestTransform(src, tgt);
    TransformationPointCloud(src, src_king, trsf_king);

    gp_Quaternion q1 = trsf.GetRotation();
    gp_Quaternion q2 = t1.GetRotation();
    gp_Quaternion q3 = trsf_king.GetRotation();

    double A1,B1,C1;
    double A2,B2,C2;
    double A3,B3,C3;

    q1.GetEulerAngles(gp_Intrinsic_ZYX, A1, B1, C1);
    q2.GetEulerAngles(gp_Intrinsic_ZYX, A2, B2, C2);
    q3.GetEulerAngles(gp_Intrinsic_ZYX, A3, B3, C3);

    A1 *= 180 / PI;
    B1 *= 180 / PI;
    C1 *= 180 / PI;

    A2 *= 180 / PI;
    B2 *= 180 / PI;
    C2 *= 180 / PI;

    A3 *= 180 / PI;
    B3 *= 180 / PI;
    C3 *= 180 / PI;

    qDebug() << "PCA" << trsf.TranslationPart().X() << trsf.TranslationPart().Y() << trsf.TranslationPart().Z() << A1 << B1 << C1;
    qDebug() << "Target" << t1.TranslationPart().X() << t1.TranslationPart().Y() << t1.TranslationPart().Z() << A2 << B2 << C2;
    qDebug() << "trsf_king" << trsf_king.TranslationPart().X() << trsf_king.TranslationPart().Y() << trsf_king.TranslationPart().Z() << A3 << B3 << C3;

    savePLYFileBinary("src_king.ply", *src_king);

    qDebug() << "src_pca_ax.Angle(tgt_pca_ax) : " << src_pca_ax.Angle(tgt_pca_ax) * 180 / PI;

    qDebug() << "Source Error";
    CalculateError(src, tgt);

    qDebug() << "Source Rot Error";
    CalculateError(src_rot, tgt);

    qDebug() << "Source King Error";
    CalculateError(src_king, tgt);

    *src = *src_king;

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
//  init();
    TestProgram();
    return 0;
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

  qDebug() << "Transformasyon hesaplama süresi : " << time.elapsed() / 1000.0;

  std::cerr << transform << std::endl;
  // Transform the data and write it to disk
  PointCloud<PointXYZ> output;

  transformPointCloud (*src, output, transform);

  qDebug() << time.elapsed() / 1000.0;

  savePCDFileBinary ("source_transformed.pcd", output);
  savePLYFileBinary ("source_transformed.ply", output);
}
