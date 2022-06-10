//
// Created by zpk on 2021/3/5.
//
#include "NewEdge.h"
#include "Thirdparty/g2o/g2o/core/factory.h"
#include "Thirdparty/g2o/g2o/stuff/macros.h"
using namespace std;
namespace ORB_SLAM2
{

    bool EdgeGroundProjectOnlyPose::read(std::istream& is)
    {
        is>>_measurement;

        is>>information()(0,0);

        return true;
    }

    bool EdgeGroundProjectOnlyPose::write(std::ostream& os) const
    {
        os<<_measurement<<" ";

        os<<" "<<information()(0,0);

        return os.good();
    }

    void EdgeGroundProjectOnlyPose::linearizeOplus() {
        VertexSE3Expmap *vi = static_cast<VertexSE3Expmap *>(_vertices[0]);

        Vector3d xyz = vi->estimate().inverse().map(Po);
        Matrix3d xyz_trans = Matrix3d::Zero();
        xyz_trans(0, 1) = -xyz(2, 0);
        xyz_trans(0, 2) = xyz(1, 0);
        xyz_trans(1, 0) = xyz(2, 0);
        xyz_trans(1, 2) = -xyz(0, 0);
        xyz_trans(2, 0) = -xyz(1, 0);
        xyz_trans(2, 1) = xyz(0, 0);


        Matrix<double, 4, 6> POM;
        POM.block(0, 0, 3, 3) = Matrix3d::Identity();
        POM.block(0, 3, 3, 3) = -xyz_trans;
        POM.row(3) = Vector6d::Zero().transpose();

        Vector4d M;
        M << a, b, c, d;

        _jacobianOplusXi = M.transpose() * POM;
    }




    bool EdgeGroundProjectRelPose::read(std::istream& is)
    {
        is>>_measurement;

        is>>information()(0,0);

        return true;
    }

    bool EdgeGroundProjectRelPose::write(std::ostream& os) const
    {
        os<<_measurement<<" ";

        os<<" "<<information()(0,0);

        return os.good();
    }

    void EdgeGroundProjectRelPose::linearizeOplus() {
        VertexSE3Expmap *vi = static_cast<VertexSE3Expmap *>(_vertices[0]);
        VertexSE3Expmap *vj = static_cast<VertexSE3Expmap *>(_vertices[1]);

        Vector3d xyzj = vj->estimate().inverse().map(Po);
        Vector3d xyzi = vi->estimate().map(xyzj);
        Matrix3d xyz_transi = Matrix3d::Zero();
        Matrix3d xyz_transj = Matrix3d::Zero();

        xyz_transi(0, 1) = -xyzi(2, 0);
        xyz_transi(0, 2) = xyzi(1, 0);
        xyz_transi(1, 0) = xyzi(2, 0);
        xyz_transi(1, 2) = -xyzi(0, 0);
        xyz_transi(2, 0) = -xyzi(1, 0);
        xyz_transi(2, 1) = xyzi(0, 0);


        xyz_transj(0, 1) = -xyzj(2, 0);
        xyz_transj(0, 2) = xyzj(1, 0);
        xyz_transj(1, 0) = xyzj(2, 0);
        xyz_transj(1, 2) = -xyzj(0, 0);
        xyz_transj(2, 0) = -xyzj(1, 0);
        xyz_transj(2, 1) = xyzj(0, 0);




        Matrix<double, 4, 6> pomi = Matrix<double, 4, 6>::Zero();
        Matrix<double, 4, 6> pomj = Matrix<double, 4, 6>::Zero();
        pomi.block(0, 0, 3, 3) = Matrix3d::Identity();
        pomi.block(0, 3, 3, 3) = -xyz_transi;
        pomj.block(0, 0, 3, 3) = Matrix3d::Identity();
        pomj.block(0, 3, 3, 3) = -xyz_transj;


        Vector4d M;
        M << a, b, c, d;

        _jacobianOplusXi = -M.transpose() * pomi;
        _jacobianOplusXj = M.transpose() * pomj;
    }




    void EdgeCoplanarPrior::linearizeOplus() {
        const VertexGroundPlane* v1 = static_cast<const VertexGroundPlane*>(_vertices[0]);
        const VertexSBAPointXYZ* v2 = static_cast<const VertexSBAPointXYZ*>(_vertices[1]);
        const Vector3d G = v1->estimate();
        const Vector3d P = v2->estimate();

        _jacobianOplusXi = P;
        _jacobianOplusXj = G;

        }

    void EdgeGroundProjectPlanePose::linearizeOplus() {
        VertexGroundPlane *v1 = static_cast<VertexGroundPlane *>(_vertices[0]);
        VertexSE3Expmap *v2 = static_cast<VertexSE3Expmap *>(_vertices[0]);

        Vector3d xyz = v2->estimate().inverse().map(Po);
        Vector3d vari_M = v1->estimate();
        Matrix3d xyz_trans = Matrix3d::Zero();
        xyz_trans(0, 1) = -xyz(2, 0);
        xyz_trans(0, 2) = xyz(1, 0);
        xyz_trans(1, 0) = xyz(2, 0);
        xyz_trans(1, 2) = -xyz(0, 0);
        xyz_trans(2, 0) = -xyz(1, 0);
        xyz_trans(2, 1) = xyz(0, 0);


        Matrix<double, 4, 6> POM;
        POM.block(0, 0, 3, 3) = Matrix3d::Identity();
        POM.block(0, 3, 3, 3) = -xyz_trans;
        POM.row(3) = Vector6d::Zero().transpose();

        Vector4d M;
        M << vari_M(0,0),vari_M(1,0),vari_M(2,0),1;

        _jacobianOplusXi = xyz;
        _jacobianOplusXj = - M.transpose() * POM;

    }

}//end namespace