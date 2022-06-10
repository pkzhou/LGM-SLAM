//
// Created by zpk on 2021/3/5.
//

#ifndef ORB_SLAM2_NEWEDGE_H
#define ORB_SLAM2_NEWEDGE_H
#include "Thirdparty/g2o/g2o/core/base_vertex.h"
#include "Thirdparty/g2o/g2o/core/base_binary_edge.h"
#include "Thirdparty/g2o/g2o/core/base_unary_edge.h"
#include "Thirdparty/g2o/g2o/types/se3_ops.h"
#include "Thirdparty/g2o/g2o/types/se3quat.h"
#include "Thirdparty/g2o/g2o/types/types_sba.h"
#include "Thirdparty/g2o/g2o/types/types_six_dof_expmap.h"
#include "Thirdparty/g2o/g2o/types/types_seven_dof_expmap.h"
#include <Eigen/Geometry>
using namespace Eigen;
using namespace g2o;
namespace ORB_SLAM2
{
    class  VertexGroundPlane : public BaseVertex<3, Eigen::Vector3d>{
    public:
        EIGEN_MAKE_ALIGNED_OPERATOR_NEW

        bool read(std::istream& is){}

        bool write(std::ostream& os) const {}

        virtual void setToOriginImpl() {
            _estimate <<0,0,0;
        }

        virtual void oplusImpl(const double* update_)  {
            _estimate += Vector3d(update_);
        }
    };

    class  EdgeCoplanarPrior: public  BaseBinaryEdge<1, double, VertexGroundPlane, VertexSBAPointXYZ>{
    public:
        EIGEN_MAKE_ALIGNED_OPERATOR_NEW

        EdgeCoplanarPrior(){}

        bool read(std::istream& is) { }

        bool write(std::ostream& os) const {}

        void computeError() override  {
            const VertexGroundPlane* v1 = static_cast<const VertexGroundPlane*>(_vertices[0]);
            const VertexSBAPointXYZ* v2 = static_cast<const VertexSBAPointXYZ*>(_vertices[1]);
            const Vector3d G = v1->estimate();
            const Vector3d P = v2->estimate();

            _error(0,0) =G.transpose()*P+1-_measurement;

        }

        virtual void linearizeOplus();
    };



    class  EdgeGroundProjectOnlyPose: public  BaseUnaryEdge<1, double, VertexSE3Expmap>{
    public:
        EIGEN_MAKE_ALIGNED_OPERATOR_NEW

        EdgeGroundProjectOnlyPose(){}

        bool read(std::istream& is);

        bool write(std::ostream& os) const;

        void computeError() override  {
            const VertexSE3Expmap* v1 = static_cast<const VertexSE3Expmap*>(_vertices[0]);
            //yi Matrix<double,1,1> seems perform trouble when running instead of compiling.
            Vector3d Pw = v1->estimate().inverse().map(Po);
            Vector3d M;
            M<< a,b,c;
            double axbycz = (M.transpose()*Pw)(0,0);
            _error(0,0) = _measurement-(axbycz+d);
        }

        /*bool isDepthPositive() {
            const VertexSE3Expmap* v1 = static_cast<const VertexSE3Expmap*>(_vertices[0]);
            return (v1->estimate().map(Xw))(2)>0.0;
        }*/

        virtual void linearizeOplus();

        Vector3d Po;
        double a,b,c,d;
        SE3Quat Tiw;//not need
    };

    class  EdgeGroundProjectPlanePose: public  BaseBinaryEdge<1, double, VertexGroundPlane, VertexSE3Expmap>{
    public:
        EIGEN_MAKE_ALIGNED_OPERATOR_NEW

        EdgeGroundProjectPlanePose(){}

        bool read(std::istream& is){}

        bool write(std::ostream& os) const { }

        void computeError() override  {
            const VertexGroundPlane* v1 = static_cast<const VertexGroundPlane*>(_vertices[0]);
            const VertexSE3Expmap* v2 = static_cast<const VertexSE3Expmap*>(_vertices[1]);
            //yi Matrix<double,1,1> seems perform trouble when running instead of compiling.
            //It is solved . The reason is dimension or something.
            Vector3d Pw = v2->estimate().inverse().map(Po);
            Vector3d M;
            M = v1->estimate();
            double axbycz = (M.transpose()*Pw)(0,0);
            _error(0,0) = (axbycz+1)-_measurement;
        }

        /*bool isDepthPositive() {
            const VertexSE3Expmap* v1 = static_cast<const VertexSE3Expmap*>(_vertices[0]);
            return (v1->estimate().map(Xw))(2)>0.0;
        }*/

        virtual void linearizeOplus();

        Vector3d Po;
    };

    class  EdgeGroundProjectRelPose: public  BaseBinaryEdge<1, double, VertexSE3Expmap, VertexSE3Expmap>{
    public:
        EIGEN_MAKE_ALIGNED_OPERATOR_NEW

        EdgeGroundProjectRelPose(){}

        bool read(std::istream& is);

        bool write(std::ostream& os) const;

        void computeError() override  {
            const VertexSE3Expmap* vi = static_cast<const VertexSE3Expmap*>(_vertices[0]);
            const VertexSE3Expmap* vj = static_cast<const VertexSE3Expmap*>(_vertices[1]);
            //yi Matrix<double,1,1> seems perform trouble when running instead of compiling.
            //It is solved . The reason is dimension or something.
            Vector3d Pw = vj->estimate().inverse().map(Po);
            Vector3d M;
            M<< a,b,c;
            SE3Quat Tiw = vi->estimate();
            double axbycz = (M.transpose()*(Tiw.map(Pw)))(0,0);
            _error(0,0) = _measurement-(axbycz+d);
        }

        /*bool isDepthPositive() {
            const VertexSE3Expmap* v1 = static_cast<const VertexSE3Expmap*>(_vertices[0]);
            return (v1->estimate().map(Xw))(2)>0.0;
        }*/

        virtual void linearizeOplus();

        Vector3d Po;
        double a,b,c,d;
    };
}//end namespace
#endif //ORB_SLAM2_NEWEDGE_H
