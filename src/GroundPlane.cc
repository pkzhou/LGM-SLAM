/**
* This file is part of ORB-SLAM2.
*
* Copyright (C) 2014-2016 Raúl Mur-Artal <raulmur at unizar dot es> (University of Zaragoza)
* For more information see <https://github.com/raulmur/ORB_SLAM2>
*
* ORB-SLAM2 is free software: you can redistribute it and/or modify
* it under the terms of the GNU General Public License as published by
* the Free Software Foundation, either version 3 of the License, or
* (at your option) any later version.
*
* ORB-SLAM2 is distributed in the hope that it will be useful,
* but WITHOUT ANY WARRANTY; without even the implied warranty of
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
* GNU General Public License for more details.
*
* You should have received a copy of the GNU General Public License
* along with ORB-SLAM2. If not, see <http://www.gnu.org/licenses/>.
*/


#include"KeyFrame.h"
#include"Frame.h"
#include"Map.h"
#include"GroundPlane.h"

#include<opencv2/core/core.hpp>
#include<mutex>

namespace ORB_SLAM2
{
    long unsigned int GroundPlane::nNextId=0;


    GroundPlane::GroundPlane(array<double, 4> _coefs, array<double, 3> _center, array<double,3> maxpoint, KeyFrame *_poriginKF,
                             double range,Eigen::Matrix3d cov):
            mThrange(range), mMgroundcov(cov), coefs(_coefs), center(_center), poriginKF(_poriginKF) , PLanepoint(maxpoint)
    {
        //plane ID
        mnId=nNextId++;
        nMatches = 0;

        //for visualization
        auto MaxP = PLanepoint;
        auto normal = this->GetNormal();
        array<double,3> direction1 = {MaxP[0] - center[0], MaxP[1] - center[1], MaxP[2] - center[2]};
        array<double,3> direction2 = {normal[1] * direction1[2] - normal[2] * direction1[1],
                                      normal[2] * direction1[0] - normal[0] * direction1[2],
                                      normal[0] * direction1[1] - normal[1] * direction1[0]};
        double d1 = cv::norm(direction1, array<double,3>{0, 0, 0});
        double d2 = cv::norm(direction2, array<double,3>{0, 0, 0});
        dir1 = {direction1[0] / d1, direction1[1] / d1, direction1[2] / d1};
        dir2 = {direction2[0] / d2, direction2[1] / d2, direction2[2] / d2};

    }

    void GroundPlane::Setcoefs(Eigen::Vector3d abc)
    {
        //unique_lock<mutex> lock(mMutexcoefs);
        if (coefs[3]!=1)
            cerr<<"can't use vector3d abc because the coefs is not normalized by d."<<endl;
        //assert(coefs[3] == 1);
        coefs[0] = abc(0,0);
        coefs[1] = abc(1,0);
        coefs[2] = abc(2,0);


        //for visualization
        auto MaxP = PLanepoint;
        auto normal = this->GetNormal();
        array<double,3> direction1 = {MaxP[0] - center[0], MaxP[1] - center[1], MaxP[2] - center[2]};
        array<double,3> direction2 = {normal[1] * direction1[2] - normal[2] * direction1[1],
                                      normal[2] * direction1[0] - normal[0] * direction1[2],
                                      normal[0] * direction1[1] - normal[1] * direction1[0]};
        double d1 = cv::norm(direction1, array<double,3>{0, 0, 0});
        double d2 = cv::norm(direction2, array<double,3>{0, 0, 0});
        dir1 = {direction1[0] / d1, direction1[1] / d1, direction1[2] / d1};
        dir2 = {direction2[0] / d2, direction2[1] / d2, direction2[2] / d2};
    }

    array<double,4> GroundPlane::Getcoefs() {

        //unique_lock<mutex> lock(mMutexcoefs);
        array<double,4> worldcoefs;
        worldcoefs = {coefs[0], coefs[1], coefs[2], coefs[3]};
        //Note that the return value is a pointer, and the pointer is also a protection type.
        //所以要copy一份新的指针地址。
        return worldcoefs;
    }

    array<double,3> GroundPlane::GetNormal()
    {
        array<double,3> direction = {coefs[0],coefs[1],coefs[2]};

        double length = cv::norm(direction,array<double, 3>{0,0,0});
        array<double,3> normal = {coefs[0]/length,coefs[1]/length,coefs[2]/length};


        return normal;
    }
    Eigen::Vector3d GroundPlane::GetcoefsEigen(){
        //unique_lock<mutex> lock(mMutexcoefs);
        if (coefs[3]!=1)
            cerr<<"can't use vector3d abc because the coefs is not normalized by d."<<endl;
        return Eigen::Vector3d {coefs[0],coefs[1],coefs[2]};
    }

    array<double,3> GroundPlane::GetPlaneCenter(){
        //unique_lock<mutex> lock(mMutexcenter);
        array<double,3> worldcenter = {center[0],center[1],center[2]};
        return worldcenter;
    }
    void GroundPlane::Getplanedirection(array<double,3> &d1, array<double,3> &d2) const
    {
        d1 = dir1;
        d2 = dir2;
    }
    int GroundPlane::GetpOriginKFid(){
        return poriginKF->mnId;
    }

} //namespace ORB_SLAM

