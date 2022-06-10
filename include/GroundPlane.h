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

#ifndef GROUNDPLANE_H
#define GROUNDPLANE_H

#include"KeyFrame.h"
#include"Frame.h"
#include"Map.h"
#include"MapPoint.h"
#include<Eigen/Core>
#include<opencv2/core/core.hpp>
#include<mutex>

namespace ORB_SLAM2
{

class KeyFrame;
class Map;
class Frame;


class GroundPlane
{
public:
    GroundPlane(array<double,4> _coefs, array<double,3> _center, array<double,3> maxPoint, KeyFrame* _poriginKF,
                double range, Eigen::Matrix3d cov);

    void Setcoefs(Eigen::Vector3d abc);
    array<double,4> Getcoefs();
    Eigen::Vector3d GetcoefsEigen();
    array<double,3> GetPlaneCenter();
    array<double,3> GetNormal();
    void Getplanedirection(array<double,3> &d1, array<double,3> &d2) const;
    int GetpOriginKFid();


public:

    int nMatches;//matched numbers
    long unsigned int mnId;
    double mThrange;
    long unsigned int mnBAlocalforKF = -1;
    //vector<KeyFrame*> mpKFMatches;
    Eigen::Matrix3d mMgroundcov;


    //for visualization
    array<double,3> PLanepoint;
    array<double,3> dir1,dir2;

protected:

    static long unsigned int nNextId;
    array<double,4> coefs;
    array<double,3> center;
    const KeyFrame* poriginKF;

    std::mutex mMutexcoefs;
    std::mutex mMutexcenter;
};

} //namespace ORB_SLAM

#endif // GROUNDPLANE_H
