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

#include "LocalMapping.h"
#include "LoopClosing.h"
#include "ORBmatcher.h"
#include "Optimizer.h"
#include <unistd.h>
#include <mutex>

#include "GRANSAC.hpp"
#include "PlaneModel.hpp"
#include <iostream>
namespace ORB_SLAM2
{

LocalMapping::LocalMapping(Map *pMap, const float bMonocular):
    mbMonocular(bMonocular), mbResetRequested(false), mbFinishRequested(false), mbFinished(true), mpMap(pMap),
    mbAbortBA(false), mbStopped(false), mbStopRequested(false), mbNotStop(false), mbAcceptKeyFrames(true)
{
}

void LocalMapping::SetLoopCloser(LoopClosing* pLoopCloser)
{
    mpLoopCloser = pLoopCloser;
}

void LocalMapping::SetTracker(Tracking *pTracker)
{
    mpTracker=pTracker;
}

vector<array<double, 3>> LocalMapping::GetgrounddepthfromStereo()
{
    vector<array<double, 3>> groundp;
    vector<int> id;

    const float Thdepth = mpCurrentKeyFrame->mThDepth/2;
    const float fx = mpCurrentKeyFrame->fx;
    const float fy = mpCurrentKeyFrame->fy;
    const float cx = mpCurrentKeyFrame->cx;
    const float cy = mpCurrentKeyFrame->cy;

    const auto &kp = mpCurrentKeyFrame->mvKeysUn;
    const auto &kpd = mpCurrentKeyFrame->mvDepth;
    const cv::Mat Tjw = mpCurrentKeyFrame->GetPose();
    cv::Mat Rjw = Tjw.rowRange(0, 3).colRange(0, 3);
    cv::Mat tjw = Tjw.rowRange(0, 3).col(3);
    cv::Mat Rwj = Tjw.rowRange(0, 3).colRange(0, 3).t();
    cv::Mat twj = -Rwj * Tjw.rowRange(0, 3).col(3);


    int N = kp.size();
    int groundn = 0,map1 = 0,map2 = 0,map3 = 0;
    for (int i = 0; i < N; ++i)
    {
        const auto p = kp[i];
        if (p.class_id==1)
        {

            MapPoint *pmap = mpCurrentKeyFrame->GetMapPoint(i);
            if (pmap)
            {
                map1++;
                cv::Mat wpos = pmap->GetWorldPos();
                cv::Mat localpos = Rjw*wpos+tjw;
                array<double, 3> point = {wpos.at<float>(0,0),wpos.at<float>(1,0),wpos.at<float>(2,0)};
                if(localpos.at<float>(2,0)>0&&localpos.at<float>(2,0)<Thdepth)
                    {groundp.push_back(point);
                        map2++;
                        groundn++;
                    }
//               {
//                    const float d = kpd[i];
//                    array<double, 3> point2 = {((p.pt.x - cx) / fx * d) , (p.pt.y - cy) / fy * d, d};
//                    cout<<"Mappoint :"<<point[0]<<" "<<point[1]<<" "<<point[2]<<" \n";
//                    cout<<"Keypoint :"<<point2[0]<<" "<<point2[1]<<" "<<point2[2]<<" \n";
//               }
            }
            else
            {
                const float d = kpd[i];
                map3++;
                if(d>0&&d<Thdepth)
                {
                    cv::Mat localpos = (cv::Mat_<float>(3, 1) << (p.pt.x - cx)/fx*d, (p.pt.y - cy)/fy*d, d );
                    cv::Mat wpos = Rwj*localpos+twj;
                        array<double, 3> point = {wpos.at<float>(0,0),wpos.at<float>(1,0),wpos.at<float>(2,0)};
                        groundp.push_back(point);
                        groundn++;

                }  
            }
            
        }
    }
    //cout<<"The size of all the points get: "<<groundn<<endl;
    return groundp;

}
void LocalMapping::GMapPointsMatch()
{
    const auto &kp = mpCurrentKeyFrame->mvKeysUn;
    array<double,3> center = mpCurrentKeyFrame->mpfittingPlane->GetPlaneCenter();
    cv::Mat ccenter = (cv::Mat_<float>(3, 1) << center[0],center[1],center[2]);
    int N = kp.size();
    for (int i = 0; i < N; ++i) {
        const auto p = kp[i];
        if (p.class_id == 1) {
            MapPoint *pmap = mpCurrentKeyFrame->GetMapPoint(i);
            if (pmap)
            {
                cv::Mat wpos = pmap->GetWorldPos();
                double distance = cv::norm(wpos,ccenter);
                if (distance<mlastth)
                    if (distance<(pmap->ground_dist))
                    {
                        pmap->ground_dist = distance;
                        pmap->mpMatchKF = mpCurrentKeyFrame;
                        pmap->num++;
                        //cout<<"MapPoint : "<<pmap->mnId<<" change its MatchKeyFrame "<<pmap->num<<"th"<<endl;
                    }
            }
        }
    }

}

void LocalMapping::GroundFitting()
{
    int new_N;//new points observation
    int total_N;//all of the points used
    int inliners_N;//RANSAC inliners number


    std::vector<std::array<double, 3>> groundp;
    std::array<double,3> center;
    std::array<double,4> coefs;

    double inlinerfraction = 0;

    if(mpCurrentKeyFrame->mvgroundcandis.empty())
    {
        groundp = GetgrounddepthfromStereo();
    }
    else
    {groundp = mpCurrentKeyFrame->mvgroundcandis;}
    new_N = groundp.size();
    if(new_N < 15)
    {
        mpCurrentKeyFrame->mbgroundflag = false;
        return;
    }

    for (int i = 0; i < 3; i++)
    {
        double temp = 0;
        for (int j = 0; j < new_N; ++j) {
            temp = temp+groundp[j][i];
        }
        center[i] = temp/new_N;
    }
    //cout<<"center is : "<<center[0]<<" "<<center[1]<<" "<<center[2]<<endl;

    //计算该localground可容忍的最小距离（mThreshold）
    double thre = 0.0;
    for (auto point : groundp) {
        double distance = cv::norm(center,point);
        if (distance>thre) {
            thre = distance;
        }
    }
//    cv::Mat currentTcw = mpCurrentKeyFrame->GetPose();
//    cv::Mat Rwc = currentTcw.rowRange(0, 3).colRange(0, 3).t();
//    cv::Mat twc = -Rwc * currentTcw.rowRange(0, 3).col(3);
//    cv::Mat Rcw = Rwc.t();
//    cv::Mat tcw = currentTcw.rowRange(0,3).col(3);
    int flag = 1;
    for (int i = 0; i < new_N; ++i)
    {
        if (i<mlastgpoints.size())
        {
            array<double,3> lpoint = mlastgpoints[i];
                for (int j = 0; j < new_N; ++j) {
                    double distance = cv::norm(lpoint, groundp[j]);
                    //cout << " jth distance is " << j << ", " << distance << endl;
                    if (distance < 0.4 || distance > 1.6*thre) {
                        flag = 0;
                        break;
                    }
                    else
                    {
                        if (distance<1)
                            break;
                    }
                }
                if (flag) {
                    groundp.push_back(lpoint);
                }

                flag = 1;
        }
    }
    total_N = groundp.size();
    //cout<<"lastpoint size is : "<<mlastgpoints.size()<<endl;
    //cout << "new groundpoint size is :" << new_N << " and all size of the groundpoint is : " << total_N << endl;


	std::vector<std::shared_ptr<GRANSAC::AbstractParameter>> CandPoints;

//#pragma omp parallel for num_threads(6)
	for (int i = 0; i < total_N; ++i) {
        Vector3VP p = groundp[i];
            std::shared_ptr<GRANSAC::AbstractParameter> CandPt = std::make_shared<Point3D>(p[0], p[1],p[2]);
            CandPoints.push_back(CandPt);

	}
	GRANSAC::RANSAC<PlaneModel, 3> Estimator;
    
    Estimator.Initialize(0.005, 100); // Threshold, iterations
    //int64_t start = cv::getTickCount();
	Estimator.Estimate(CandPoints, inlinerfraction);

	//Least Square solve
    const std::vector<std::shared_ptr<GRANSAC::AbstractParameter>>& inliners = Estimator.GetBestInliers();
    inliners_N = inliners.size();

    std::vector<std::array<double, 3>> currentgpoint;

    Eigen::MatrixXd A(inliners_N, 3);
    Eigen::MatrixXd b(inliners_N, 1);
    for (int i = 0; i < inliners_N; ++i) {
        auto ExtPoint3D = std::dynamic_pointer_cast<Point3D>(inliners[i]);
        currentgpoint.push_back({ExtPoint3D->m_Point3D[0],ExtPoint3D->m_Point3D[1],ExtPoint3D->m_Point3D[2]});

        Eigen::Vector3d t(ExtPoint3D->m_Point3D[0],ExtPoint3D->m_Point3D[1],ExtPoint3D->m_Point3D[2]);
        A.block(i,0,1,3) = t.transpose();
        b(i,0) = - 1;

    }
    Eigen::JacobiSVD<Eigen ::MatrixXd> svd(A, Eigen ::ComputeThinU | Eigen ::ComputeThinV);
    Eigen::MatrixXd temp = svd.solve(b);
    Eigen::Vector4d  newcoefs {temp(0,0),temp(1,0),temp(2,0),1};

    //int64_t end = cv::getTickCount();
    //std::cout << "RANSAC took: " << GRANSAC::VPFloat(end - start) / GRANSAC::VPFloat(cv::getTickFrequency()) * 1000.0 << " ms." << std::endl;
    //cout<<"inlinerfraction is :"<<inlinerfraction<<endl;
    auto BestPlane = Estimator.GetBestModel();

	if (BestPlane == nullptr)
	{
        mpCurrentKeyFrame->mbgroundflag = false;
	    return;
	}

    for (int i = 0; i < 4; i++) {
        coefs[i] = newcoefs(i,0);
    }

    //     //cout<<"-----------------vertical computing-------------------"<<endl;
    // if (mpCurrentKeyFrame->mpReferenceKF)
    // {
    //     cv::Mat ti = mpCurrentKeyFrame->GetTranslation();
    //     cv::Mat tj = mpCurrentKeyFrame->mpReferenceKF->GetTranslation();
    //     cv::Mat Rwc = mpCurrentKeyFrame->GetRotation().t();
    //     cv::Mat t = ti-tj;
    //     double tnor = sqrt(pow(t.at<float>(0,0),2)+pow(t.at<float>(1,0),2)+pow(t.at<float>(2,0),2));
    //     t = t/tnor;
    //     cv::Mat r1 = Rwc.col(0);
    //     double r1nor = sqrt(pow(r1.at<float>(0,0),2)+pow(r1.at<float>(1,0),2)+pow(r1.at<float>(2,0),2));
    //     r1 = r1/r1nor;

    //     double cunor = sqrt(coefs[0] * coefs[0] + coefs[1] * coefs[1] + coefs[2] * coefs[2]);
    //     cv::Mat curr_normal =(cv::Mat_<float>(3, 1) << coefs[0] / cunor, coefs[1] / cunor, coefs[2] / cunor );
    //     cv::Mat verticality1 = t.t()*curr_normal;
    //     cv::Mat verticality2 = r1.t()*curr_normal;
    //     //cos similarity compute, no morethan +-5angle
    //     //cout<<"verticality is : "<<verticality<<endl;
    //     if (abs(verticality1.at<float>(0,0))>0.173||abs(verticality2.at<float>(0,0))>0.1)
    //     //if (abs(verticality1.at<float>(0,0))>0.173)
    //     {
    //         mpCurrentKeyFrame->mbgroundflag = false;
    //         return;
    //     }
    // }
    std::vector<std::array<double, 4>> mlast4normal;
    vector<KeyFrame*> vpKFs = mpMap->GetAllKeyFrames();
    sort(vpKFs.begin(), vpKFs.end(), KeyFrame::lId);
    if(vpKFs.size()>4)
    {
        int i=0;
    vector<KeyFrame *>::iterator pKF = vpKFs.end() - 2;
    while (i<4)
   {

       if ((*pKF)->mbgroundflag) {
                auto plane = (*pKF)->mpfittingPlane;
                if (plane)
                {
                mlast4normal.push_back(plane->Getcoefs());
                i++;
                }
       }
       if (pKF> vpKFs.begin())
       {
          pKF--;
       }
       else 
       break;
   }}
   
   
    //cout<<"-----------------parallelism computing-------------------"<<endl;
      cv::Mat Rwc = mpCurrentKeyFrame->GetRotation().t();
      cv::Mat v_n=-Rwc.col(1);
      double tnor = sqrt(pow(v_n.at<float>(0,0),2)+pow(v_n.at<float>(1,0),2)+pow(v_n.at<float>(2,0),2));
      v_n=v_n/tnor;
    double cunor = sqrt(coefs[0] * coefs[0] + coefs[1] * coefs[1] + coefs[2] * coefs[2]);
    cv::Mat curr_normal =(cv::Mat_<float>(3, 1) << coefs[0] / cunor, coefs[1] / cunor, coefs[2] / cunor );
    cv::Mat Parallelism = v_n.t()*curr_normal;
    if (abs(Parallelism.at<float>(0,0))<0.866)
         {
             mpCurrentKeyFrame->mbgroundflag = false;
             return;
         }
         else
         {
             if (abs(Parallelism.at<float>(0,0))<0.98)
             {
                 if(mlast4normal.size()==4)
                 {std::array<double,4> temp;
                 for (size_t i = 0; i < 4; i++)
                 {
                     //cout<<"加权平均start"<<endl;
                     temp[i]=mlast4normal[0][i]+mlast4normal[1][i]+mlast4normal[2][i]+mlast4normal[3][i]+2*coefs[i];
                    temp[i]=temp[i]/6;
                    //cout<<"加权平均end"<<endl;
                 }
                 coefs=temp;}
             }
             
         }
         //cout<<"succeed."<<endl;

    Eigen::MatrixXd AtA = A.transpose()*A;
    Eigen::MatrixXd Sigma = AtA.inverse();
    //cout<<"worldcoefs is : "<<coefs[0]<<" "<<coefs[1]<<" "<<coefs[2]<<" "<<coefs[3]<<endl;

    //存入Keyframe
    mpCurrentKeyFrame->mpfittingPlane = new GroundPlane(coefs, center, currentgpoint[0], mpCurrentKeyFrame, thre, Sigma);
    mpCurrentKeyFrame->mbgroundflag = true;
    //update the last groundnts
    mlastgpoints.clear();
    mlastgpoints = groundp;
    mlastth = thre;
    //mlgroundparameter = coefs;

    GMapPointsMatch();

    cv::Mat wcenter = (cv::Mat_<float>(3, 1) << center[0], center[1], center[2]);
    cv::Mat wcoefs = (cv::Mat_<float>(1, 4) << coefs[0], coefs[1], coefs[2], coefs[3]);
    cv::Mat wpoints = cv::Mat_<float>(inliners_N, 3);//inliner points

    //wcenter = Rwc * wcenter + twc;
    //wcoefs = wcoefs * currentTcw;



    //无回环情况下，Tcw0应该是单位的，直接从keyframe中获取
	if(GROUND_RECORD)
    {
        groundout << wcenter.at<float>(0, 0) << " " << wcenter.at<float>(1, 0) << " " << wcenter.at<float>(2, 0)
                  << " ";
        groundout << wcoefs.at<float>(0, 0) << " " << wcoefs.at<float>(0, 1) << " " << wcoefs.at<float>(0, 2) << " "
                  << wcoefs.at<float>(0, 3);
        groundout << endl;

        for(int i=0; i < inliners_N; i++)
        {
            auto point3d = currentgpoint[i];
            pointout<< point3d[0] <<" "<<point3d[1]<<" "<<point3d[2]<<" ";
        }
        pointout<<endl;
    }
}

void LocalMapping::Run()
{

    mbFinished = false;
    if(GROUND_RECORD)
    {
        groundout.open( "ground.txt" ,ios_base::out);
        pointout.open( "points.txt" ,ios_base::out);
    }


    while(1)
    {
        // Tracking will see that Local Mapping is busy
        SetAcceptKeyFrames(false);
        
        // Check if there are keyframes in the queue
        if(CheckNewKeyFrames())
        {
            // BoW conversion and insertion in Map
            //首先将新的关键帧Ki作为新的节点Ki加入Covibility Graph，
            //并且更新与那些能够共享地图点的关键帧节点相连接的边。
            ProcessNewKeyFrame();
            


            // Check recent MapPoints
            //为了保存地图点，必须在创建该点云的前三帧测试通过约束，
            //才能真正被保存，这样才能保证可跟踪且不容易在三角化时出现较大误差
            /*
            一个点要被加入Map，需要满足下面条件： 
（1）这个点要在可预测到能够观察到该点的关键帧中，有超过25%的关键帧能够跟踪到这个点； 
（2）如果一个地图点被构建，它必须被超过三个关键帧观察到（在代码中，可以发现如果是单摄像头，这个阈值被设置为2）。
            */
            MapPointCulling();

            // Triangulate new MapPoints
            CreateNewMapPoints();



            if(!CheckNewKeyFrames())
            {
                // Find more matches in neighbor keyframes and fuse point duplications
                SearchInNeighbors();
            }

            mbAbortBA = false;

            if(!CheckNewKeyFrames() && !stopRequested())
            {
                // Local BA
                if(mpMap->KeyFramesInMap()>2)
                    Optimizer::LocalBundleAdjustment(mpCurrentKeyFrame,&mbAbortBA, mpMap);
                //计算该关键帧对应的地面方程和中心位置，并返回inliner
                //地面拟合,结果在inliner和参数向量
                GroundFitting();
                // Check redundant local Keyframes
                KeyFrameCulling();
            }

            mpLoopCloser->InsertKeyFrame(mpCurrentKeyFrame);
        }
        else if(Stop())
        {
            // Safe area to stop
            while(isStopped() && !CheckFinish())
            {
                usleep(3000);
            }
            if(CheckFinish())
                break;
        }

        ResetIfRequested();

        // Tracking will see that Local Mapping is busy
        //完成LocalMapping工作，将标志设为SetAcceptKeyFrames(true)，
        //以允许Tracking线程继续得到关键帧
        SetAcceptKeyFrames(true);

        if(CheckFinish())
            break;

        usleep(3000);
    }

    SetFinish();
}

void LocalMapping::InsertKeyFrame(KeyFrame *pKF)
{
    unique_lock<mutex> lock(mMutexNewKFs);
    mlNewKeyFrames.push_back(pKF);
    mbAbortBA=true;
}


bool LocalMapping::CheckNewKeyFrames()
{
    unique_lock<mutex> lock(mMutexNewKFs);
    return(!mlNewKeyFrames.empty());
}

void LocalMapping::ProcessNewKeyFrame()
{
    {
        unique_lock<mutex> lock(mMutexNewKFs);
        mpCurrentKeyFrame = mlNewKeyFrames.front();
        mlNewKeyFrames.pop_front();
    }

    // Compute Bags of Words structures
    mpCurrentKeyFrame->ComputeBoW();

    // Associate MapPoints to the new keyframe and update normal and descriptor
    const vector<MapPoint*> vpMapPointMatches = mpCurrentKeyFrame->GetMapPointMatches();

    for(size_t i=0; i<vpMapPointMatches.size(); i++)
    {
        MapPoint* pMP = vpMapPointMatches[i];
        if(pMP)
        {
            if(!pMP->isBad())
            {   //如果不在关键帧内
                if(!pMP->IsInKeyFrame(mpCurrentKeyFrame))
                {
                    pMP->AddObservation(mpCurrentKeyFrame, i);
                    pMP->UpdateNormalAndDepth();
                    pMP->ComputeDistinctiveDescriptors();
                }
                else // this can only happen for new stereo points inserted by the Tracking
                {
                    mlpRecentAddedMapPoints.push_back(pMP);
                }
            }
        }
    }    

    // Update links in the Covisibility Graph
    mpCurrentKeyFrame->UpdateConnections();



    // Insert Keyframe in Map
    mpMap->AddKeyFrame(mpCurrentKeyFrame);
}

void LocalMapping::MapPointCulling()
{
    // Check Recent Added MapPoints
    list<MapPoint*>::iterator lit = mlpRecentAddedMapPoints.begin();
    const unsigned long int nCurrentKFid = mpCurrentKeyFrame->mnId;

    int nThObs;
    if(mbMonocular)
        nThObs = 2;
    else
        nThObs = 3;
    const int cnThObs = nThObs;

    while(lit!=mlpRecentAddedMapPoints.end())
    {
        MapPoint* pMP = *lit;
        if(pMP->isBad())
        {
            lit = mlpRecentAddedMapPoints.erase(lit);
        }
        else if(pMP->GetFoundRatio()<0.25f )
        {
            pMP->SetBadFlag();
            lit = mlpRecentAddedMapPoints.erase(lit);
        }
        else if(((int)nCurrentKFid-(int)pMP->mnFirstKFid)>=2 && pMP->Observations()<=cnThObs)
        {
            pMP->SetBadFlag();
            lit = mlpRecentAddedMapPoints.erase(lit);
        }
        else if(((int)nCurrentKFid-(int)pMP->mnFirstKFid)>=3)
            lit = mlpRecentAddedMapPoints.erase(lit);
        else
            lit++;
    }
}
//单目视觉地图点的生成主要在这里，对关键帧三角化地图点
//所以跟踪时，一般都是跟踪参考关键帧而非上一帧，因为上一帧不一定是关键帧，可能没有地图点
void LocalMapping::CreateNewMapPoints()
{
    // Retrieve neighbor keyframes in covisibility graph
    int nn = 10;
    if(mbMonocular)
        nn=20;
    const vector<KeyFrame*> vpNeighKFs = mpCurrentKeyFrame->GetBestCovisibilityKeyFrames(nn);

    ORBmatcher matcher(0.6,false);

    cv::Mat Rcw1 = mpCurrentKeyFrame->GetRotation();
    cv::Mat Rwc1 = Rcw1.t();
    cv::Mat tcw1 = mpCurrentKeyFrame->GetTranslation();
    cv::Mat Tcw1(3,4,CV_32F);
    Rcw1.copyTo(Tcw1.colRange(0,3));
    tcw1.copyTo(Tcw1.col(3));
    cv::Mat Ow1 = mpCurrentKeyFrame->GetCameraCenter();

    const float &fx1 = mpCurrentKeyFrame->fx;
    const float &fy1 = mpCurrentKeyFrame->fy;
    const float &cx1 = mpCurrentKeyFrame->cx;
    const float &cy1 = mpCurrentKeyFrame->cy;
    const float &invfx1 = mpCurrentKeyFrame->invfx;
    const float &invfy1 = mpCurrentKeyFrame->invfy;

    const float ratioFactor = 1.5f*mpCurrentKeyFrame->mfScaleFactor;

    int nnew=0;

    // Search matches with epipolar restriction and triangulate
    for(size_t i=0; i<vpNeighKFs.size(); i++)
    {
        if(i>0 && CheckNewKeyFrames())
            return;

        KeyFrame* pKF2 = vpNeighKFs[i];

        // Check first that baseline is not too short
        cv::Mat Ow2 = pKF2->GetCameraCenter();
        cv::Mat vBaseline = Ow2-Ow1;
        const float baseline = cv::norm(vBaseline);

        if(!mbMonocular)
        {
            if(baseline<pKF2->mb)
            continue;
        }
        else
        {
            const float medianDepthKF2 = pKF2->ComputeSceneMedianDepth(2);
            const float ratioBaselineDepth = baseline/medianDepthKF2;

            if(ratioBaselineDepth<0.01)
                continue;
        }

        // Compute Fundamental Matrix
        cv::Mat F12 = ComputeF12(mpCurrentKeyFrame,pKF2);

        // Search matches that fullfil epipolar constraint
        vector<pair<size_t,size_t> > vMatchedIndices;
        matcher.SearchForTriangulation(mpCurrentKeyFrame,pKF2,F12,vMatchedIndices,false);

        cv::Mat Rcw2 = pKF2->GetRotation();
        cv::Mat Rwc2 = Rcw2.t();
        cv::Mat tcw2 = pKF2->GetTranslation();
        cv::Mat Tcw2(3,4,CV_32F);
        Rcw2.copyTo(Tcw2.colRange(0,3));
        tcw2.copyTo(Tcw2.col(3));

        const float &fx2 = pKF2->fx;
        const float &fy2 = pKF2->fy;
        const float &cx2 = pKF2->cx;
        const float &cy2 = pKF2->cy;
        const float &invfx2 = pKF2->invfx;
        const float &invfy2 = pKF2->invfy;

        // Triangulate each match
        const int nmatches = vMatchedIndices.size();
        for(int ikp=0; ikp<nmatches; ikp++)
        {
            const int &idx1 = vMatchedIndices[ikp].first;
            const int &idx2 = vMatchedIndices[ikp].second;

            const cv::KeyPoint &kp1 = mpCurrentKeyFrame->mvKeysUn[idx1];
            const float kp1_ur=mpCurrentKeyFrame->mvuRight[idx1];
            bool bStereo1 = kp1_ur>=0;

            const cv::KeyPoint &kp2 = pKF2->mvKeysUn[idx2];
            const float kp2_ur = pKF2->mvuRight[idx2];
            bool bStereo2 = kp2_ur>=0;

            // Check parallax between rays
            cv::Mat xn1 = (cv::Mat_<float>(3,1) << (kp1.pt.x-cx1)*invfx1, (kp1.pt.y-cy1)*invfy1, 1.0);
            cv::Mat xn2 = (cv::Mat_<float>(3,1) << (kp2.pt.x-cx2)*invfx2, (kp2.pt.y-cy2)*invfy2, 1.0);

            cv::Mat ray1 = Rwc1*xn1;
            cv::Mat ray2 = Rwc2*xn2;
            const float cosParallaxRays = ray1.dot(ray2)/(cv::norm(ray1)*cv::norm(ray2));

            float cosParallaxStereo = cosParallaxRays+1;
            float cosParallaxStereo1 = cosParallaxStereo;
            float cosParallaxStereo2 = cosParallaxStereo;

            if(bStereo1)
                cosParallaxStereo1 = cos(2*atan2(mpCurrentKeyFrame->mb/2,mpCurrentKeyFrame->mvDepth[idx1]));
            else if(bStereo2)
                cosParallaxStereo2 = cos(2*atan2(pKF2->mb/2,pKF2->mvDepth[idx2]));

            cosParallaxStereo = min(cosParallaxStereo1,cosParallaxStereo2);

            cv::Mat x3D;
            if(cosParallaxRays<cosParallaxStereo && cosParallaxRays>0 && (bStereo1 || bStereo2 || cosParallaxRays<0.9998))
            {
                // Linear Triangulation Method
                cv::Mat A(4,4,CV_32F);
                A.row(0) = xn1.at<float>(0)*Tcw1.row(2)-Tcw1.row(0);
                A.row(1) = xn1.at<float>(1)*Tcw1.row(2)-Tcw1.row(1);
                A.row(2) = xn2.at<float>(0)*Tcw2.row(2)-Tcw2.row(0);
                A.row(3) = xn2.at<float>(1)*Tcw2.row(2)-Tcw2.row(1);

                cv::Mat w,u,vt;
                cv::SVD::compute(A,w,u,vt,cv::SVD::MODIFY_A| cv::SVD::FULL_UV);

                x3D = vt.row(3).t();

                if(x3D.at<float>(3)==0)
                    continue;

                // Euclidean coordinates
                x3D = x3D.rowRange(0,3)/x3D.at<float>(3);

            }
            else if(bStereo1 && cosParallaxStereo1<cosParallaxStereo2)
            {
                x3D = mpCurrentKeyFrame->UnprojectStereo(idx1);                
            }
            else if(bStereo2 && cosParallaxStereo2<cosParallaxStereo1)
            {
                x3D = pKF2->UnprojectStereo(idx2);
            }
            else
                continue; //No stereo and very low parallax

            cv::Mat x3Dt = x3D.t();

            //Check triangulation in front of cameras
            float z1 = Rcw1.row(2).dot(x3Dt)+tcw1.at<float>(2);
            if(z1<=0)
                continue;

            float z2 = Rcw2.row(2).dot(x3Dt)+tcw2.at<float>(2);
            if(z2<=0)
                continue;

            //Check reprojection error in first keyframe
            const float &sigmaSquare1 = mpCurrentKeyFrame->mvLevelSigma2[kp1.octave];
            const float x1 = Rcw1.row(0).dot(x3Dt)+tcw1.at<float>(0);
            const float y1 = Rcw1.row(1).dot(x3Dt)+tcw1.at<float>(1);
            const float invz1 = 1.0/z1;

            if(!bStereo1)
            {
                float u1 = fx1*x1*invz1+cx1;
                float v1 = fy1*y1*invz1+cy1;
                float errX1 = u1 - kp1.pt.x;
                float errY1 = v1 - kp1.pt.y;
                if((errX1*errX1+errY1*errY1)>5.991*sigmaSquare1)
                    continue;
            }
            else
            {
                float u1 = fx1*x1*invz1+cx1;
                float u1_r = u1 - mpCurrentKeyFrame->mbf*invz1;
                float v1 = fy1*y1*invz1+cy1;
                float errX1 = u1 - kp1.pt.x;
                float errY1 = v1 - kp1.pt.y;
                float errX1_r = u1_r - kp1_ur;
                if((errX1*errX1+errY1*errY1+errX1_r*errX1_r)>7.8*sigmaSquare1)
                    continue;
            }

            //Check reprojection error in second keyframe
            const float sigmaSquare2 = pKF2->mvLevelSigma2[kp2.octave];
            const float x2 = Rcw2.row(0).dot(x3Dt)+tcw2.at<float>(0);
            const float y2 = Rcw2.row(1).dot(x3Dt)+tcw2.at<float>(1);
            const float invz2 = 1.0/z2;
            if(!bStereo2)
            {
                float u2 = fx2*x2*invz2+cx2;
                float v2 = fy2*y2*invz2+cy2;
                float errX2 = u2 - kp2.pt.x;
                float errY2 = v2 - kp2.pt.y;
                if((errX2*errX2+errY2*errY2)>5.991*sigmaSquare2)
                    continue;
            }
            else
            {
                float u2 = fx2*x2*invz2+cx2;
                float u2_r = u2 - mpCurrentKeyFrame->mbf*invz2;
                float v2 = fy2*y2*invz2+cy2;
                float errX2 = u2 - kp2.pt.x;
                float errY2 = v2 - kp2.pt.y;
                float errX2_r = u2_r - kp2_ur;
                if((errX2*errX2+errY2*errY2+errX2_r*errX2_r)>7.8*sigmaSquare2)
                    continue;
            }

            //Check scale consistency
            cv::Mat normal1 = x3D-Ow1;
            float dist1 = cv::norm(normal1);

            cv::Mat normal2 = x3D-Ow2;
            float dist2 = cv::norm(normal2);

            if(dist1==0 || dist2==0)
                continue;

            const float ratioDist = dist2/dist1;
            const float ratioOctave = mpCurrentKeyFrame->mvScaleFactors[kp1.octave]/pKF2->mvScaleFactors[kp2.octave];

            /*if(fabs(ratioDist-ratioOctave)>ratioFactor)
                continue;*/
            if(ratioDist*ratioFactor<ratioOctave || ratioDist>ratioOctave*ratioFactor)
                continue;

            // Triangulation is succesfull
            MapPoint* pMP = new MapPoint(x3D,mpCurrentKeyFrame,mpMap);

            pMP->AddObservation(mpCurrentKeyFrame,idx1);            
            pMP->AddObservation(pKF2,idx2);

            mpCurrentKeyFrame->AddMapPoint(pMP,idx1);
            pKF2->AddMapPoint(pMP,idx2);

            pMP->ComputeDistinctiveDescriptors();

            pMP->UpdateNormalAndDepth();

            mpMap->AddMapPoint(pMP);
            mlpRecentAddedMapPoints.push_back(pMP);

            nnew++;
        }
    }
}

void LocalMapping::SearchInNeighbors()
{
    // Retrieve neighbor keyframes
    int nn = 10;
    if(mbMonocular)
        nn=20;
    const vector<KeyFrame*> vpNeighKFs = mpCurrentKeyFrame->GetBestCovisibilityKeyFrames(nn);
    vector<KeyFrame*> vpTargetKFs;
    for(vector<KeyFrame*>::const_iterator vit=vpNeighKFs.begin(), vend=vpNeighKFs.end(); vit!=vend; vit++)
    {
        KeyFrame* pKFi = *vit;
        if(pKFi->isBad() || pKFi->mnFuseTargetForKF == mpCurrentKeyFrame->mnId)
            continue;
        vpTargetKFs.push_back(pKFi);
        pKFi->mnFuseTargetForKF = mpCurrentKeyFrame->mnId;

        // Extend to some second neighbors
        const vector<KeyFrame*> vpSecondNeighKFs = pKFi->GetBestCovisibilityKeyFrames(5);
        for(vector<KeyFrame*>::const_iterator vit2=vpSecondNeighKFs.begin(), vend2=vpSecondNeighKFs.end(); vit2!=vend2; vit2++)
        {
            KeyFrame* pKFi2 = *vit2;
            if(pKFi2->isBad() || pKFi2->mnFuseTargetForKF==mpCurrentKeyFrame->mnId || pKFi2->mnId==mpCurrentKeyFrame->mnId)
                continue;
            vpTargetKFs.push_back(pKFi2);
        }
    }


    // Search matches by projection from current KF in target KFs
    ORBmatcher matcher;
    vector<MapPoint*> vpMapPointMatches = mpCurrentKeyFrame->GetMapPointMatches();
    for(vector<KeyFrame*>::iterator vit=vpTargetKFs.begin(), vend=vpTargetKFs.end(); vit!=vend; vit++)
    {
        KeyFrame* pKFi = *vit;

        matcher.Fuse(pKFi,vpMapPointMatches);
    }

    // Search matches by projection from target KFs in current KF
    vector<MapPoint*> vpFuseCandidates;
    vpFuseCandidates.reserve(vpTargetKFs.size()*vpMapPointMatches.size());

    for(vector<KeyFrame*>::iterator vitKF=vpTargetKFs.begin(), vendKF=vpTargetKFs.end(); vitKF!=vendKF; vitKF++)
    {
        KeyFrame* pKFi = *vitKF;

        vector<MapPoint*> vpMapPointsKFi = pKFi->GetMapPointMatches();

        for(vector<MapPoint*>::iterator vitMP=vpMapPointsKFi.begin(), vendMP=vpMapPointsKFi.end(); vitMP!=vendMP; vitMP++)
        {
            MapPoint* pMP = *vitMP;
            if(!pMP)
                continue;
            if(pMP->isBad() || pMP->mnFuseCandidateForKF == mpCurrentKeyFrame->mnId)
                continue;
            pMP->mnFuseCandidateForKF = mpCurrentKeyFrame->mnId;
            vpFuseCandidates.push_back(pMP);
        }
    }

    matcher.Fuse(mpCurrentKeyFrame,vpFuseCandidates);


    // Update points
    vpMapPointMatches = mpCurrentKeyFrame->GetMapPointMatches();
    for(size_t i=0, iend=vpMapPointMatches.size(); i<iend; i++)
    {
        MapPoint* pMP=vpMapPointMatches[i];
        if(pMP)
        {
            if(!pMP->isBad())
            {
                pMP->ComputeDistinctiveDescriptors();
                pMP->UpdateNormalAndDepth();
            }
        }
    }

    // Update connections in covisibility graph
    mpCurrentKeyFrame->UpdateConnections();
}

cv::Mat LocalMapping::ComputeF12(KeyFrame *&pKF1, KeyFrame *&pKF2)
{
    cv::Mat R1w = pKF1->GetRotation();
    cv::Mat t1w = pKF1->GetTranslation();
    cv::Mat R2w = pKF2->GetRotation();
    cv::Mat t2w = pKF2->GetTranslation();

    cv::Mat R12 = R1w*R2w.t();
    cv::Mat t12 = -R1w*R2w.t()*t2w+t1w;

    cv::Mat t12x = SkewSymmetricMatrix(t12);

    const cv::Mat &K1 = pKF1->mK;
    const cv::Mat &K2 = pKF2->mK;


    return K1.t().inv()*t12x*R12*K2.inv();
}

void LocalMapping::RequestStop()
{
    unique_lock<mutex> lock(mMutexStop);
    mbStopRequested = true;
    unique_lock<mutex> lock2(mMutexNewKFs);
    mbAbortBA = true;
}

bool LocalMapping::Stop()
{
    unique_lock<mutex> lock(mMutexStop);
    if(mbStopRequested && !mbNotStop)
    {
        mbStopped = true;
        cout << "Local Mapping STOP" << endl;
        return true;
    }

    return false;
}

bool LocalMapping::isStopped()
{
    unique_lock<mutex> lock(mMutexStop);
    return mbStopped;
}

bool LocalMapping::stopRequested()
{
    unique_lock<mutex> lock(mMutexStop);
    return mbStopRequested;
}

void LocalMapping::Release()
{
    unique_lock<mutex> lock(mMutexStop);
    unique_lock<mutex> lock2(mMutexFinish);
    if(mbFinished)
        return;
    mbStopped = false;
    mbStopRequested = false;
    for(list<KeyFrame*>::iterator lit = mlNewKeyFrames.begin(), lend=mlNewKeyFrames.end(); lit!=lend; lit++)
        delete *lit;
    mlNewKeyFrames.clear();

    cout << "Local Mapping RELEASE" << endl;
}

bool LocalMapping::AcceptKeyFrames()
{
    unique_lock<mutex> lock(mMutexAccept);
    return mbAcceptKeyFrames;
}

void LocalMapping::SetAcceptKeyFrames(bool flag)
{
    unique_lock<mutex> lock(mMutexAccept);
    mbAcceptKeyFrames=flag;
}

bool LocalMapping::SetNotStop(bool flag)
{
    unique_lock<mutex> lock(mMutexStop);

    if(flag && mbStopped)
        return false;

    mbNotStop = flag;

    return true;
}

void LocalMapping::InterruptBA()
{
    mbAbortBA = true;
}

void LocalMapping::KeyFrameCulling()
{
    // Check redundant keyframes (only local keyframes)
    // A keyframe is considered redundant if the 90% of the MapPoints it sees, are seen
    // in at least other 3 keyframes (in the same or finer scale)
    // We only consider close stereo points
    vector<KeyFrame*> vpLocalKeyFrames = mpCurrentKeyFrame->GetVectorCovisibleKeyFrames();

    for(vector<KeyFrame*>::iterator vit=vpLocalKeyFrames.begin(), vend=vpLocalKeyFrames.end(); vit!=vend; vit++)
    {
        KeyFrame* pKF = *vit;
        if(pKF->mnId==0)
            continue;
        const vector<MapPoint*> vpMapPoints = pKF->GetMapPointMatches();

        int nObs = 3;
        const int thObs=nObs;
        int nRedundantObservations=0;
        int nMPs=0;
        for(size_t i=0, iend=vpMapPoints.size(); i<iend; i++)
        {
            MapPoint* pMP = vpMapPoints[i];
            if(pMP)
            {
                if(!pMP->isBad())
                {
                    if(!mbMonocular)
                    {
                        if(pKF->mvDepth[i]>pKF->mThDepth || pKF->mvDepth[i]<0)
                            continue;
                    }

                    nMPs++;
                    if(pMP->Observations()>thObs)
                    {
                        const int &scaleLevel = pKF->mvKeysUn[i].octave;
                        const map<KeyFrame*, size_t> observations = pMP->GetObservations();
                        int nObs=0;
                        for(map<KeyFrame*, size_t>::const_iterator mit=observations.begin(), mend=observations.end(); mit!=mend; mit++)
                        {
                            KeyFrame* pKFi = mit->first;
                            if(pKFi==pKF)
                                continue;
                            const int &scaleLeveli = pKFi->mvKeysUn[mit->second].octave;

                            if(scaleLeveli<=scaleLevel+1)
                            {
                                nObs++;
                                if(nObs>=thObs)
                                    break;
                            }
                        }
                        if(nObs>=thObs)
                        {
                            nRedundantObservations++;
                        }
                    }
                }
            }
        }  

        if(nRedundantObservations>0.9*nMPs)
            pKF->SetBadFlag();
    }
}

cv::Mat LocalMapping::SkewSymmetricMatrix(const cv::Mat &v)
{
    return (cv::Mat_<float>(3,3) <<             0, -v.at<float>(2), v.at<float>(1),
            v.at<float>(2),               0,-v.at<float>(0),
            -v.at<float>(1),  v.at<float>(0),              0);
}

void LocalMapping::RequestReset()
{
    {
        unique_lock<mutex> lock(mMutexReset);
        mbResetRequested = true;
    }

    while(1)
    {
        {
            unique_lock<mutex> lock2(mMutexReset);
            if(!mbResetRequested)
                break;
        }
        usleep(3000);
    }
}

void LocalMapping::ResetIfRequested()
{
    unique_lock<mutex> lock(mMutexReset);
    if(mbResetRequested)
    {
        mlNewKeyFrames.clear();
        mlpRecentAddedMapPoints.clear();
        mbResetRequested=false;
    }
}

void LocalMapping::RequestFinish()
{
    unique_lock<mutex> lock(mMutexFinish);
    mbFinishRequested = true;
}

bool LocalMapping::CheckFinish()
{
    unique_lock<mutex> lock(mMutexFinish);
    return mbFinishRequested;
}

void LocalMapping::SetFinish()
{
    unique_lock<mutex> lock(mMutexFinish);
    mbFinished = true;    
    unique_lock<mutex> lock2(mMutexStop);
    mbStopped = true;
}

bool LocalMapping::isFinished()
{
    unique_lock<mutex> lock(mMutexFinish);
    return mbFinished;
}

} //namespace ORB_SLAM
