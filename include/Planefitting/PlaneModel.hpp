#pragma once

#include "AbstractModel.hpp"

typedef std::array<GRANSAC::VPFloat, 3> Vector3VP;
class Point3D
	: public GRANSAC::AbstractParameter
{
public:
	Point3D(GRANSAC::VPFloat x, GRANSAC::VPFloat y, GRANSAC::VPFloat z)
	{
		m_Point3D[0] = x;
		m_Point3D[1] = y;
		m_Point3D[2] = z;
	};

	Vector3VP m_Point3D;
};

class PlaneModel
	: public GRANSAC::AbstractModel<3>
{
protected:
	// Parametric form
	GRANSAC::VPFloat m_a, m_b, m_c, m_d; // ax + by +cz+d = 0
	GRANSAC::VPFloat m_DistDenominator; // = sqrt(a^2 + b^2+c^2). Stored for efficiency reasons

										// Another parametrization y = mx + d
	GRANSAC::VPFloat a; // Slope
	GRANSAC::VPFloat b; // Intercept
	GRANSAC::VPFloat d; // Intercept

    //std::ofstream debug;

	virtual GRANSAC::VPFloat ComputeDistanceMeasure(std::shared_ptr<GRANSAC::AbstractParameter> Param) override
	{
		auto ExtPoint3D = std::dynamic_pointer_cast<Point3D>(Param);
		if (ExtPoint3D == nullptr)
			throw std::runtime_error("Line2DModel::ComputeDistanceMeasure() - Passed parameter are not of type Point2D.");

		// Return distance between passed "point" and this line
		// http://mathworld.wolfram.com/Point-LineDistance2-Dimensional.html
		GRANSAC::VPFloat Numer = fabs(m_a * ExtPoint3D->m_Point3D[0] + m_b * ExtPoint3D->m_Point3D[1]
			+ m_c*	ExtPoint3D->m_Point3D[2] + m_d);
		GRANSAC::VPFloat Dist = Numer / m_DistDenominator;

		// Debug
		//debug << "Point: " << ExtPoint3D->m_Point3D[0] << ", " << ExtPoint3D->m_Point3D[1] << ", " << ExtPoint3D->m_Point3D[2] << " ";
		//debug << "Distance: " << Dist << std::endl<< std::endl ;

		return Dist;
	};


public:
	PlaneModel(const std::vector<std::shared_ptr<GRANSAC::AbstractParameter>> &InputParams)
	{
		Initialize(InputParams);
	};

	virtual void Initialize(const std::vector<std::shared_ptr<GRANSAC::AbstractParameter>> &InputParams) override
	{
		if (InputParams.size() != 3)
			throw std::runtime_error("Line2DModel - Number of input parameters does not match minimum number required for this model.");

		// Check for AbstractParamter types
		auto Point1 = std::dynamic_pointer_cast<Point3D>(InputParams[0]);
		auto Point2 = std::dynamic_pointer_cast<Point3D>(InputParams[1]);
		auto Point3 = std::dynamic_pointer_cast<Point3D>(InputParams[2]);
		if (Point1 == nullptr || Point2 == nullptr || Point3 == nullptr)
			throw std::runtime_error("Line2DModel - InputParams type mismatch. It is not a Point2D.");

		std::copy(InputParams.begin(), InputParams.end(), m_MinModelParams.begin());

		// Compute the line parameters
		//m_m = (Point2->m_Point2D[1] - Point1->m_Point2D[1]) / (Point2->m_Point2D[0] - Point1->m_Point2D[0]); // Slope
		//m_d = Point1->m_Point2D[1] - m_m * Point1->m_Point2D[0]; // Intercept
																 // m_d = Point2->m_Point2D[1] - m_m * Point2->m_Point2D[0]; // Intercept - alternative should be the same as above

		// Compute the plane parameters    ax +by+ d = z    ax+by+cz+1=0;
		double x0 = Point1->m_Point3D[0];
		double y0 = Point1->m_Point3D[1];
		double z0 = Point1->m_Point3D[2];

		double x1 = Point2->m_Point3D[0];
		double y1 = Point2->m_Point3D[1];
		double z1 = Point2->m_Point3D[2];

		double x2 = Point3->m_Point3D[0];
		double y2 = Point3->m_Point3D[1];
		double z2 = Point3->m_Point3D[2];
		a = (z0*(y1 - y2)) / (x0*y1 - x1*y0 - x0*y2 + x2*y0 + x1*y2 - x2*y1) - (z1*(y0 - y2)) / (x0*y1 - x1*y0 - x0*y2 + x2*y0 + x1*y2 - x2*y1) + (z2*(y0 - y1)) / (x0*y1 - x1*y0 - x0*y2 + x2*y0 + x1*y2 - x2*y1);
		b = (z1*(x0 - x2)) / (x0*y1 - x1*y0 - x0*y2 + x2*y0 + x1*y2 - x2*y1) - (z0*(x1 - x2)) / (x0*y1 - x1*y0 - x0*y2 + x2*y0 + x1*y2 - x2*y1) - (z2*(x0 - x1)) / (x0*y1 - x1*y0 - x0*y2 + x2*y0 + x1*y2 - x2*y1);
		d = (z2*(x0*y1 - x1*y0)) / (x0*y1 - x1*y0 - x0*y2 + x2*y0 + x1*y2 - x2*y1) - (z1*(x0*y2 - x2*y0)) / (x0*y1 - x1*y0 - x0*y2 + x2*y0 + x1*y2 - x2*y1) + (z0*(x1*y2 - x2*y1)) / (x0*y1 - x1*y0 - x0*y2 + x2*y0 + x1*y2 - x2*y1);
		m_a = a/d;
		m_b = b/d;
		m_c = -1.0/d;
		m_d = d/d;

		m_DistDenominator = sqrt(m_a * m_a + m_b * m_b + m_c*m_c ); // Cache square root for efficiency
		m_PlaneCoefs[0] = m_a;
		m_PlaneCoefs[1] = m_b;
		m_PlaneCoefs[2] = m_c;
		m_PlaneCoefs[3] = m_d;

		m_PointCenter[0] = (x0 + x1 + x2) / 3;
		m_PointCenter[1] = (y0 + y1 + y2) / 3;
		m_PointCenter[2] = (z0 + z1 + z2) / 3;
	};
public:
	GRANSAC::VPFloat m_PlaneCoefs[4] = { 0 };
	GRANSAC::VPFloat m_PointCenter[3] = { 0 };


	virtual std::pair<GRANSAC::VPFloat, std::vector<std::shared_ptr<GRANSAC::AbstractParameter>>> Evaluate(const std::vector<std::shared_ptr<GRANSAC::AbstractParameter>>& EvaluateParams, GRANSAC::VPFloat Threshold)
	{
		std::vector<std::shared_ptr<GRANSAC::AbstractParameter>> Inliers;
		
		int nTotalParams = EvaluateParams.size();
		int nInliers = 0;
		//Inliers.reserve(nTotalParams);
        //debug.open( "debug.txt" ,ios_base::ate);
        //debug << "Plane: " << m_a << " x + " << m_b << " y + "  << m_c << "z + " <<m_d<< " = 0"<< std::endl;
		for (auto& Param : EvaluateParams)
		{
			if (ComputeDistanceMeasure(Param) < Threshold)
			{
				Inliers.push_back(Param);
				nInliers++;
			}
		}
		GRANSAC::VPFloat InlierFraction = GRANSAC::VPFloat(nInliers) / GRANSAC::VPFloat(nTotalParams); // This is the inlier fraction
		return std::make_pair(InlierFraction, Inliers);
	};
};

