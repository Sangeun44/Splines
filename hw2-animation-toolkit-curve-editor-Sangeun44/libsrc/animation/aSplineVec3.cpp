#include "aSplineVec3.h"
#include <algorithm>
#include <Eigen/Dense>

#pragma warning(disable:4018)
#pragma warning(disable:4244)

ASplineVec3::ASplineVec3() : mInterpolator(new ALinearInterpolatorVec3())
{
}

ASplineVec3::~ASplineVec3()
{
  delete mInterpolator;
}

void ASplineVec3::setFramerate(double fps)
{
  mInterpolator->setFramerate(fps);
}

double ASplineVec3::getFramerate() const
{
  return mInterpolator->getFramerate();
}

void ASplineVec3::setLooping(bool loop)
{
  mLooping = loop;
}

bool ASplineVec3::getLooping() const
{
  return mLooping;
}

void ASplineVec3::setInterpolationType(ASplineVec3::InterpolationType type)
{
  double fps = getFramerate();

  delete mInterpolator;
  switch (type)
  {
  case LINEAR: mInterpolator = new ALinearInterpolatorVec3();
    break;
  case CUBIC_BERNSTEIN: mInterpolator = new ABernsteinInterpolatorVec3();
    break;
  case CUBIC_CASTELJAU: mInterpolator = new ACasteljauInterpolatorVec3();
    break;
  case CUBIC_MATRIX: mInterpolator = new AMatrixInterpolatorVec3();
    break;
  case CUBIC_HERMITE: mInterpolator = new AHermiteInterpolatorVec3();
    break;
  case CUBIC_BSPLINE: mInterpolator = new ABSplineInterpolatorVec3();
    break;
  };

  mInterpolator->setFramerate(fps);
  computeControlPoints();
  cacheCurve();
}

ASplineVec3::InterpolationType ASplineVec3::getInterpolationType() const
{
  return mInterpolator->getType();
}

void ASplineVec3::editKey(int keyID, const vec3& value)
{
  assert(keyID >= 0 && keyID < mKeys.size());
  mKeys[keyID].second = value;
  computeControlPoints();
  cacheCurve();
}

void ASplineVec3::editControlPoint(int ID, const vec3& value)
{
  assert(ID >= 0 && ID < mCtrlPoints.size()+2);
  if (ID == 0)
  {
    mStartPoint = value;
    computeControlPoints();
  }
  else if (ID == mCtrlPoints.size() + 1)
  {
    mEndPoint = value;
    computeControlPoints();
  }
  else mCtrlPoints[ID - 1] = value;
  cacheCurve();
}

void ASplineVec3::appendKey(double time, const vec3& value, bool updateCurve)
{
  mKeys.push_back(Key(time, value));

  if (mKeys.size() >= 2)
  {
    int totalPoints = mKeys.size();

    //If there are more than 1 interpolation point, set up the 2 end points to help determine the curve.
    //They lie on the tangent of the first and last interpolation points.
    vec3 tmp = mKeys[0].second - mKeys[1].second;
    double n = tmp.Length();
    mStartPoint = mKeys[0].second + (tmp / n) * n * 0.25;
    // distance to endpoint is 25% of distance between first 2 points

    tmp = mKeys[totalPoints - 1].second - mKeys[totalPoints - 2].second;
    n = tmp.Length();
    mEndPoint = mKeys[totalPoints - 1].second + (tmp / n) * n * 0.25;
  }

  if (updateCurve)
  {
    computeControlPoints();
    cacheCurve();
  }
}

void ASplineVec3::appendKey(const vec3& value, bool updateCurve)
{
  if (mKeys.size() == 0)
  {
    appendKey(0, value, updateCurve);
  }
  else
  {
    double lastT = mKeys[mKeys.size() - 1].first;
    appendKey(lastT + 1, value, updateCurve);
  }
}

void ASplineVec3::deleteKey(int keyID)
{
  assert(keyID >= 0 && keyID < mKeys.size());
  mKeys.erase(mKeys.begin() + keyID);
  computeControlPoints();
  cacheCurve();
}

vec3 ASplineVec3::getKey(int keyID)
{
  assert(keyID >= 0 && keyID < mKeys.size());
  return mKeys[keyID].second;
}

int ASplineVec3::getNumKeys() const
{
  return mKeys.size();
}

vec3 ASplineVec3::getControlPoint(int ID)
{
  assert(ID >= 0 && ID < mCtrlPoints.size()+2);
  if (ID == 0) return mStartPoint;
  else if (ID == mCtrlPoints.size() + 1) return mEndPoint;
  else return mCtrlPoints[ID - 1];
}

int ASplineVec3::getNumControlPoints() const
{
  return mCtrlPoints.size() + 2; // include endpoints
}

void ASplineVec3::clear()
{
  mKeys.clear();
}

double ASplineVec3::getDuration() const
{
  return mKeys[mKeys.size() - 1].first;
}

double ASplineVec3::getNormalizedTime(double t) const
{
  return (t / getDuration());
}

vec3 ASplineVec3::getValue(double t)
{
  if (mCachedCurve.size() == 0) return vec3();

  double dt = mInterpolator->getDeltaTime();
  int rawi = (int)(t / dt); // assumes uniform spacing
  int i = rawi % mCachedCurve.size();
  double frac = t - rawi * dt;
  int inext = i + 1;
  if (!mLooping) inext = std::min<int>(inext, mCachedCurve.size() - 1);
  else inext = inext % mCachedCurve.size();

  vec3 v1 = mCachedCurve[i];
  vec3 v2 = mCachedCurve[inext];
  vec3 v = v1 * (1 - frac) + v2 * frac;
  return v;
}

void ASplineVec3::cacheCurve()
{
  mInterpolator->interpolate(mKeys, mCtrlPoints, mCachedCurve);
}

void ASplineVec3::computeControlPoints()
{
  mInterpolator->computeControlPoints(mKeys, mCtrlPoints, mStartPoint, mEndPoint);
}

int ASplineVec3::getNumCurveSegments() const
{
  return mCachedCurve.size();
}

vec3 ASplineVec3::getCurvePoint(int i) const
{
  return mCachedCurve[i];
}

//---------------------------------------------------------------------
AInterpolatorVec3::AInterpolatorVec3(ASplineVec3::InterpolationType t) : mDt(1.0 / 120.0), mType(t)
{
}

void AInterpolatorVec3::setFramerate(double fps)
{
  mDt = 1.0 / fps;
}

double AInterpolatorVec3::getFramerate() const
{
  return 1.0 / mDt;
}

double AInterpolatorVec3::getDeltaTime() const
{
  return mDt;
}

void AInterpolatorVec3::interpolate(const std::vector<ASplineVec3::Key>& keys,
                                    const std::vector<vec3>& ctrlPoints, std::vector<vec3>& curve)
{
  vec3 val = 0.0;
  double u = 0.0;

  curve.clear();

  int numSegments = keys.size() - 1;
  for (int segment = 0; segment < numSegments; segment++)
  {
    for (double t = keys[segment].first; t < keys[segment + 1].first - FLT_EPSILON; t += mDt)
    {
      // dTODO: Compute u, fraction of duration between segment and segmentnext, for example,
      // u = 0.0 when t = keys[segment-1].first  
      // u = 1.0 when t = keys[segment].first
		double frac = (t - keys[segment].first) / (keys[segment + 1].first - keys[segment].first);
		u = frac;

	//linear interpolate with u as lerp percentage
      val = interpolateSegment(keys, ctrlPoints, segment, u);
      curve.push_back(val);
    }
  }
  // add last point
  if (keys.size() > 1)
  {
    u = 1.0;
    val = interpolateSegment(keys, ctrlPoints, numSegments - 1, u);
    curve.push_back(val);
  }
}

vec3 lerp(vec3 p0, vec3 p1, double lerp) {
	vec3 final1 = p0 + (p1 - p0) * lerp;
	return final1;
}

vec3 ALinearInterpolatorVec3::interpolateSegment(
  const std::vector<ASplineVec3::Key>& keys, 
	const std::vector<vec3>& ctrlPoints, int segment, double u)
{
  vec3 curveValue(0, 0, 0);
  vec3 key0 = keys[segment].second;
  vec3 key1 = keys[segment + 1].second;

  // dTODO: 
  //Step 1: Create a Lerp helper function
  //Step 2: Linear interpolate between key0 and key1 so that u = 0 returns key0 and u = 1 returns key1
  curveValue = lerp(key0, key1, u);

  return curveValue;
}

vec3 ABernsteinInterpolatorVec3::interpolateSegment(
  const std::vector<ASplineVec3::Key>& keys,
  const std::vector<vec3>& ctrlPoints,
  int segment, double u)
{
  vec3 b0;
  vec3 b1;
  vec3 b2;
  vec3 b3;
  vec3 curveValue(0, 0, 0);
  // TODO: 
  // Step1: Get the 4 control points, b0, b1, b2 and b3 from the ctrlPoints vector
  // Step2: Compute the interpolated value f(u) point using  Bernstein polynomials
  //4 control points
  b0 = ctrlPoints[0 + segment * 4];
  b1 = ctrlPoints[1 + segment * 4];
  b2 = ctrlPoints[2 + segment * 4];
  b3 = ctrlPoints[3 + segment * 4];

  //compute the f(u) point using Bernstein

  vec3 final2 = b0 + (3 * b1 - 3 * b0) * u + (3 * b2 - 6 * b1 + 3 * b0) * u * u + (b3 - 3*b2 + 3*b1 - b0) * u * u * u;

  curveValue = final2;
  return curveValue;
}


vec3 ACasteljauInterpolatorVec3::interpolateSegment(
  const std::vector<ASplineVec3::Key>& keys,
  const std::vector<vec3>& ctrlPoints,
  int segment, double u)
{
  vec3 b0;
  vec3 b1;
  vec3 b2;
  vec3 b3;
  vec3 curveValue(0, 0, 0);

  // TODO: 
  // Step1: Get the 4 control points, b0, b1, b2 and b3 from the ctrlPoints vector
  // Step2: Compute the interpolated value f(u) point using  deCsteljau alogithm
  b0 = ctrlPoints[0 + segment * 4];
  b1 = ctrlPoints[1 + segment * 4];
  b2 = ctrlPoints[2 + segment * 4];
  b3 = ctrlPoints[3 + segment * 4];

  vec3 b01 = lerp(b0, b1, u);
  vec3 b11 = lerp(b1, b2, u);
  vec3 b12 = lerp(b2, b3, u);

  vec3 b20 = lerp(b01, b12, u);
  vec3 b21 = lerp(b11, b12, u);

  vec3 b30 = lerp(b20, b21, u);
	
  curveValue = b30;

  return curveValue;
}

vec3 AMatrixInterpolatorVec3::interpolateSegment(
  const std::vector<ASplineVec3::Key>& keys,
  const std::vector<vec3>& ctrlPoints,
  int segment, double u)
{
  vec3 b0;
  vec3 b1;
  vec3 b2;
  vec3 b3;
  vec3 curveValue(0, 0, 0);

  // TODO: 
  // Step1: Get the 4 control points, b0, b1, b2 and b3 from the ctrlPoints vector
  // Step2: Compute the interpolated value f(u) point using  matrix method f(u) = GMU
  // Hint: Using Eigen::MatrixXd data representations for a matrix operations
  b0 = ctrlPoints[0 + segment * 4];
  b1 = ctrlPoints[1 + segment * 4];
  b2 = ctrlPoints[2 + segment * 4];
  b3 = ctrlPoints[3 + segment * 4];

  // x and y 3x4
  Eigen::MatrixXf matB(3, 4);
  matB(0, 0) = b0[0]; matB(0, 1) = b1[0]; matB(0, 2) = b2[0]; matB(0, 3) = b3[0];
  matB(1, 0) = b0[1]; matB(1, 1) = b1[1]; matB(1, 2) = b2[1]; matB(1, 3) = b3[1];
  matB(2, 0) = b0[2]; matB(2, 1) = b1[2]; matB(2, 2) = b2[2]; matB(2, 3) = b3[2];

  // 4x4
  Eigen::Matrix4f matG(4, 4);
  matG(0, 0) = 1; matG(0, 1) = -3; matG(0, 2) = 3; matG(0, 3) = -1;
  matG(1, 0) = 0; matG(1, 1) = 3; matG(1, 2) = -6; matG(1, 3) = 3;
  matG(2, 0) = 0; matG(2, 1) = 0; matG(2, 2) = 3; matG(2, 3) = -3;
  matG(3, 0) = 0; matG(3, 1) = 0; matG(3, 2) = 0; matG(3, 3) = 1;

  Eigen::Vector4f vecU(4);
  vecU(0) = 1.;
  vecU(1) = u;
  vecU(2) = u * u;
  vecU(3) = u * u * u;
	  
  Eigen::Vector3f vecfinal = matB * matG * vecU;
  curveValue = vec3(vecfinal(0), vecfinal(1), vecfinal(2));
  return curveValue;
}

vec3 AHermiteInterpolatorVec3::interpolateSegment(
  const std::vector<ASplineVec3::Key>& keys,
  const std::vector<vec3>& ctrlPoints,
  int segment, double u)
{
  vec3 p0 = keys[segment].second;
  vec3 p1 = keys[segment + 1].second;
  vec3 q0 = ctrlPoints[segment]; // slope at p0
  vec3 q1 = ctrlPoints[segment + 1]; // slope at p1
  vec3 curveValue(0, 0, 0);

  // TODO: Compute the interpolated value h(u) using a cubic Hermite polynomial  
  double u3 = u * u * u;
  double u2 = u * u;

  vec3 b1 = vec3(
				(2 * u3 - 3 * u2 + 1) * p0[0],
				(2 * u3 - 3 * u2 + 1) * p0[1],
				(2 * u3 - 3 * u2 + 1) * p0[2]
	  );
  vec3 b2 = vec3(
		(-2 * u3 + 3 * u2) * p1[0],
	  (-2 * u3 + 3 * u2) * p1[1],
	  (-2 * u3 + 3 * u2) * p1[2]
	  );
  vec3 b3 = vec3(
	  (u3 - 2 * u2 + u) * q0[0],
	  (u3 - 2 * u2 + u) * q0[1],
	  (u3 - 2 * u2 + u) * q0[2]
	  );
  vec3 b4 = vec3(
	  (u3 - u2) * q1[0],
	  (u3 - u2) * q1[1],
	  (u3 - u2) * q1[2]
	  );

  vec3 hu = b1 + b2 + b3 + b4;
  curveValue = hu;
  return curveValue;
}


vec3 ABSplineInterpolatorVec3::interpolateSegment(
  const std::vector<ASplineVec3::Key>& keys,
  const std::vector<vec3>& ctrlPoints,
  int segment, double u)
{
  vec3 curveValue(0, 0, 0);

  // Hint: Create a recursive helper function N(knots,n,j,t) to calculate BSpline basis function values at t, where
  //     knots = knot array
  //	   n = degree of the spline curves (n =3 for cubic)
  //     j = curve interval on knot vector in which to interpolate
  //     t = time value	

  // Step 1: determine the index j
  // Step 2: compute the n nonzero Bspline Basis functions N given j
  // Step 3: get the corresponding control points from the ctrlPoints vector
  // Step 4: compute the Bspline curveValue at time t


  return curveValue;
}

void ACubicInterpolatorVec3::computeControlPoints(
  const std::vector<ASplineVec3::Key>& keys,
  std::vector<vec3>& ctrlPoints,
  vec3& startPoint, vec3& endPoint)
{
  ctrlPoints.clear();
  if (keys.size() <= 1) return;

  for (int i = 0; i < keys.size() - 1; i++)
  {
	  vec3 b0, b1, b2, b3;
	  // TODO: compute b0, b1, b2, b3

	 if (keys.size() == 2) {
		 b0 = keys[0].second;
		 b3 = keys[1].second;

		 vec3 s0 = (b3 - startPoint) / 2.;
		 vec3 s1 = (endPoint - b0) / 2.;

		 b1 = b0 + s0 / 3.;
		 b2 = b3 - s1 / 3.;

	 } else if (i == 0) {
		b0 = keys[i].second;
		b3 = keys[i + 1].second;

		vec3 s0 = (b3 - startPoint) / 2.;// keys[i + 1].second;
		vec3 s1 = (keys[i + 2].second - b0) / 2.;

		b1 = b0 + s0 / 3.;
		b2 = b3 - s1 / 3.;
		
	} else if (i == keys.size() - 2) {
		b0 = keys[i].second;
		b3 = keys[i + 1].second;

		vec3 s0 = (b3 - keys[i - 1].second) / 2.;// keys[i + 1].second;
		vec3 s1 = (endPoint - b0) / 2.;

		b1 = b0 + s0 / 3.;
		b2 = b3 - s1 / 3.;

	}  else {
		b0 = keys[i].second;
		b3 = keys[i + 1].second;

		vec3 s0 = (b3 - keys[i - 1].second) / 2.;// keys[i + 1].second;
		vec3 s1 = (keys[i + 2].second - b0) / 2.;

		b1 = b0 + s0 / 3.;
		b2 = b3 - s1 / 3.;
	 }

    ctrlPoints.push_back(b0);
    ctrlPoints.push_back(b1);
    ctrlPoints.push_back(b2);
    ctrlPoints.push_back(b3);
  }
}

void AHermiteInterpolatorVec3::computeControlPoints(
  const std::vector<ASplineVec3::Key>& keys,
  std::vector<vec3>& ctrlPoints,
  vec3& startPoint, vec3& endPoint)
{
  ctrlPoints.clear();
  if (keys.size() <= 1) return;

  int numKeys = keys.size();

  // TODO: 
  // For each key point pi, compute the corresonding value of the slope pi_prime.

  // Hints: Using Eigen::MatrixXd for a matrix data structures, 
  // this can be accomplished by solving the system of equations AC=D for C.

  // Don't forget to save the values computed for C in ctrlPoints

  // For clamped endpoint conditions, 
  // set 1st derivative at first and last points (p0 and pm) to s0 and s1, respectively
  // For natural endpoints, 
  // set 2nd derivative at first and last points (p0 and pm) equal to 0

  // Step 1: Initialize A
  // Step 2: Initialize D
  // Step 3: Solve AC=D for C
  // Step 4: Save control points in ctrlPoints

  Eigen::MatrixXd matA(numKeys, numKeys);
  Eigen::MatrixXd matC(numKeys, 3);

  Eigen::MatrixXd matD(numKeys, 3);

  for (int i = 0; i < numKeys; i++) {
	  for (int j = 0; j < numKeys; j++) {
		  matA(i, j) = 0;
	  }
  }
  //matA initialization
  for (int i = 0; i < numKeys; i++) {
	  for (int j = 0; j < numKeys; j++) {
		if (i == 0 && j == 0) {
			 matA(i, j) = 2;
		} else if (i == 0 && j == 1) {
			 matA(i, j) = 1;
		} else if (i == numKeys - 1 && j == numKeys - 1) {
			matA(i, j) = 2;
		} else if (i == numKeys - 1 && j == numKeys - 2) {
			matA(i, j) = 1;
		} else if(i > 0 && i < numKeys - 1) {
			if (j == i) {
				matA(i, j) = 4;
			} else if (i == j + 1) {
				matA(i, j) = 1;
			} else if (i == j - 1) {
				matA(i, j) = 1;
			} 
		} 
	  }
  }

  //matD initialization
  //get s1 and s0
  vec3 p1 = keys[1].second;
  vec3 p0 = keys[0].second;
  vec3 s0 = vec3(3 * (p1[0] - p0[0]),
				3 * (p1[1] - p0[1]),
				3 * (p1[2] - p0[2])
  );

  vec3 pnm1 = keys[numKeys - 1].second;
  vec3 pnm2 = keys[numKeys - 2].second;
  vec3 s1 = 
	  vec3(3 * (pnm1[0] - pnm2[0]),
			3 * (pnm1[1] - pnm2[1]),
			3 * (pnm1[2] - pnm2[2])
  );

  matD(0, 0) = s0[0]; 
  matD(0, 1) = s0[1]; 
  matD(0, 2) = s0[2];
  matD(numKeys - 1, 0) = s1[0]; 
  matD(numKeys - 1, 1) = s1[1];
  matD(numKeys - 1, 2) = s1[2];

  //initalize matD row 2 to n-2
  for (int i = 1; i < numKeys - 1; i++) {
	  vec3 pi = keys[i + 1].second;
	  vec3 pnm1 = keys[i - 1].second;
	  vec3 result = 
			vec3(3 * (pi[0] - pnm1[0]),
				3 * (pi[1] - pnm1[1]),
				3 * (pi[2] - pnm1[2])
	  );
	  matD(i, 0) = result[0];
	  matD(i, 1) = result[1];
	  matD(i, 2) = result[2];
  }

  Eigen::MatrixXd matAInv = matA.inverse();
  Eigen::MatrixXd finalMat = matAInv * matD;

  for (int i = 0; i < numKeys; i++) {
	  vec3 newCoord = vec3(finalMat(i,0), finalMat(i, 1), finalMat(i, 2));
	  ctrlPoints.push_back(newCoord);
  }

}


void ABSplineInterpolatorVec3::computeControlPoints(
  const std::vector<ASplineVec3::Key>& keys,
  std::vector<vec3>& ctrlPoints,
  vec3& startPt, vec3& endPt)
{
  ctrlPoints.clear();
  if (keys.size() <= 1) return;

  // TODO: c
  // Hints: 
  // 1. use Eigen::MatrixXd to calculate the control points by solving the system of equations AC=D for C

  // 2. Create a recursive helper function dN(knots,n,t,l) to calculate derivative BSpline values at t, where
  //     knots = knot array
  //	   n = degree of the spline curves (n =3 for cubic)
  //     j = interval on knot vector in which to interpolate
  //     t = time value
  //     l = derivative (l = 1 => 1st derivative)

  // Step 1: Calculate knot vector using a uniform BSpline
  //         (assume knots are evenly spaced 1 apart and the start knot is at time = 0.0)

  // Step 2: Calculate A matrix  for a natural BSpline
  //         (Set 2nd derivative at t0 and tm to zero, where tm is the last point knot; m = #segments)

  // Step 3: Calculate  D matrix composed of our target points to interpolate

  // Step 4: Solve AC=D for C 

  // Step 5: save control points in ctrlPoints
}
