// Copyright(c) 2023-2026 aslze
// Licensed under the MIT License (http://opensource.org/licenses/MIT)

#pragma once
#include "points.h"

namespace asl
{
template<class T, class TP, class R, class F, class F2>
R fitR(const Array<TP>& points, T threshold, int ms, F fit, F2 dist, R deflt)
{
	Array<TP> sample;
	R         best;
	int       mi = 0;
	if (points.length() == 0)
		return deflt;

	Random random;

	for (int i = 0; i < 200; i++)
	{
		sample.clear();
		for (int j = 0; j < ms; j++)
			sample << points[random(points.length() - 1)];
		auto model = fit(sample);
		int  ni = 0;
		for (int j = 0; j < points.length(); j++)
		{
			if (dist(model, points[j]) < threshold)
				ni++;
		}
		if (ni > mi)
		{
			mi = ni;
			best = model;
		}
		if ((i > 10 && ni > points.length() * 0.9) || (i > 2 && ni > points.length() * 0.995))
			break;
	}

	Array<TP> inliers;
	for (int i = 0; i < points.length(); i++)
	{
		if (dist(best, points[i]) < threshold)
			inliers << points[i];
	}

	if (mi < points.length() * 0.10)
		best = fit(points);
	else
		best = fit(inliers);

	return best;
}

template<class T, class TP, class R, class F, class F2>
R fit2R(const Array<TP>& points1, const Array<TP>& points2, T threshold, int ms, F fit, F2 dist, R deflt)
{
	Array<TP> sample1, sample2;
	R         best;
	int       mi = 0;
	if (points1.length() == 0 || points1.length() != points2.length())
		return deflt;

	Random random;

	for (int i = 0; i < 200; i++)
	{
		sample1.clear();
		sample2.clear();
		for (int j = 0; j < ms; j++)
		{
			int k = random(points1.length() - 1);
			sample1 << points1[k];
			sample2 << points2[k];
		}
		auto model = fit(sample1, sample2);
		int  ni = 0;
		for (int j = 0; j < points1.length(); j++)
		{
			if (dist(model, points1[j], points2[j]) < threshold)
				ni++;
		}
		if (ni > mi)
		{
			mi = ni;
			best = model;
		}
		if ((i > 10 && ni > points1.length() * 0.9) || (i > 2 && ni > points1.length() * 0.995))
			break;
	}

	Array<TP> inliers1, inliers2;
	inliers1.reserve(mi);
	inliers2.reserve(mi);
	for (int i = 0; i < points1.length(); i++)
	{
		if (dist(best, points1[i], points2[i]) < threshold)
		{
			inliers1 << points1[i];
			inliers2 << points2[i];
		}
	}

	if (mi < points1.length() * 0.10)
		best = fit(points1, points2);
	else
		best = fit(inliers1, inliers2);

	return best;
}

/**
 * Fits a plane over the XY plane to a set of 3D points using RANSAC. The plane is represented as a 3D vector (a, b, c)
 * such that z = ax + by + c.
 */
template<class T>
Array<T> fitPlaneXYR(const Array<Vec3_<T>>& points, T threshold)
{
	auto r = fitR(
	    points, threshold, 3, [](const Array<Vec3_<T>>& sample) { return fitPlaneXY(sample); },
	    [](const Array<T>& plane, const Vec3_<T>& p) { return fabs(plane[0] * p.x + plane[1] * p.y + plane[2] - p.z); },
	    Array<T>{ T(0), T(0), T(0) });

	return r;
}

/**
 * Fits a plane to a set of 3D points using RANSAC and returns [cx, cy, cz, nx, ny, nz].
 */
template<class T>
Array<T> fitPlaneR(const Array<Vec3_<T>>& points, T threshold)
{
	auto r = fitR(
	    points, threshold, 3, [](const Array<Vec3_<T>>& sample) { return fitPlane(sample); },
	    [](const Array<T>& plane, const Vec3_<T>& p) {
		    return fabs((p - Vec3_<T>(plane[0], plane[1], plane[2])) * Vec3_<T>(plane[3], plane[4], plane[5]));
	    },
	    Array<T>{ T(0), T(0), T(0), T(0), T(0), T(1) });

	return r;
}

/**
 * Fits a polynomial to a set of 2D points using RANSAC. The polynomial is represented as a matrix of coefficients.
 */
template<class T>
Matrix_<T> fitPolyR(const Array<Vec2_<T>>& points, int deg, T threshold)
{
	auto r = fitR(
	    points, threshold, deg + 1, [deg](const Array<Vec2_<T>>& sample) { return fitPoly(sample, deg); },
	    [](const Matrix_<T>& poly, const Vec2_<T>& p) { return fabs(polynomial(poly, p.x) - p.y); },
	    Matrix_<T>(deg + 1, 1, T(0)));
	return r;
}

/**
 * Fits a polynomial surface to a set of 3D points using RANSAC. The polynomial is represented as a matrix of
 * coefficients.
 */
template<class T>
Matrix_<T> fitPolyR(const Array<Vec3_<T>>& points, int deg, T threshold)
{
	auto r = fitR(
	    points, threshold, 2 * (deg + 1) * (deg + 1),
	    [deg](const Array<Vec3_<T>>& sample) { return fitPoly(sample, deg); },
	    [](const Matrix_<T>& poly, const Vec3_<T>& p) { return fabs(polynomial(poly, p.x, p.y) - p.z); },
	    Matrix_<T>((deg + 1) * (deg + 1), 1, T(0)));
	return r;
}

/**
 * Fits a circle to a set of 2D points using RANSAC. The circle is represented as a 3D vector (x, y, radius).
 */
template<class T>
Vec3_<T> fitCircleR(const Array<Vec2_<T>>& points, T threshold)
{
	auto r = fitR(
	    points, threshold, 3, [](const Array<Vec2_<T>>& sample) { return fitCircle(sample); },
	    [](const Vec3_<T>& circle, const Vec2_<T>& p) {
		    return fabs(sqrt((p.x - circle.x) * (p.x - circle.x) + (p.y - circle.y) * (p.y - circle.y)) - circle.z);
	    },
	    Vec3_<T>(0, 0, 0));

	return r;
}

/**
 * Fits an ellipse to a set of 2D points using RANSAC. The ellipse is represented {cx, cy, a, b, angle}.
 */
template<class T>
Array<T> fitEllipseR(const Array<Vec2_<T>>& points, T threshold)
{
	auto r = fitR(
	    points, threshold, 5, [](const Array<Vec2_<T>>& sample) { return fitEllipse(sample); },
	    [](const Array<T>& e, const Vec2_<T>& p) {
		    auto q = (Matrix3::rotate(-e[4]) * Matrix3::translate(-e[0], -e[1])) * p;
		    T    r = sqrt(sqr(q.x) / sqr(e[2]) + sqr(q.y) / sqr(e[3]));
		    return fabs(1 - 1 / r) * sqrt(sqr(q.x) + sqr(q.y));
	    },
	    Array<T>{ T(0), T(0), T(0), T(0), T(0) });

	return r;
}

/**
 * Fits a sphere to a set of 3D points using RANSAC. The sphere is represented as a 4D vector (x, y, z, radius).
 */
template<class T>
Vec4_<T> fitSphereR(const Array<Vec3_<T>>& points, T threshold)
{
	auto r = fitR(
	    points, threshold, 4, [](const Array<Vec3_<T>>& sample) { return fitSphere(sample); },
	    [](const Vec4_<T>& sphere, const Vec3_<T>& p) {
		    return fabs(sqrt((p.x - sphere.x) * (p.x - sphere.x) + (p.y - sphere.y) * (p.y - sphere.y) +
		                     (p.z - sphere.z) * (p.z - sphere.z)) -
		                sphere.w);
	    },
	    Vec4_<T>(0, 0, 0, 0));
	return r;
}

/**
 * Fits a line to a set of 2D points using RANSAC. The line is represented as a pair of points (start and end).
 */
template<class T>
Pair<Vec2_<T>> fitLineR(const Array<Vec2_<T>>& points, T threshold)
{
	auto r = fitR(
	    points, threshold, 2, [](const Array<Vec2_<T>>& sample) { return fitLine(sample); },
	    [](const Pair<Vec2_<T>>& line, const Vec2_<T>& p) { return distancePointLine(p, line.first, line.second); },
	    Pair<Vec2_<T>>(Vec2_<T>(0, 0), Vec2_<T>(1, 0)));
	return r;
}

/**
 * Fits a line to a set of 3D points using RANSAC. The line is represented as a pair of points (start and end).
 */
template<class T>
Pair<Vec3_<T>> fitLineR(const Array<Vec3_<T>>& points, T threshold)
{
	auto r = fitR(
	    points, threshold, 2, [](const Array<Vec3_<T>>& sample) { return fitLine(sample); },
	    [](const Pair<Vec3_<T>>& line, const Vec3_<T>& p) { return distancePointLine(p, line.first, line.second); },
	    Pair<Vec3_<T>>(Vec3_<T>(0, 0, 0), Vec3_<T>(1, 0, 0)));

	return r;
}

/**
 * Finds a rigid transform between two sets of 2D points using RANSAC
 */
template<class T>
Matrix3_<T> findRigidTransformR(const Array<Vec2_<T>>& points1, const Array<Vec2_<T>>& points2, T threshold)
{
	auto r = fit2R(
	    points1, points2, threshold, 2,
	    [](const Array<Vec2_<T>>& sample1, const Array<Vec2_<T>>& sample2) {
		    return findRigidTransform(sample1, sample2);
	    },
	    [](const Matrix3_<T>& transform, const Vec2_<T>& p1, const Vec2_<T>& p2) { return (p2 - transform * p1).norm(); },
	    Matrix3_<T>::identity());
	return r;
}

/**
 * Finds a rigid transform between two sets of 3D points using RANSAC
 */
template<class T>
Matrix4_<T> findRigidTransformR(const Array<Vec3_<T>>& points1, const Array<Vec3_<T>>& points2, T threshold)
{
	auto r = fit2R(
	    points1, points2, threshold, 3,
	    [](const Array<Vec3_<T>>& sample1, const Array<Vec3_<T>>& sample2) {
		    return findRigidTransform(sample1, sample2);
	    },
	    [](const Matrix4_<T>& transform, const Vec3_<T>& p1, const Vec3_<T>& p2) { return (p2 - transform * p1).norm(); },
	    Matrix4_<T>::identity());
	return r;
}

/**
 * Finds an affine transform between two sets of 3D points using RANSAC
 */
template<class T>
Matrix4_<T> findAffineTransformR(const Array<Vec3_<T>>& points1, const Array<Vec3_<T>>& points2, T threshold)
{
	auto r = fit2R(
	    points1, points2, threshold, 3,
	    [](const Array<Vec3_<T>>& sample1, const Array<Vec3_<T>>& sample2) {
		    return findAffineTransform(sample1, sample2);
	    },
	    [](const Matrix4_<T>& transform, const Vec3_<T>& p1, const Vec3_<T>& p2) { return (p2 - transform * p1).norm(); },
	    Matrix4_<T>::identity());
	return r;
}

/**
 * Finds an affine transform between two sets of 2D points using RANSAC
 */
template<class T>
Matrix3_<T> findAffineTransformR(const Array<Vec2_<T>>& points1, const Array<Vec2_<T>>& points2, T threshold)
{
	auto r = fit2R(
	    points1, points2, threshold, 3,
	    [](const Array<Vec2_<T>>& sample1, const Array<Vec2_<T>>& sample2) {
		    return findAffineTransform(sample1, sample2);
	    },
	    [](const Matrix3_<T>& transform, const Vec2_<T>& p1, const Vec2_<T>& p2) { return (p2 - transform * p1).norm(); },
	    Matrix3_<T>::identity());
	return r;
}

/**
 * Finds an homography transform between two sets of 2D points using RANSAC
 */
template<class T>
Matrix3_<T> findHomographyR(const Array<Vec2_<T>>& points1, const Array<Vec2_<T>>& points2, T threshold)
{
	auto r = fit2R(
	    points1, points2, threshold, 3,
	    [](const Array<Vec2_<T>>& sample1, const Array<Vec2_<T>>& sample2) { return findHomography(sample1, sample2); },
	    [](const Matrix3_<T>& transform, const Vec2_<T>& p1, const Vec2_<T>& p2) { return (p2 - transform ^ p1).norm(); },
	    Matrix3_<T>::identity());
	return r;
}

}
