// Copyright(c) 1999-2026 aslze
// Licensed under the MIT License (http://opensource.org/licenses/MIT)

#pragma once
#include "points.h"

namespace asl
{
template<class T, class TP, class R, class F, class F2>
R fitR(const asl::Array<TP>& points, T threshold, int ms, F fit, F2 dist, R deflt)
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

template<class T>
Array<T> fitPlaneXYR(const Array<Vec3_<T>>& points, T threshold)
{
	auto r = fitR(
	    points, threshold, 3, [](const Array<Vec3_<T>>& sample) { return fitPlaneXY(sample); },
	    [](const Array<T>& plane, const Vec3_<T>& p) { return fabs(plane[0] * p.x + plane[1] * p.y + plane[2] - p.z); },
	    Array<T>{ T(0), T(0), T(0) });

	return r;
}

template<class T>
Matrix_<T> fitPolyR(const Array<Vec2_<T>>& points, int deg, T threshold)
{
	auto r = fitR(
	    points, threshold, deg + 1, [deg](const Array<Vec2_<T>>& sample) { return fitPoly(sample, deg); },
	    [](const Matrix_<T>& poly, const Vec2_<T>& p) { return fabs(polynomial(poly, p.x) - p.y); },
	    Matrix_<T>(deg + 1, 1, T(0)));
	return r;
}

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

template<class T>
Pair<Vec2_<T>> fitLineR(const Array<Vec2_<T>>& points, T threshold)
{
	auto r = fitR(
	    points, threshold, 2, [](const Array<Vec2_<T>>& sample) { return fitLine(sample); },
	    [](const Pair<Vec2_<T>>& line, const Vec2_<T>& p) { return distancePointLine(p, line.first, line.second); },
	    Pair<Vec2_<T>>(Vec2_<T>(0, 0), Vec2_<T>(1, 0)));
	return r;
}

template<class T>
Pair<Vec3_<T>> fitLineR(const Array<Vec3_<T>>& points, T threshold)
{
	auto r = fitR(
	    points, threshold, 2, [](const Array<Vec3_<T>>& sample) { return fitLine(sample); },
	    [](const Pair<Vec3_<T>>& line, const Vec3_<T>& p) { return distancePointLine(p, line.first, line.second); },
	    Pair<Vec3_<T>>(Vec3_<T>(0, 0, 0), Vec3_<T>(1, 0, 0)));

	return r;
}

}
