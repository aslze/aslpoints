// Copyright(c) 1999-2023 aslze
// Licensed under the MIT License (http://opensource.org/licenses/MIT)

#ifndef ASL_POINTS_H
#define ASL_POINTS_H

#include <asl/Array.h>
#include <asl/Matrix.h>
#include <asl/Matrix3.h>
#include <asl/Matrix4.h>
#include <asl/Vec3.h>

namespace asl
{
/**
 * Fits a circle to a set of points and returns its center and radius as Vec3(cx, cy, r)
 */
template<class T>
asl::Vec3_<T> fitCircle(const asl::Array<asl::Vec2_<T>>& points)
{
	if (points.length() == 2)
	{
		asl::Vec2_<T> c = (points[0] + points[1]) / 2;
		return asl::Vec3_<T>(c.x, c.y, (points[0] - c).length());
	}
	asl::Matrix_<T> A(points.length(), 3);
	asl::Matrix_<T> b(points.length(), 1);
	for (int i = 0; i < points.length(); i++)
	{
		A(i, 0) = 2 * points[i].x;
		A(i, 1) = 2 * points[i].y;
		A(i, 2) = 1;
		b(i, 0) = sqr(points[i].x) + sqr(points[i].y);
	}

	asl::Matrix_<T> a = solve(A, b);

	return asl::Vec3_<T>(a[0], a[1], sqrt(a[2] + sqr(a[0]) + sqr(a[1])));
}

/**
 * Fits a plane to set of 3D points, in the form z = a*x + b*y + c
 * \return the parameters [a, b, c]
 */
template<class T>
asl::Array<T> fitPlaneXY(const asl::Array<asl::Vec3_<T>>& points)
{
	asl::Matrix_<T> A(points.length(), 3);
	asl::Matrix_<T> b(points.length(), 1);
	for (int i = 0; i < points.length(); i++)
	{
		A(i, 0) = points[i].x;
		A(i, 1) = points[i].y;
		A(i, 2) = 1;
		b(i, 0) = points[i].z;
	}

	return solve(A, b).data();
}

/**
 * Fits a plane to a set of 3D points (returns a point and a normal [x, y, z, nx, ny, nz])
 */
template<class T>
asl::Array<T> fitPlane(const asl::Array<asl::Vec3_<T>>& points)
{
	Matrix_<T> A(points.length(), 3);
	Matrix_<T> b(points.length(), 1);

	for (int i = 0; i < points.length(); i++)
	{
		A(i, 0) = points[i].y;
		A(i, 1) = points[i].z;
		A(i, 2) = 1;
		b(i, 0) = points[i].x;
	}

	Matrix_<T> p1 = solve(A, b);
	Vec3_<T>   n1 = Vec3_<T>(1, -p1[0], -p1[1]);

	for (int i = 0; i < points.length(); i++)
	{
		A(i, 0) = points[i].x;
		A(i, 1) = points[i].z;
		A(i, 2) = 1;
		b(i, 0) = points[i].y;
	}

	Matrix_<T> p2 = solve(A, b);

	Vec3_<T> n2 = Vec3_<T>(-p2[0], 1, -p2[1]);

	for (int i = 0; i < points.length(); i++)
	{
		A(i, 0) = points[i].x;
		A(i, 1) = points[i].y;
		A(i, 2) = 1;
		b(i, 0) = points[i].z;
	}

	Matrix_<T> p3 = solve(A, b);
	Vec3_<T>   n3 = Vec3_<T>(-p3[0], -p3[1], 1);

	T        s1 = 0, s2 = 0, s3 = 0;
	Vec3_<T> c(0, 0, 0);

	for (const auto& p : points)
	{
		s1 += sqr(p1[0] * p.y + p1[1] * p.z + p1[2] - p.x);
		s2 += sqr(p2[0] * p.x + p2[1] * p.z + p2[2] - p.y);
		s3 += sqr(p3[0] * p.x + p3[1] * p.y + p3[2] - p.z);
		c += p;
	}

	c /= (T)points.length();

	Vec3_<T> n;
	Vec3_<T> o;

	if (fabs(s1) < fabs(s2) && fabs(s1) < fabs(s3))
	{
		n = n1.normalized();
		o = Vec3_<T>(p1[2], 0, 0);
	}
	else if (fabs(s2) < fabs(s1) && fabs(s2) < fabs(s3))
	{
		n = n2.normalized();
		o = Vec3_<T>(0, p2[2], 0);
	}
	else
	{
		n = n3.normalized();
		o = Vec3_<T>(0, 0, p3[2]);
	}

	s1 = 0;
	for (int i = 0; i < points.length(); i++)
		s1 += (n * (points[i] - c) * n).length2();

	c -= n * (c - o) * n;

	s1 = 0;
	for (int i = 0; i < points.length(); i++)
		s1 += (n * (points[i] - c) * n).length2();

	return { c.x, c.y, c.z, n.x, n.y, n.z };
}

template<class T>
Pair<Vec3_<T>> getOrthonormalBase(const asl::Vec3_<T>& v)
{
	Vec3_<T> a(1, 0, 0);
	if (fabs(a * v.normalized()) > T(0.8))
		a = { 0, 1, 0 };
	return { (a ^ v).normalized(), ((a ^ v) ^ v).normalized() };
}

/**
 * Fits a circle to a set of 3D points and returns its center and its plane normal, with length equal to its radius
 */
template <class T>
static Pair<Vec3_<T> > fitCircle(const Array<Vec3_<T> >& points)
{
	auto     plane = fitPlane(points);
	Vec3_<T> pbase(plane[0], plane[1], plane[2]);
	Vec3_<T> normal(plane[3], plane[4], plane[5]);

	auto base = getOrthonormalBase(normal);

	auto points2 = points.map_<Vec2_<T>>(
	    [=](const Vec3_<T>& p) { return Vec2_<T>(base.first * (p - pbase), base.second * (p - pbase)); });

	auto circle = fitCircle(points2);

	auto center = pbase + circle.x * base.first + circle.y * base.second;

	return { center, normal * circle.z };
}


/**
 * Estimates the rigid transform between two sets of 2D points
 */
template<class T>
asl::Matrix3_<T> findRigidTransform(const asl::Array<asl::Vec2_<T>>& points1, const asl::Array<asl::Vec2_<T>>& points2)
{
	asl::Matrix_<T> A(points1.length() * 2, 4);
	asl::Matrix_<T> b(points1.length() * 2, 1);
	for (int i = 0; i < points1.length(); i++)
	{
		A(2 * i, 0) = points1[i].x;
		A(2 * i, 1) = points1[i].y;
		A(2 * i, 2) = 1;
		A(2 * i, 3) = 0;
		A(2 * i + 1, 0) = points1[i].y;
		A(2 * i + 1, 1) = -points1[i].x;
		A(2 * i + 1, 2) = 0;
		A(2 * i + 1, 3) = 1;
		b(2 * i, 0) = points2[i].x;
		b(2 * i + 1, 0) = points2[i].y;
	}

	asl::Matrix_<T> x = solve(A, b);
	return asl::Matrix3_<T>(x[0], x[1], x[2], //
	                        -x[1], x[0], x[3]);
}

/**
 * Computes the rigid transform between two sets of 3D points
 */
template<class T>
asl::Matrix4_<T> findRigidTransform(const asl::Array<asl::Vec3_<T>>& points1, const asl::Array<asl::Vec3_<T>>& points2)
{
	asl::Matrix_<T> A(points1.length() * 3, 12, T(0));
	asl::Matrix_<T> b(points1.length() * 3);
	for (int i = 0; i < points1.length(); i++)
	{
		A(3 * i, 0) = points1[i].x;
		A(3 * i, 1) = points1[i].y;
		A(3 * i, 2) = points1[i].z;
		A(3 * i, 3) = 1;
		A(3 * i + 1, 4) = points1[i].x;
		A(3 * i + 1, 5) = points1[i].y;
		A(3 * i + 1, 6) = points1[i].z;
		A(3 * i + 1, 7) = 1;
		A(3 * i + 2, 8) = points1[i].x;
		A(3 * i + 2, 9) = points1[i].y;
		A(3 * i + 2, 10) = points1[i].z;
		A(3 * i + 2, 11) = 1;
		b(3 * i, 0) = points2[i].x;
		b(3 * i + 1, 0) = points2[i].y;
		b(3 * i + 2, 0) = points2[i].z;
	}

	asl::Matrix_<T> x = solve(A, b);
	return orthonormalize(asl::Matrix4_<T>(x[0], x[1], x[2], x[3], //
	                                       x[4], x[5], x[6], x[7], //
	                                       x[8], x[9], x[10], x[11]));
}

/**
 * Computes a homogrphy (perspective transform) between two sets of 4 2D points. Points from source to dest can be
 * computed as: p2 = (H * Vec3(p1, 1)).h2c();
 */
template<class T>
asl::Matrix3_<T> findHomography(const asl::Array<asl::Vec2_<T>>& points1, const asl::Array<asl::Vec2_<T>>& points2)
{
	if (points1.length() != points2.length() || points1.length() < 4)
		return asl::Matrix3_<T>::identity();
	asl::Matrix_<T> A(2 * points1.length() + 1, 9, T(0));
	asl::Matrix_<T> b(2 * points1.length() + 1, 1, T(0));
	for (int i = 0; i < points1.length(); i++)
	{
		A(2 * i, 0) = -points1[i].x;
		A(2 * i, 1) = -points1[i].y;
		A(2 * i, 2) = -1;
		A(2 * i, 6) = points1[i].x * points2[i].x;
		A(2 * i, 7) = points1[i].y * points2[i].x;
		A(2 * i, 8) = points2[i].x;
		A(2 * i + 1, 3) = -points1[i].x;
		A(2 * i + 1, 4) = -points1[i].y;
		A(2 * i + 1, 5) = -1;
		A(2 * i + 1, 6) = points1[i].x * points2[i].y;
		A(2 * i + 1, 7) = points1[i].y * points2[i].y;
		A(2 * i + 1, 8) = points2[i].y;
	}
	A(2 * points1.length(), 8) = 1;
	b(2 * points1.length(), 0) = 1;
	asl::Matrix_<T> x = solve(A, b);
	x *= 1 / x[8];
	return asl::Matrix3_<T>(x[0], x[1], x[2], //
	                        x[3], x[4], x[5], //
	                        x[6], x[7], x[8]);
}

/**
 * Computes the coefficients of polynomial of the given degree that best fits the points p
 * \return Array of coefficients in increasing exponent order [a0, a1, a2, ...] for a0 + a1*x + a2*x^2 ...
 */
template<class T>
asl::Matrix_<T> fitPoly(const asl::Array<asl::Vec2_<T>>& p, int deg = 1)
{
	asl::Matrix_<T> A(p.length(), deg + 1);
	asl::Matrix_<T> b(p.length());

	for (int i = 0; i < p.length(); i++)
	{
		T xx = 1;
		for (int j = 0; j < deg + 1; j++)
		{
			A(i, j) = xx;
			xx *= p[i].x;
		}
		b[i] = p[i].y;
	}

	return solve(A, b);
}

/**
 * Computes the coefficients of bivariate polynomial that best fits the points p
 * \return Array of coefficients [a00, a01, ..., a10, a11, ...] f(x,y) = sum(a[i,j] * x^i * y^j)
 */
template<class T>
asl::Matrix_<T> fitPoly(const asl::Array<asl::Vec3_<T>>& p, int deg = 2)
{
	asl::Matrix_<T> A(p.length(), (deg + 1) * (deg + 1));
	asl::Matrix_<T> b(p.length());

	for (int i = 0; i < p.length(); i++)
	{
		T xx = 1;
		for (int j = 0; j < deg + 1; j++)
		{
			T yy = 1;
			for (int k = 0; k < deg + 1; k++)
			{
				A(i, j * (deg + 1) + k) = xx * yy;
				yy *= p[i].y;
			}
			xx *= p[i].x;
		}
		b(i, 0) = p[i].z;
	}

	return solve(A, b);
}

/**
 * Evaluates polynomial at (x)
 */
template<class T>
T polynomial(const asl::Matrix_<T>& a, T x)
{
	int deg = a.length() - 1;
	T   z = 0, xx = 1;
	for (int j = 0; j < deg + 1; j++)
	{
		z += a[j] * xx;
		xx *= x;
	}

	return z;
}

/**
 * Evaluates a bivariate polynomial at (x, y)
 */
template<class T>
T polynomial(const asl::Matrix_<T>& a, T x, T y)
{
	int deg = (int)sqrt(a.length()) - 1;
	T   z = 0, xx = 1;
	for (int j = 0; j < deg + 1; j++)
	{
		T yy = 1;
		for (int k = 0; k < deg + 1; k++)
		{
			z += a[j * (deg + 1) + k] * xx * yy;
			yy *= y;
		}
		xx *= x;
	}

	return z;
}

}

#endif
