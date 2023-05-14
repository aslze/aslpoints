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

	return (A.pseudoinverse() * b).data();
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

	Matrix_<T> p1 = A.pseudoinverse() * b;
	Vec3_<T>   n1 = Vec3_<T>(1, -p1[0], -p1[1]);

	for (int i = 0; i < points.length(); i++)
	{
		A(i, 0) = points[i].x;
		A(i, 1) = points[i].z;
		A(i, 2) = 1;
		b(i, 0) = points[i].y;
	}

	Matrix_<T> p2 = A.pseudoinverse() * b;

	Vec3_<T> n2 = Vec3_<T>(-p2[0], 1, -p2[1]);

	for (int i = 0; i < points.length(); i++)
	{
		A(i, 0) = points[i].x;
		A(i, 1) = points[i].y;
		A(i, 2) = 1;
		b(i, 0) = points[i].z;
	}

	Matrix_<T> p3 = A.pseudoinverse() * b;
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

	// optionally refine with distances to plane ?

	return { c.x, c.y, c.z, n.x, n.y, n.z };
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

	asl::Matrix_<T> x = A.pseudoinverse() * b;
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

	asl::Matrix_<T> x = A.pseudoinverse() * b;
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
	asl::Matrix_<T> x = points1.length() == 4 ? solve(A, b) : A.pseudoinverse() * b;
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
		for (int j = 0; j < deg + 1; j++)
			A(i, j) = pow(p[i].x, T(j));
		b[i] = p[i].y;
	}

	return A.pseudoinverse() * b;
}

/**
 * Computes the coefficients of multivariate polynomial that best fits the points p
 * \return Array of coefficients in a0 + a1*x + a2*y + a3*x^2 + a4*y^2 + a5
 */
template<class T>
asl::Matrix_<T> fitPoly(const asl::Array<asl::Vec3_<T>>& p, int deg = 2)
{
	asl::Matrix_<T> A(p.length(), 2 + deg * 2);
	asl::Matrix_<T> b(p.length());

	for (int i = 0; i < p.length(); i++)
	{
		A(i, 0) = 1;
		for (int j = 0; j < deg; j++)
		{
			A(i, 1 + j * 2) = pow(p[i].x, T(j + 1));
			A(i, 2 + j * 2) = pow(p[i].y, T(j + 1));
		}
		A(i, deg * 2 + 1) = p[i].x * p[i].y;
		b(i, 0) = p[i].z;
	}

	return A.pseudoinverse() * b;
}

}

#endif
