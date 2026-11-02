#include <asl/Array.h>
#include <asl/Matrix.h>
#include <asl/Matrix3.h>
#include <asl/Matrix4.h>
#include <asl/Vec3.h>
#include <asl/String.h>
#include <asl/points.h>
#include <asl/pointsR.h>
#include <asl/testing.h>

#pragma warning(disable : 4244)

using namespace asl;

ASL_TEST_ENABLE_MAIN()

inline Vec3 randomVec3(float k)
{
	return Vec3(asl::random(-k, k), asl::random(-k, k), asl::random(-k, k));
}

ASL_TEST(fiPlaneXY)
{
	Array<Vec3>  points;
	Array<float> planeXY = { 0.1f, -0.05f, 10.3f };
	for (int i = 0; i < 2000; i++)
	{
		float x = asl::random(-500.f, 500.f);
		float y = asl::random(-500.f, 500.f);
		float z = planeXY[0] * x + planeXY[1] * y + planeXY[2] + asl::random(-1.f, 1.f);
		points << Vec3(x, y, z);
	}

	auto plane = fitPlaneXY(points);

	ASL_EXPECT_NEAR(plane[0], planeXY[0], 0.01f);
	ASL_EXPECT_NEAR(plane[1], planeXY[1], 0.01f);
	ASL_EXPECT_NEAR(plane[2], planeXY[2], 0.1f);

	// add some outliers
	for (int i = 0; i < 20; i++)
		points << Vec3(asl::random(50.f, 500.f), asl::random(0.f, 500.f), asl::random(100.f, 1000.f));

	auto planeR = fitPlaneXYR(points, 2.0f);

	ASL_EXPECT_NEAR(planeR[0], planeXY[0], 0.01f);
	ASL_EXPECT_NEAR(planeR[1], planeXY[1], 0.01f);
	ASL_EXPECT_NEAR(planeR[2], planeXY[2], 0.1f);
}

ASL_TEST(fitSphere)
{
	Array<Vec3> points;
	Vec3        center = randomVec3(100);
	float       radius = asl::random(50.f, 100.f);
	for (int i = 0; i < 2000; i++)
	{
		Vec3 p = center + randomVec3(1.0f).normalized() * radius;
		points << p;
	}
	auto sphere = fitSphere(points);
	ASL_EXPECT_NEAR(sphere.xyz(), center, 0.02f);
	ASL_EXPECT_NEAR(sphere.w, radius, 0.01f);
	// add some outliers
	for (int i = 0; i < 20; i++)
		points << Vec3(asl::random(50.f, 500.f), asl::random(0.f, 500.f), asl::random(100.f, 1000.f));
	auto sphereR = fitSphereR(points, 2.0f);

	ASL_EXPECT_NEAR(sphereR.xyz(), center, 0.001f);
	ASL_EXPECT_NEAR(sphereR.w, radius, 0.01f);
}

ASL_TEST(fitLine)
{
	Array<Vec3> points;
	asl::Vec3   origin = randomVec3(100);
	asl::Vec3   dir = randomVec3(1.0f).normalized();

	for (int i = 0; i < 2000; i++)
	{
		points << origin + dir * asl::random(-500.f, 500.f) + randomVec3(1.0f);
	}

	auto line = fitLine(points);

	ASL_EXPECT(distancePointLine(origin, line.first, line.second), <, 1.0f);
	ASL_ASSERT(line.second.angle(dir) < deg2rad(1.f) || line.second.angle(dir) > deg2rad(179.f));

	for (int i = 0; i < 20; i++)
		points << randomVec3(1000.f);

	auto lineR = fitLineR(points, 2.0f);

	ASL_EXPECT(distancePointLine(origin, lineR.first, lineR.second), <, 1.0f);
	ASL_ASSERT(lineR.second.angle(dir) < deg2rad(1.f) || lineR.second.angle(dir) > deg2rad(179.f));
}

ASL_TEST(fitCircle2D)
{
	Array<Vec2> points;
	for (int i = 0; i < 32; i++)
	{
		float t = float(2 * PI * i / 32);
		points << Vec2(7 + 5 * cos(t), -3 + 5 * sin(t));
	}
	ASL_EXPECT_NEAR(fitCircle(points), Vec3(7, -3, 5), 0.001f);
	Array<Vec2> diameter = { Vec2(4, -7), Vec2(10, 1) };
	ASL_EXPECT_NEAR(fitCircle(diameter), Vec3(7, -3, 5), 0.001f);
}

ASL_TEST(fitEllipse)
{
	// Check both axis-aligned and rotated ellipses; orientation is modulo pi.
	for (float angle : { 0.0f, 0.6f })
	{
		Array<Vec2> points;
		for (int i = 0; i < 64; i++)
		{
			float t = float(2 * PI * i / 64);
			float x = 5 * cos(t), y = 2 * sin(t);
			points << Vec2(3 + x * cos(angle) - y * sin(angle), -4 + x * sin(angle) + y * cos(angle));
		}
		auto e = fitEllipse(points);
		ASL_ASSERT(e.length() == 5);
		ASL_EXPECT_NEAR(e[0], 3.f, 0.001f);
		ASL_EXPECT_NEAR(e[1], -4.f, 0.001f);
		ASL_EXPECT_NEAR(e[2], 5.f, 0.001f);
		ASL_EXPECT_NEAR(e[3], 2.f, 0.001f);
		ASL_EXPECT_NEAR(float(fabs(cos(e[4] - angle))), 1.0f, 0.001f);
	}
	for (int n = 0; n < 5; n++)
	{
		auto e = fitEllipse(Array<Vec2>(n, Vec2(1, 2)));
		ASL_ASSERT(e.length() == 5);
		for (int i = 0; i < 5; i++)
			ASL_EXPECT_NEAR(e[i], 0.f, 0.0001f);
	}
}

ASL_TEST(fitPlane)
{
	Array<Vec3> points;
	for (int x = -2; x <= 2; x++)
		for (int y = -2; y <= 2; y++)
			points << Vec3(x, y, 2 * x - 3 * y + 4);
	Vec3 expectedNormal = Vec3(-2, 3, 1).normalized();
	for (auto plane : { fitPlane(points), fitPlane3D(points), fitPlaneXYZ(points) })
	{
		ASL_ASSERT(plane.length() == 6);
		Vec3 center(plane[0], plane[1], plane[2]);
		Vec3 normal(plane[3], plane[4], plane[5]);
		ASL_EXPECT_NEAR(center, Vec3(0, 0, 4), 0.001f);
		ASL_EXPECT_NEAR(normal.length(), 1.f, 0.001f);
		ASL_EXPECT_NEAR(float(fabs(normal * expectedNormal)), 1.0f, 0.001f);
	}
	// A vertical plane cannot be represented as z = ax + by + c.
	Array<Vec3> vertical = { Vec3(5, -2, -3), Vec3(5, 2, -3), Vec3(5, 2, 3), Vec3(5, -2, 3) };
	auto        plane = fitPlane(vertical);
	ASL_EXPECT_NEAR(Vec3(plane[0], plane[1], plane[2]), Vec3(5, 0, 0), 0.001f);
	ASL_EXPECT_NEAR(float(fabs(plane[3])), 1.0f, 0.001f);
	ASL_EXPECT_NEAR(plane[4], 0.f, 0.001f);
	ASL_EXPECT_NEAR(plane[5], 0.f, 0.001f);
}

ASL_TEST(fitPlaneSmallInputs)
{
	for (int n = 0; n < 3; n++)
	{
		auto p = fitPlane(Array<Vec3>(n, Vec3(1, 2, 3)));
		ASL_ASSERT(p.length() == 6);
		for (int i = 0; i < 6; i++)
			ASL_EXPECT_NEAR(p[i], 0.f, 0.0001f);
	}
	Array<Vec3> triangle = { Vec3(1, 2, 3), Vec3(1, 3, 3), Vec3(2, 2, 3) };
	auto        p = fitPlane(triangle);
	ASL_EXPECT_NEAR(Vec3(p[0], p[1], p[2]), triangle[0], 0.001f);
	ASL_EXPECT_NEAR(Vec3(p[3], p[4], p[5]), Vec3(0, 0, 1), 0.001f);
	Array<Vec3> collinear = { Vec3(0, 0, 0), Vec3(1, 0, 0), Vec3(2, 0, 0) };
	p = fitPlane(collinear);
	ASL_EXPECT_NEAR(Vec3(p[3], p[4], p[5]), Vec3(0, 0, 1), 0.001f);
}

ASL_TEST(getOrthonormalBase)
{
	for (Vec3 v : { Vec3(4, 0, 0), Vec3(0, -3, 0), Vec3(1, 2, 3) })
	{
		auto base = getOrthonormalBase(v);
		ASL_EXPECT_NEAR(base.first.length(), 1.f, 0.0001f);
		ASL_EXPECT_NEAR(base.second.length(), 1.f, 0.0001f);
		ASL_EXPECT_NEAR(base.first * base.second, 0.f, 0.0001f);
		ASL_EXPECT_NEAR(base.first * v, 0.f, 0.0001f);
		ASL_EXPECT_NEAR(base.second * v, 0.f, 0.0001f);
	}
}

ASL_TEST(fitCircle3D)
{
	Vec3        center(3, -4, 7);
	Vec3        u = Vec3(1, 1, 0).normalized(), v(0, 0, 1);
	Array<Vec3> points;
	for (int i = 0; i < 32; i++)
	{
		float t = float(2 * PI * i / 32);
		points << center + 6.f * (u * float(cos(t)) + v * float(sin(t)));
	}
	auto circle = fitCircle(points);
	ASL_EXPECT_NEAR(circle.first, center, 0.001f);
	ASL_EXPECT_NEAR(circle.second.length(), 6.f, 0.001f);
	ASL_EXPECT_NEAR(float(fabs(circle.second.normalized() * (u ^ v))), 1.0f, 0.001f);
}

ASL_TEST(fitLine2D)
{
	for (Vec2 dir : { Vec2(0.6f, 0.8f), Vec2(0, 1), Vec2(1, 0) })
	{
		Array<Vec2> points;
		for (int i = -5; i <= 5; i++)
			points << Vec2(3, -4) + float(i) * dir;
		auto line = fitLine(points);
		ASL_EXPECT_NEAR(line.first, Vec2(3, -4), 0.001f);
		ASL_EXPECT_NEAR(line.second.length(), 1.f, 0.001f);
		ASL_EXPECT_NEAR(float(fabs(line.second * dir)), 1.0f, 0.001f);
	}
}

ASL_TEST(fitLineSmallInputs)
{
	for (int n = 0; n < 2; n++)
	{
		auto l2 = fitLine(Array<Vec2>(n, Vec2(3, 4)));
		auto l3 = fitLine(Array<Vec3>(n, Vec3(3, 4, 5)));
		ASL_EXPECT_NEAR(l2.first, Vec2(0, 0), 0.0001f);
		ASL_EXPECT_NEAR(l2.second, Vec2(1, 0), 0.0001f);
		ASL_EXPECT_NEAR(l3.first, Vec3(0, 0, 0), 0.0001f);
		ASL_EXPECT_NEAR(l3.second, Vec3(1, 0, 0), 0.0001f);
	}
	Array<Vec3> points = { Vec3(1, 2, 3), Vec3(1, 5, 7) };
	auto        line = fitLine(points);
	ASL_EXPECT_NEAR(line.first, points[0], 0.0001f);
	ASL_EXPECT_NEAR(line.second, Vec3(0, 0.6f, 0.8f), 0.0001f);
}

ASL_TEST(fitSphereTwoPoints)
{
	Array<Vec3> points = { Vec3(1, 2, 3), Vec3(1, 5, 7) };
	auto        sphere = fitSphere(points);
	ASL_EXPECT_NEAR(sphere.xyz(), Vec3(1, 3.5f, 5), 0.0001f);
	ASL_EXPECT_NEAR(sphere.w, 2.5f, 0.0001f);
}

ASL_TEST(distancePointLine)
{
	// Directions are unit vectors; projections may lie on either side of the origin.
	for (float t : { -4.f, 0.f, 7.f })
	{
		ASL_EXPECT_NEAR(distancePointLine(Vec2(2 + t, 8), Vec2(2, 3), Vec2(1, 0)), 5.f, 0.0001f);
		ASL_EXPECT_NEAR(distancePointLine(Vec2(2 + t, 3), Vec2(2, 3), Vec2(1, 0)), 0.f, 0.0001f);
		ASL_EXPECT_NEAR(distancePointLine(Vec3(2 + t, 6, 8), Vec3(2, 3, 4), Vec3(1, 0, 0)), 5.f, 0.0001f);
		ASL_EXPECT_NEAR(distancePointLine(Vec3(2 + t, 3, 4), Vec3(2, 3, 4), Vec3(1, 0, 0)), 0.f, 0.0001f);
	}
}

ASL_TEST(transforms2D)
{
	Array<Vec2> source = { Vec2(-2, -1), Vec2(3, -1), Vec2(3, 4), Vec2(-2, 4), Vec2(0, 0) };
	Matrix3     affine(2, 0.5f, 3, -0.25f, 1.5f, -4);
	Matrix3     rigid(0.6f, -0.8f, 3, 0.8f, 0.6f, -4);
	for (bool isRigid : { false, true })
	{
		Matrix3     expected = isRigid ? rigid : affine;
		Array<Vec2> dest;
		for (auto p : source)
			dest << (expected * Vec3(p, 1)).h2c();
		auto actual = isRigid ? findRigidTransform(source, dest) : findAffineTransform(source, dest);
		for (auto p : { Vec2(0, 0), Vec2(1, 0), Vec2(0, 1), Vec2(7, -3) })
			ASL_EXPECT_NEAR((actual * Vec3(p, 1)).h2c(), (expected * Vec3(p, 1)).h2c(), 0.001f);
	}
}

ASL_TEST(transforms3D)
{
	Array<Vec3> source = { Vec3(0, 0, 0), Vec3(2, 0, 0), Vec3(0, 3, 0), Vec3(0, 0, 4), Vec3(1, -2, 3) };
	Matrix4     affine(2, 0.5f, -0.2f, 3, -0.25f, 1.5f, 0.4f, -4, 0.3f, -0.1f, 0.8f, 5);
	Matrix4     rigid = Matrix4::translate(Vec3(3, -4, 5)) * Matrix4::rotate(Vec3(0.3f, -0.2f, 0.4f));
	for (int mode = 0; mode < 3; mode++)
	{
		Matrix4     expected = mode == 0 ? affine : rigid;
		Array<Vec3> dest;
		for (auto p : source)
			dest << expected * p;
		auto actual = mode == 0   ? findAffineTransform(source, dest)
		              : mode == 1 ? findRigidTransform(source, dest)
		                          : findRigidTransform(source, dest, 0);
		for (auto p : { Vec3(0, 0, 0), Vec3(1, 0, 0), Vec3(0, 1, 0), Vec3(0, 0, 1), Vec3(7, -3, 2) })
			ASL_EXPECT_NEAR(actual * p, expected * p, 0.002f);
	}
	Array<Vec3> shortDest = { Vec3(1, 2, 3) };
	for (auto m : { findAffineTransform(source, shortDest), findRigidTransform(source, shortDest) })
		for (auto p : source)
			ASL_EXPECT_NEAR(m * p, p, 0.0001f);
}

ASL_TEST(findHomography)
{
	Matrix3     expected(1.2f, 0.3f, 2, -0.2f, 0.9f, -1, 0.04f, -0.03f, 1);
	Array<Vec2> source = { Vec2(-2, -1), Vec2(3, -1), Vec2(3, 4), Vec2(-2, 4) };
	Array<Vec2> dest;
	for (auto p : source)
		dest << (expected * Vec3(p, 1)).h2c();
	for (int n = 4; n <= 5; n++)
	{
		auto h = findHomography(source, dest);
		for (auto p : { Vec2(0, 0), Vec2(1, 2), Vec2(-1, 3) })
			ASL_EXPECT_NEAR((h * Vec3(p, 1)).h2c(), (expected * Vec3(p, 1)).h2c(), 0.001f);
		source << Vec2(1, 1);
		dest << (expected * Vec3(1, 1, 1)).h2c();
	}
	Array<Vec2> few = { Vec2(0, 0), Vec2(1, 0), Vec2(0, 1) };
	for (auto h : { findHomography(few, few), findHomography(source, few) })
		ASL_EXPECT_NEAR(h * Vec3(2, 3, 1), Vec3(2, 3, 1), 0.0001f);
}

ASL_TEST(polynomial1D)
{
	Matrix      coefficients = { 2.f, -3.f, 4.f, -0.5f };
	Array<Vec2> points, linear;
	for (int x = -3; x <= 3; x++)
	{
		float y = 2 - 3 * x + 4 * x * x - 0.5f * x * x * x;
		points << Vec2(x, y);
		linear << Vec2(x, 5 - 2 * x);
		ASL_EXPECT_NEAR(polynomial(coefficients, float(x)), y, 0.001f);
	}
	auto fitted = fitPoly(points, 3);
	ASL_ASSERT(fitted.length() == 4);
	for (int i = 0; i < 4; i++)
		ASL_EXPECT_NEAR(fitted[i], coefficients[i], 0.001f);
	auto line = fitPoly(linear);
	ASL_ASSERT(line.length() == 2);
	ASL_EXPECT_NEAR(line[0], 5.f, 0.001f);
	ASL_EXPECT_NEAR(line[1], -2.f, 0.001f);
	Matrix constant = { 7.f };
	ASL_EXPECT_NEAR(polynomial(constant, -2.f), 7.f, 0.001f);
}

ASL_TEST(polynomial2D)
{
	Matrix      coefficients = { 1.f, 2.f, 3.f, -4.f, 5.f, -6.f, 7.f, -8.f, 9.f };
	Array<Vec3> points;
	for (int x = -2; x <= 2; x++)
		for (int y = -2; y <= 2; y++)
		{
			float z =
			    1 + 2 * y + 3 * y * y - 4 * x + 5 * x * y - 6 * x * y * y + 7 * x * x - 8 * x * x * y + 9 * x * x * y * y;
			points << Vec3(x, y, z);
			ASL_EXPECT_NEAR(polynomial(coefficients, float(x), float(y)), z, 0.001f);
		}
	auto fitted = fitPoly(points);
	ASL_ASSERT(fitted.length() == 9);
	for (int i = 0; i < 9; i++)
		ASL_EXPECT_NEAR(fitted[i], coefficients[i], 0.001f);
	Matrix constant = { 7.f };
	ASL_EXPECT_NEAR(polynomial(constant, -2.f, 3.f), 7.f, 0.001f);
}

ASL_TEST(polynomial3D)
{
	// Distinct coefficients expose swapped variable order, including mixed terms.
	Matrix coefficients(27);
	for (int i = 0; i < 27; i++)
		coefficients[i] = float(i + 1);
	Array<Vec4> points;
	for (int x = -2; x <= 2; x++)
		for (int y = -2; y <= 2; y++)
			for (int z = -2; z <= 2; z++)
			{
				float w = 0;
				for (int i = 0; i < 3; i++)
					for (int j = 0; j < 3; j++)
						for (int k = 0; k < 3; k++)
							w += float((9 * i + 3 * j + k + 1) * pow(double(x), i) * pow(double(y), j) *
							           pow(double(z), k));
				points << Vec4(x, y, z, w);
				ASL_EXPECT_NEAR(polynomial(coefficients, float(x), float(y), float(z)), w, 0.001f);
			}
	auto fitted = fitPoly(points);
	ASL_ASSERT(fitted.length() == 27);
	for (int i = 0; i < 27; i++)
		ASL_EXPECT_NEAR(fitted[i], coefficients[i], 0.002f);
	Matrix constant = { 7.f };
	ASL_EXPECT_NEAR(polynomial(constant, -2.f, 3.f, 4.f), 7.f, 0.001f);
}
