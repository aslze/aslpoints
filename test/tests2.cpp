#include <asl/pointsR.h>
#include <asl/testing.h>

#pragma warning(disable : 4244)

using namespace asl;

ASL_TEST_ENABLE_MAIN()

// Fixed data keeps the expected answers reproducible. RANSAC itself samples randomly.
// Outliers are well outside the threshold and make up about 10% of each data set.
ASL_TEST(fitPlaneXYR)
{
	Array<Vec3> points;
	for (int x = -5; x <= 5; x++)
		for (int y = -5; y <= 5; y++)
			points << Vec3(x, y, 2 * x - 3 * y + 4);
	for (int i = 0; i < 12; i++)
		points << Vec3(i - 5, i % 7 - 3, 100 + i * i);
	auto plane = fitPlaneXYR(points, 0.1f);
	ASL_ASSERT(plane.length() == 3);
	ASL_EXPECT_NEAR(plane[0], 2.f, 0.001f);
	ASL_EXPECT_NEAR(plane[1], -3.f, 0.001f);
	ASL_EXPECT_NEAR(plane[2], 4.f, 0.001f);
}

ASL_TEST(fitPlaneR)
{
	for (bool vertical : { false, true })
	{
		Array<Vec3> points;
		for (int x = -5; x <= 5; x++)
			for (int y = -5; y <= 5; y++)
				points << (vertical ? Vec3(5, x, y) : Vec3(x, y, 2 * x - 3 * y + 4));
		for (int i = 0; i < 12; i++)
			points << Vec3(30 + i, i % 7, 100 + i * i);
		auto plane = fitPlaneR(points, 0.1f);
		ASL_ASSERT(plane.length() == 6);
		Vec3 center(plane[0], plane[1], plane[2]);
		Vec3 normal(plane[3], plane[4], plane[5]);
		Vec3 expected = vertical ? Vec3(1, 0, 0) : Vec3(-2, 3, 1).normalized();
		ASL_EXPECT_NEAR(center, vertical ? Vec3(5, 0, 0) : Vec3(0, 0, 4), 0.001f);
		ASL_EXPECT_NEAR(normal.length(), 1.f, 0.001f);
		ASL_EXPECT_NEAR(float(fabs(normal * expected)), 1.f, 0.001f);
	}
}

ASL_TEST(fitCircleR)
{
	Array<Vec2> points;
	for (int i = 0; i < 100; i++)
	{
		float t = float(2 * PI * i / 100);
		points << Vec2(7 + 5 * cos(t), -3 + 5 * sin(t));
	}
	for (int i = 0; i < 10; i++)
		points << Vec2(30 + i, -40 - i * i);
	ASL_EXPECT_NEAR(fitCircleR(points, 0.1f), Vec3(7, -3, 5), 0.001f);
}

ASL_TEST(fitEllipseR)
{
	for (float angle : { 0.f, 0.6f })
	{
		Array<Vec2> points;
		for (int i = 0; i < 100; i++)
		{
			float t = float(2 * PI * i / 100);
			float x = 5 * cos(t), y = 2 * sin(t);
			points << Vec2(3 + x * cos(angle) - y * sin(angle), -4 + x * sin(angle) + y * cos(angle));
		}
		for (int i = 0; i < 10; i++)
			points << Vec2(30 + i, -40 - i * i);
		auto ellipse = fitEllipseR(points, 0.1f);
		ASL_ASSERT(ellipse.length() == 5);
		ASL_EXPECT_NEAR(ellipse[0], 3.f, 0.001f);
		ASL_EXPECT_NEAR(ellipse[1], -4.f, 0.001f);
		ASL_EXPECT_NEAR(ellipse[2], 5.f, 0.001f);
		ASL_EXPECT_NEAR(ellipse[3], 2.f, 0.001f);
		ASL_EXPECT_NEAR(float(fabs(cos(ellipse[4] - angle))), 1.f, 0.001f);
	}
}

ASL_TEST(fitSphereR)
{
	Vec3        center(3, -4, 7);
	Array<Vec3> points;
	for (int i = 0; i < 100; i++)
	{
		float z = -0.99f + 0.02f * i;
		float t = 2.399963f * i;
		float r = sqrt(1 - z * z);
		points << center + 6.f * Vec3(r * cos(t), r * sin(t), z);
	}
	for (int i = 0; i < 10; i++)
		points << Vec3(30 + i, -40 - i * i, 50 + i);
	auto sphere = fitSphereR(points, 0.1f);
	ASL_EXPECT_NEAR(sphere.xyz(), center, 0.001f);
	ASL_EXPECT_NEAR(sphere.w, 6.f, 0.001f);
}

ASL_TEST(fitLineR2D)
{
	for (Vec2 dir : { Vec2(0.6f, 0.8f), Vec2(0, 1), Vec2(1, 0) })
	{
		Array<Vec2> points;
		for (int i = -50; i <= 50; i++)
			points << Vec2(3, -4) + float(i) * dir;
		for (int i = 0; i < 10; i++)
			points << Vec2(100 + i, -100 - i * i);
		auto line = fitLineR(points, 0.1f);
		ASL_EXPECT_NEAR(line.first, Vec2(3, -4), 0.001f);
		ASL_EXPECT_NEAR(line.second.length(), 1.f, 0.001f);
		ASL_EXPECT_NEAR(float(fabs(line.second * dir)), 1.f, 0.001f);
	}
}

ASL_TEST(fitLineR3D)
{
	Vec3        origin(3, -4, 7), dir = Vec3(1, 2, 3).normalized();
	Array<Vec3> points;
	for (int i = -50; i <= 50; i++)
		points << origin + float(i) * dir;
	for (int i = 0; i < 10; i++)
		points << Vec3(100 + i, -100 - i * i, 50 + i);
	auto line = fitLineR(points, 0.1f);
	ASL_EXPECT_NEAR(line.first, origin, 0.001f);
	ASL_EXPECT_NEAR(line.second.length(), 1.f, 0.001f);
	ASL_EXPECT_NEAR(float(fabs(line.second * dir)), 1.f, 0.001f);
}

ASL_TEST(fitPolyR1D)
{
	Array<Vec2> points;
	for (int i = -50; i <= 50; i++)
	{
		float x = i / 20.f;
		points << Vec2(x, 2 - 3 * x + 4 * x * x - 0.5f * x * x * x);
	}
	for (int i = 0; i < 10; i++)
		points << Vec2(i / 5.f - 1, 100 + i * i);
	auto   poly = fitPolyR(points, 3, 0.1f);
	Matrix expected = { 2.f, -3.f, 4.f, -0.5f };
	ASL_ASSERT(poly.length() == 4);
	for (int i = 0; i < 4; i++)
		ASL_EXPECT_NEAR(poly[i], expected[i], 0.001f);
}

ASL_TEST(fitPolyR2D)
{
	Array<Vec3> points;
	for (int i = -5; i <= 5; i++)
		for (int j = -5; j <= 5; j++)
		{
			float x = i / 3.f, y = j / 3.f;
			points << Vec3(x, y,
			               1 + 2 * y + 3 * y * y - 4 * x + 5 * x * y - 6 * x * y * y + 7 * x * x - 8 * x * x * y +
			                   9 * x * x * y * y);
		}
	for (int i = 0; i < 12; i++)
		points << Vec3(i / 10.f - 0.5f, (i % 5) / 4.f, 500 + i * i);
	auto   poly = fitPolyR(points, 2, 0.1f);
	Matrix expected = { 1.f, 2.f, 3.f, -4.f, 5.f, -6.f, 7.f, -8.f, 9.f };
	ASL_ASSERT(poly.length() == 9);
	for (int i = 0; i < 9; i++)
		ASL_EXPECT_NEAR(poly[i], expected[i], 0.002f);
}

ASL_TEST(transformsR2D)
{
	Matrix3 affine(2, 0.5f, 3, -0.25f, 1.5f, -4);
	Matrix3 rigid(0.6f, -0.8f, 3, 0.8f, 0.6f, -4);
	for (bool isRigid : { false, true })
	{
		Matrix3     expected = isRigid ? rigid : affine;
		Array<Vec2> source, dest;
		for (int x = -5; x <= 5; x++)
			for (int y = -5; y <= 5; y++)
			{
				Vec2 p(x, y);
				source << p;
				dest << expected * p;
			}
		for (int i = 0; i < 12; i++)
		{
			source << Vec2(i / 3.f, i % 5);
			dest << Vec2(100 + i, -100 - i * i);
		}
		auto actual = isRigid ? findRigidTransformR(source, dest, 0.1f) : findAffineTransformR(source, dest, 0.1f);
		for (Vec2 p : { Vec2(0, 0), Vec2(1, 0), Vec2(0, 1), Vec2(7, -3) })
			ASL_EXPECT_NEAR(actual * p, expected * p, 0.001f);
	}
}

ASL_TEST(transformsR3D)
{
	Matrix4 affine(2, 0.5f, -0.2f, 3, -0.25f, 1.5f, 0.4f, -4, 0.3f, -0.1f, 0.8f, 5);
	Matrix4 rigid = Matrix4::translate(Vec3(3, -4, 5)) * Matrix4::rotate(Vec3(0.3f, -0.2f, 0.4f));
	for (bool isRigid : { false, true })
	{
		Matrix4     expected = isRigid ? rigid : affine;
		Array<Vec3> source, dest;
		for (int x = -2; x <= 2; x++)
			for (int y = -2; y <= 2; y++)
				for (int z = -2; z <= 2; z++)
				{
					Vec3 p(x, y, z);
					source << p;
					dest << expected * p;
				}
		for (int i = 0; i < 12; i++)
		{
			source << Vec3(i / 3.f, i % 5, i % 7);
			dest << Vec3(100 + i, -100 - i * i, 50 + i);
		}
		auto actual = isRigid ? findRigidTransformR(source, dest, 0.1f) : findAffineTransformR(source, dest, 0.1f);
		for (Vec3 p : { Vec3(0, 0, 0), Vec3(1, 0, 0), Vec3(0, 1, 0), Vec3(0, 0, 1), Vec3(7, -3, 2) })
			ASL_EXPECT_NEAR(actual * p, expected * p, 0.002f);
	}
}

ASL_TEST(findHomographyR)
{
	Matrix3     expected(1.2f, 0.3f, 2, -0.2f, 0.9f, -1, 0.04f, -0.03f, 1);
	Array<Vec2> source, dest;
	for (int x = -5; x <= 5; x++)
		for (int y = -5; y <= 5; y++)
		{
			Vec2 p(x, y);
			source << p;
			dest << (expected * Vec3(p, 1)).h2c();
		}
	for (int i = 0; i < 12; i++)
	{
		source << Vec2(i / 3.f, i % 5);
		dest << Vec2(100 + i, -100 - i * i);
	}
	auto actual = findHomographyR(source, dest, 0.1f);
	for (Vec2 p : { Vec2(0, 0), Vec2(1, 2), Vec2(-1, 3), Vec2(7, -3) })
		ASL_EXPECT_NEAR((actual * Vec3(p, 1)).h2c(), (expected * Vec3(p, 1)).h2c(), 0.001f);
}

ASL_TEST(emptyFitsR)
{
	Array<Vec2> p2;
	Array<Vec3> p3;
	ASL_EXPECT_NEAR(fitCircleR(p2, 0.1f), Vec3(0, 0, 0), 0.0001f);
	auto sphere = fitSphereR(p3, 0.1f);
	ASL_EXPECT_NEAR(sphere.xyz(), Vec3(0, 0, 0), 0.0001f);
	ASL_EXPECT_NEAR(sphere.w, 0.f, 0.0001f);
	auto l2 = fitLineR(p2, 0.1f);
	auto l3 = fitLineR(p3, 0.1f);
	ASL_EXPECT(l2.second.norm(), ==, 0.0f);
	ASL_EXPECT_NEAR(l3.second.norm(), 0.0f, 0.0001f);
	auto xy = fitPlaneXYR(p3, 0.1f);
	auto plane = fitPlaneR(p3, 0.1f);
	auto ellipse = fitEllipseR(p2, 0.1f);
	ASL_ASSERT(xy.length() == 3 && plane.length() == 6 && ellipse.length() == 5);
	for (int i = 0; i < 3; i++)
		ASL_EXPECT_NEAR(xy[i], 0.f, 0.0001f);
	for (int i = 0; i < 6; i++)
		ASL_EXPECT_NEAR(plane[i], 0.f, 0.0001f);
	for (int i = 0; i < 5; i++)
		ASL_EXPECT_NEAR(ellipse[i], 0.f, 0.0001f);
	auto poly1 = fitPolyR(p2, 3, 0.1f);
	auto poly2 = fitPolyR(p3, 2, 0.1f);
	ASL_ASSERT(poly1.length() == 4 && poly2.length() == 9);
	for (int i = 0; i < poly1.length(); i++)
		ASL_EXPECT_NEAR(poly1[i], 0.f, 0.0001f);
	for (int i = 0; i < poly2.length(); i++)
		ASL_EXPECT_NEAR(poly2[i], 0.f, 0.0001f);
}

ASL_TEST(invalidTransformsR)
{
	// Empty and unequal correspondence lists return null transforms.
	float t = 0.1f;
	for (int n : { 0, 2 })
	{
		Array<Vec2> p2(n, Vec2(1, 2)), q2;
		Array<Vec3> p3(n, Vec3(1, 2, 3)), q3;
		for (auto m : { findAffineTransformR(p2, q2, t), findRigidTransformR(p2, q2, t), findHomographyR(p2, q2, t) })
			ASL_EXPECT(m.norm(), ==, 0);
		for (auto m : { findAffineTransformR(p3, q3, t), findRigidTransformR(p3, q3, t) })
			ASL_EXPECT(m.norm(), ==, 0);
	}
}
