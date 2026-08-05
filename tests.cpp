#include <asl/Array.h>
#include <asl/Matrix.h>
#include <asl/Matrix3.h>
#include <asl/Matrix4.h>
#include <asl/Vec3.h>
#include <asl/String.h>
#include <asl/points.h>
#include <asl/pointsR.h>
#include <asl/testing.h>

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
