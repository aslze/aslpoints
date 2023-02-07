#include <asl/Array.h>
#include <asl/Matrix.h>
#include <asl/Matrix3.h>
#include <asl/Matrix4.h>
#include <asl/Vec3.h>
#include <asl/String.h>
#include <asl/points.h>

using namespace asl;

int main()
{
	Array<Vec3> points;
	for (int i = 0; i < 1000; i++)
	{
		points << Matrix4::rotateX((float)PI / 4) * Vec3(asl::random(200.0f), asl::random(200.0f), asl::random(10.0f));
	}

	double t01 = now();

	auto plane = fitPlane(points);

	double t02 = now();

	printf("%s (t plane = %f)\n", *plane.join(", "), t02 - t01);

	Array<Vec2> points1 = points.map_<Vec2>([](Vec3 p) { return p.xy(); });

	auto tr = Matrix3::translate(50, 30) * Matrix3::rotate(0.5f);

	float k = 10.5f;

	Array<Vec2> points2 = points1.map([=](Vec2 p) { return (tr * p) + Vec2(asl::random(-k, k), asl::random(-k, k)); });

	double t1 = now();

	auto tr2 = findRigidTransform(points1, points2);

	double t2 = now();

	printf("t solve rigid 2D: %f s (%f)\n", t2 - t1, (tr2.inverse() * tr).trace() - 3);


	auto trb = Matrix4::translate(50, 30, 10) * Matrix4::rotateZ(0.5f) * Matrix4::rotateX(0.7f);

	Array<Vec3> points2b = points.map([=](Vec3 p) { return (trb * p) + 0.25f * Vec3(asl::random(-k, k), asl::random(-k, k), asl::random(-k, k)); });

	t1 = now();

	auto tr2b = findRigidTransform(points, points2b);

	t2 = now();

	printf("t solve rigid 3D: %f s (%f)\n", t2 - t1, (tr2b.inverse() * trb).norm() - 2);


	Array<Vec2> pts1, pts2;
	for (int i = 0; i < 2; i++)
		for (int j = 0; j < 2; j++)
		{
			pts1 << Vec2{ j * 100.0f + random(10.0f), i * 100.0f + random(30.0f) };
			pts2 << Vec2{ j * 200.0f + random(20.0f), i * 150.0f + random(50.0f) };
		}

	Matrix3 h = findHomography(pts1, pts2);

	for (int i = 0; i < 4; i++)
	{
		Vec2 p = (h * Vec3(pts1[i], 1)).h2c();
		printf("%f, %f -> %f, %f\n", p.x, p.y, pts2[i].x, pts2[i].y);
	}
}
