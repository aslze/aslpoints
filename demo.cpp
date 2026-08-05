#include <asl/Array.h>
#include <asl/Matrix.h>
#include <asl/Matrix3.h>
#include <asl/Matrix4.h>
#include <asl/Vec3.h>
#include <asl/String.h>
#include <asl/points.h>
#include <asl/pointsR.h>
#include <asl/TextFile.h>

using namespace asl;

int main()
{
	Random random;
	random.seed(0);
	float  k = 9.5f;
	double t1, t2;

	if (1)
	{
		float k = 5.0f;
		// fill Array <Vec2> points with random points around a line
		asl::Array<Vec2> points;
		asl::Vec2        dir = Vec2::polar(1, deg2rad(30.0f));
		for (int i = 0; i < 1000; i++)
		{
			points << Vec2(100, 50) + dir * random(-500.f, 500.f) + Vec2(asl::random(-k, k), asl::random(-k, k));
		}

		auto line = fitLine(points);

		printf("%f %f (t line = %f)\n", line.first.x, line.first.y, rad2deg(atan2(line.second.y, line.second.x)));

		auto line2 = fitLineR(points, 5.f);

		printf("%f %f (t lineR = %f)\n", line2.first.x, line2.first.y, rad2deg(atan2(line2.second.y, line2.second.x)));

		points << Vec2(0, 90900); // add an outlier

		line = fitLine(points);
		
		printf("%f %f (t line = %f)\n", line.first.x, line.first.y, rad2deg(atan2(line.second.y, line.second.x)));

		line2 = fitLineR(points, 5.f);

		printf("%f %f (t lineR = %f)\n", line2.first.x, line2.first.y, rad2deg(atan2(line2.second.y, line2.second.x)));

		// compute distance of points to that line
		float dmax = 0.0f;
		for (int i = 0; i < points.length(); i++)
		{
			float d = distancePointLine(points[i], line.first, line.second);
			if (d > dmax)
				dmax = d;
		}

		printf("max distance to line: %f\n", dmax);

		auto circle = fitCircleR(points, 5.0f);

		printf("circle: %f %f %f\n", circle.x, circle.y, circle.z);

		auto poly = fitPolyR(points, 3, 5.0f);

		return 0;
	}

#if 1
	if(1){
		float k = 1.0f;
		asl::random.seed(0);
		// fill Array <Vec2> points with random points around a line
		for (int j = 0; j < 60; j++)
		{
			asl::Array<Vec3> points;
			asl::Vec3        dir(asl::random(-1.0f, 1.0f), asl::random(-1.0f, 1.0f), asl::random(-1.0f, 1.0f));
			dir = dir.normalized();
			for (int i = 0; i < 10000; i++)
			{
				points << Vec3(100, 50, 20) + dir * random(-500.f, 500.f) +
				              Vec3(asl::random(-k, k), asl::random(-k, k), asl::random(-k, k));
			}

			auto line = fitLine(points);

			//printf("%f %f %f (dir: %f %f %f)\n", line.first.x, line.first.y, line.first.z, line.second.x, line.second.y,
			  //     line.second.z);

			// compute distance of points to that line
			float dmax = 0.0f;
			for (int i = 0; i < points.length(); i++)
			{
				float d = distancePointLine(points[i], line.first, line.second);
				if (d > dmax)
					dmax = d;
			}

			printf("max distance to line: %f\n", dmax);

			auto plane = fitPlaneXYR(points, 2.0f);

			auto line2 = fitLineR(points, 2.0f);

			auto poly = fitPolyR(points, 3, 2.0f);

			auto sphere = fitSphereR(points, 2.0f);
		}

		return 0;
	}
#endif

	if (0)
	{
		Array<Vec3d> points;
		Array<Vec3d> normals;

		String   dir = "I:\\work\\sariki\\data3d\\";
		TextFile file(dir + "n06.xyz");
		while (!file.end())
		{
			Vec3d p;
			Vec3d n;
			file >> p.x >> p.y >> p.z >> n.x >> n.y >> n.z;
			points << p;
			normals << n;
		}

		Array<Vec4d> dataX, dataY, dataZ;
		for (int i = 0; i < points.length(); i++)
		{
			dataX << Vec4d(points[i], normals[i].x);
			dataY << Vec4d(points[i], normals[i].y);
			dataZ << Vec4d(points[i], normals[i].z);
		}

		auto polyX = fitPoly(dataX, 3);
		auto polyY = fitPoly(dataY, 3);
		auto polyZ = fitPoly(dataZ, 3);

		// print polynomials
		printf("polyX: %s\n", *polyX.array().join(", "));
		printf("polyY: %s\n", *polyY.array().join(", "));
		printf("polyZ: %s\n", *polyZ.array().join(", "));

		TextFile file2(dir + "n06-3.xyz");
		for (int i = 0; i < points.length(); i++)
		{
			auto  nx = polynomial(polyX, points[i].x, points[i].y, points[i].z);
			auto  ny = polynomial(polyY, points[i].x, points[i].y, points[i].z);
			auto  nz = polynomial(polyZ, points[i].x, points[i].y, points[i].z);
			Vec3d n = Vec3d(nx, ny, nz).normalized();
			file2 << points[i].x << " " << points[i].y << " " << points[i].z << " " << n.x << " " << n.y << " " << n.z
			      << "\n";
		}

		// return 0;
	}
	asl::TextFile file("normals.xyz", asl::TextFile::WRITE);
	asl::random.seed(0);
	int good = 0, bad = 0, verybad = 0;
	for (int j = 0; j < 1000; j++)
	{
		Array<Vec3> points;
		float       rx = asl::random(-6.3f, 6.3f);
		float       ry = asl::random(-6.3f, 6.3f);
		float       rz = asl::random(-3.3f, 3.3f);
		Vec3        center(asl::random(-400.0f, 400.0f), asl::random(-400.0f, 400.0f), asl::random(-100.0f, 100.0f));
		// print rx,ry
		printf("rx = %f ry = %f rz = %f\n", rx, ry, rz);
		for (int i = 0; i < 2000; i++)
		{
			points << center + Matrix4::rotateX(rx) * Matrix4::rotateY(ry) *
			             Matrix4::rotateZ(rz) *
			              Vec3(asl::random(150.0f), asl::random(8.0f), asl::random(0.1f));
		}

		double t01 = now();

		//auto plane = fitPlaneXYZ(points.with<Vec3d>()).with<float>();
		auto plane = fitPlaneXYZ(points);

		double t02 = now();

		printf("%s (t plane = %f)\n", *plane.join(", "), t02 - t01);

		t01 = now();

		//auto plane3D = fitPlane(points.with<Vec3d>()).with<float>();
		auto plane3D = fitPlane(points);

		t02 = now();

		printf("%s (t plane = %f)\n", *plane3D.join(", "), t02 - t01);

		auto normal1 = Vec3(plane[3], plane[4], plane[5]);
		auto normal2 = Vec3(plane3D[3], plane3D[4], plane3D[5]);

		auto planepoint1 = Vec3(plane[0], plane[1], plane[2]);
		auto planepoint2 = Vec3(plane3D[0], plane3D[1], plane3D[2]);

		if (normal1 * normal2 < 0)
			normal1 = -normal1;

		file.printf("%f %f %f\n", normal2.x, normal2.y, normal2.z);

		printf("angle mismatch = %.9f\n", rad2deg(acos(clamp(normal1 * normal2, -1.0f, 1.0f))));
		// printf("center mismatch = %.3f\n", (planepoint1 - planepoint2).length());

		// compute max distance of points to that plane
		// and sum of squared distances
		float dmax1 = 0.0f;
		float dsum1 = 0.0f;
		for (int i = 0; i < points.length(); i++)
		{
			float d = (points[i] - planepoint1) * normal1;
			if (d > dmax1)
				dmax1 = d;
			dsum1 += d * d;
		}

		printf("distance to plane1: max %5.3f   sumsq %10.3f\n", dmax1, dsum1);
		float dmax2 = 0;
		float dsum2 = 0;
		for (int i = 0; i < points.length(); i++)
		{
			float d = (points[i] - planepoint2) * normal2;
			if (d > dmax2)
				dmax2 = d;
			dsum2 += d * d;
		}

		printf("distance to plane2: max %5.3f   sumsq %10.3f\n", dmax2, dsum2);

		if (dsum2 > dsum1 * 1.01)
		{
			verybad++;

			printf("*********** fitPlane is much worse than fitPlaneXYZ ***********\n");
		}
		else if (dsum2 > dsum1 * 1.00001)
		{
			bad++;
			printf("******* fitPlane is worse than fitPlaneXYZ *******\n");
		}
		else
			good++;
	}

	printf("good: %d bad: %d verybad: %d / total: %d\n", good, bad, verybad, good + bad + verybad);

	return 0;

	/*
	Array<Vec2> points1 = points.map_<Vec2>([](Vec3 p) { return p.xy(); });

	auto tr = Matrix3::translate(50, 30) * Matrix3::rotate(0.5f);

	Array<Vec2> points2 = points1.map([=](Vec2 p) { return (tr * p) + Vec2(asl::random(-k, k), asl::random(-k, k)); });

	double t1 = now();

	auto tr2 = findRigidTransform(points1, points2);

	double t2 = now();

	printf("t solve rigid 2D: %f s (%f)\n", t2 - t1, (tr2.inverse() * tr).trace() - 3);

	// fill points with a grid of points in X and Y
	points.clear();
	*/
	Array<Vec3> points;
	for (int i = 0; i < 200; i++)
		for (int j = 0; j < 200; j++)
		{
			points << Vec3{ j * 10.0f, i * 10.0f, 0 } +
			              Vec3{ asl::random(-1.0f, 1.0f), asl::random(-1.0f, 1.0f), asl::random(-1.0f, 1.0f) };
		}

	// repeat this block 20 times with different translation and rotation

	float posErrorMax = 0.0f;
	float posErrorSum = 0.0f;
	float rotErrorMax = 0.0f;
	float rotErrorSum = 0.0f;

	for (int iter = 0; iter < 20; iter++)
	{
		Vec3 rot0 = { asl::random(-1.0f, 1.0f), asl::random(-1.0f, 1.0f), asl::random(-1.0f, 1.0f) };
		Vec3 pos0 = { asl::random(-100.0f, 100.0f), asl::random(-100.0f, 100.0f), asl::random(-50.0f, 50.0f) };

		auto trb = Matrix4::translate(pos0) * Matrix4::rotateE(rot0, "XYZ*");

		// auto trb = Matrix4::translate(trans0) * Matrix4::rotateE({ 0.7f, 0.0f, 0.5f }, "XYZ*");

		Array<Vec3> points2b = points.map(
		    [=](Vec3 p) { return (trb * p) + 1.0f * Vec3(asl::random(-k, k), asl::random(-k, k), asl::random(-k, k)); });

		t1 = now();

		auto tr2b = findRigidTransform(points2b, points, 3).inverse();

		t2 = now();

		// printf("t solve rigid 3D: %f s (%f)\n", t2 - t1, (tr2b.inverse() * trb).norm() - 2);

		auto rot = tr2b.eulerAngles("XYZ*");
		auto pos = tr2b.translation();
		printf("pos0 %f %f %f rot0 = %f %f %f\n", pos0.x, pos0.y, pos0.z, rot0.x, rot0.y, rot0.z);
		printf("pos  %f %f %f rot  = %f %f %f P err %g R err %g\n", pos.x, pos.y, pos.z, rot.x, rot.y, rot.z,
		       (pos - pos0).length(), (rot - rot0).length());

		float pe = (pos - pos0).length();
		float re = (rot - rot0).length();
		posErrorSum += pe;
		rotErrorSum += re;
		if (pe > posErrorMax)
			posErrorMax = pe;
		if (re > rotErrorMax)
			rotErrorMax = re;
	}

	printf("pos err avg %g max %g rot err avg %g max %g\n", posErrorSum / 20.0f, posErrorMax, rotErrorSum / 20.0f,
	       rotErrorMax);

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
