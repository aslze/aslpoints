#pragma once
#include <asl/Array.h>
#include <asl/Vec2.h>

// 2D Signed distance functions

namespace asl
{

/**
 * SDF of a circle of radius r centred at (0,0)
 */
template<class T>
T sdCircle(const Vec2_<T>& p, T r)
{
	return p.length() - r;
}

/**
 * SDF of a line segment from a to b
 */
template<class T>
T sdSegment(const Vec2_<T>& p, const Vec2_<T>& a, const Vec2_<T>& b)
{
	auto pa = p - a, ba = b - a;
	T    h = clamp((pa * ba) / (ba * ba), T(0), T(1));
	return (pa - ba * h).length();
}

/**
 * SDF of an arc of radius ra, with optional thickness rb, centred at (0,0) spanning angle defined by sc = (sin,cos)
 */
template<class T>
T sdArc(Vec2_<T> p, const Vec2_<T>& sc, T ra, T rb = T(0))
{
	p.x = fabs(p.x);
	return ((sc.y * p.x > sc.x * p.y) ? (p - sc * ra).length() : fabs(p.length() - ra)) - rb;
}

/**
 * SDF of a capsule aligned with Y axis and centred at (0, 0)
 */
template<class T>
T sdCapsuleY(Vec2_<T> p, T r, T h)
{
	p.x = abs(p.x);
	p.y += h / 2; // symmetric around X

	return (p.y < 0) ? (sqrt(p * p) - r) : (p.y > h) ? (sqrt(p * p + h * h - T(2) * h * p.y) - r) : (p.x - r);
}

/**
 * SDF of a rectangle half size bx, by
 */
template<class T>
inline T sdBox(const Vec2_<T>& p, T bx, T by)
{
	auto d = Vec2_<T>(fabs(p.x) - bx, fabs(p.y) - by);
	return max(d, Vec2_<T>(0, 0)).length() + min(max(d.x, d.y), T(0));
}

/**
 * SDF of a rounded rectangle of half size bx, by and rounding r
 */
template<class T>
inline T sdBoxR(const Vec2_<T>& p, T bx, T by, T r)
{
	return sdBox(p, bx - r, by - r) - r;
}

template<class T>
T sdTrapezoid(Vec2_<T> p, T r1, T r2, T he)
{
	Vec2_<T> k1(r2, he);
	Vec2_<T> k2(r2 - r1, T(2.0) * he);
	p.x = abs(p.x);
	Vec2_<T> ca(p.x - min(p.x, (p.y < 0.0) ? r1 : r2), abs(p.y) - he);
	Vec2_<T> cb = p - k1 + k2 * clamp(((k1 - p) * k2) / k2.length2(), T(0.0), T(1.0));
	T        s = (cb.x < 0.0 && ca.y < 0.0) ? T(-1.0) : T(1.0);
	return s * sqrt(min(ca.length2(), cb.length2()));
}

/**
 * SDF of a triangle defined by points p0, p1, p2
 */
template<class T>
T sdTriangle(const Vec2_<T>& p, const Vec2_<T>& p0, const Vec2_<T>& p1, const Vec2_<T>& p2)
{
	auto e0 = p1 - p0, e1 = p2 - p1, e2 = p0 - p2;
	auto v0 = p - p0, v1 = p - p1, v2 = p - p2;
	auto pq0 = v0 - e0 * clamp((v0 * e0) / (e0 * e0), T(0), T(1));
	auto pq1 = v1 - e1 * clamp((v1 * e1) / (e1 * e1), T(0), T(1));
	auto pq2 = v2 - e2 * clamp((v2 * e2) / (e2 * e2), T(0), T(1));
	T    s = (e0.x * e2.y - e0.y * e2.x) > 0 ? T(1) : T(-1);
	auto d = min(min(Vec2_<T>((pq0 * pq0), s * (v0.x * e0.y - v0.y * e0.x)),
	                 Vec2_<T>((pq1 * pq1), s * (v1.x * e1.y - v1.y * e1.x))),
	             Vec2_<T>((pq2 * pq2), s * (v2.x * e2.y - v2.y * e2.x)));
	return -sqrt(d.x) * (d.y > 0 ? T(1) : T(-1));
}

/**
 * SDF of a polygon defined by vertices v
 */
template<class T>
T sdPolygon(const Vec2_<T>& p, const Array<Vec2_<T>>& v)
{
	T   d = (p - v[0]) * (p - v[0]);
	T   s = T(1);
	int N = v.length();
	for (int i = 0, j = N - 1; i < N; j = i, i++)
	{
		auto e = v[j] - v[i];
		auto w = p - v[i];
		auto b = w - e * clamp((w * e) / (e * e), T(0), T(1));
		d = min(d, b * b);
		bool c1 = p.y >= v[i].y, c2 = p.y<v[j].y, c3 = e.x * w.y> e.y * w.x;
		if (c1 && c2 && c3 || !c1 && !c2 && !c3)
			s *= T(-1);
	}
	return s * sqrt(d);
}

/**
 * Union operation between two SDFs
 */
template<class T>
T opUnion(T d1, T d2)
{
	return min(d1, d2);
}

/**
 * Subtraction operation between two SDFs
 */
template<class T>
T opSubtraction(T d1, T d2)
{
	return max(-d1, d2);
}

/**
 * Smooth union operation between two SDFs with smoothing factor k
 */
template<class T>
T opSmoothUnion(T d1, T d2, T k)
{
	k *= T(4.0);
	T h = max(k - abs(d1 - d2), T(0));
	return min(d1, d2) - h * h * T(0.25) / k;
}

/**
 * Smooth subtraction operation between two SDFs with smoothing factor k
 */
template<class T>
T opSmoothSubtraction(T d1, T d2, T k)
{
	k *= T(4.0);
	T h = max(k - abs(-d1 - d2), T(0));
	return max(-d1, d2) + h * h * T(0.25) / k;
}

}
