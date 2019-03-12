#pragma once

enum shape_type
{
	st_sphere, st_plane, st_ellipsoid, st_cone, st_capsule, st_cylinder, st_triangle, st_intersection
};

struct shape
{
	float material1;
	float material2;
	float material3;
	float material4;
	float material5;

	float r;
	float g;
	float b;

	shape_type primitive;
};