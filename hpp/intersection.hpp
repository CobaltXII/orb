#pragma once

struct intersection: shape
{
	shape* sa;
	shape* sb;

	intersection
	(
		shape* sa,
		shape* sb,

		float r,
		float g,
		float b,

		float material1 = 0.0f,
		float material2 = 0.0f,
		float material3 = 0.0f,
		float material4 = 0.0f,
		float material5 = 0.0f
	)
	{
		this->primitive = shape_type::st_intersection;

		this->sa = sa;
		this->sb = sb;

		this->r = r;
		this->g = g;
		this->b = b;

		this->material1 = material1;
		this->material2 = material2;
		this->material3 = material3;
		this->material4 = material4;
		this->material5 = material5;
	}
};

inline intersection TO_INTERSECTION(shape* __victim)
{
	return *((intersection*)__victim);
}