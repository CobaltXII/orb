#pragma once

struct light
{
	float x;
	float y;
	float z;

	float r;
	float g;
	float b;

	float radius;

	light
	(
		float x,
		float y,
		float z,

		float r,
		float g,
		float b,

		float radius
	)
	{
		this->x = x;
		this->y = y;
		this->z = z;

		this->r = r;
		this->g = g;
		this->b = b;

		this->radius = radius;
	}
};