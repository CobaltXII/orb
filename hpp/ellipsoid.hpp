#pragma once

struct ellipsoid: shape
{
	float x;
	float y;
	float z;

	float radius_x;
	float radius_y;
	float radius_z;

	ellipsoid
	(
		float x,
		float y,
		float z,

		float r,
		float g,
		float b,

		float radius_x,
		float radius_y,
		float radius_z,

		float material1 = 0.0f,
		float material2 = 0.0f,
		float material3 = 0.0f,
		float material4 = 0.0f,
		float material5 = 0.0f
	)
	{
		this->primitive = shape_type::st_ellipsoid;

		this->x = x;
		this->y = y;
		this->z = z;

		this->r = r;
		this->g = g;
		this->b = b;

		this->radius_x = radius_x;
		this->radius_y = radius_y;
		this->radius_z = radius_z;

		this->material1 = material1;
		this->material2 = material2;
		this->material3 = material3;
		this->material4 = material4;
		this->material5 = material5;
	}
};

inline ellipsoid TO_ELLIPSOID(shape* __victim)
{
	return *((ellipsoid*)__victim);
}

inline float ellipsoid_intersect
(
	ellipsoid ellipsoid1,

	float ray_ox, float ray_oy, float ray_oz,
	float ray_dx, float ray_dy, float ray_dz,

	float* norm_x,
	float* norm_y,
	float* norm_z
)
{
	float oc_x = ray_ox - ellipsoid1.x;
	float oc_y = ray_oy - ellipsoid1.y;
	float oc_z = ray_oz - ellipsoid1.z;

	float a =
	(
		(ray_dx * ray_dx) / (ellipsoid1.radius_x * ellipsoid1.radius_x) +
		(ray_dy * ray_dy) / (ellipsoid1.radius_y * ellipsoid1.radius_y) +
		(ray_dz * ray_dz) / (ellipsoid1.radius_z * ellipsoid1.radius_z)
	);

	float b =
	(
		(2.0f * oc_x * ray_dx) / (ellipsoid1.radius_x * ellipsoid1.radius_x) +
		(2.0f * oc_y * ray_dy) / (ellipsoid1.radius_y * ellipsoid1.radius_y) +
		(2.0f * oc_z * ray_dz) / (ellipsoid1.radius_z * ellipsoid1.radius_z)
	);

	float c =
	(
		(oc_x * oc_x) / (ellipsoid1.radius_x * ellipsoid1.radius_x) +
		(oc_y * oc_y) / (ellipsoid1.radius_y * ellipsoid1.radius_y) +
		(oc_z * oc_z) / (ellipsoid1.radius_z * ellipsoid1.radius_z)

		- 1.0f
	);

	float d = b * b - 4.0f * a * c;

	if (d < 0.0f || a == 0.0f || b == 0.0f || c == 0.0f)
	{
		return -1.0f;
	}

	d = sqrtf(d);

	float t1 = (-b + d) / (2.0f * a);
	float t2 = (-b - d) / (2.0f * a);

	float eps = 0.0f;

	if (t1 <= eps && t2 <= eps)
	{
		return -1.0f;
	}

	float t = 0.0f;

	if (t1 <= eps)
	{
		t = t2;
	}
	else
	{
		if (t2 <= eps)
		{
			t = t1;
		}
		else
		{
			if (t1 < t2)
			{
				t = t1;
			}
			else
			{
				t = t2;
			}
		}
	}

	if (t < eps)
	{
		return -1.0f;
	}

	set_ptr(norm_x, ((ray_ox + ray_dx * t) - ellipsoid1.x) / ellipsoid1.radius_x);
	set_ptr(norm_y, ((ray_oy + ray_dy * t) - ellipsoid1.y) / ellipsoid1.radius_y);
	set_ptr(norm_z, ((ray_oz + ray_dz * t) - ellipsoid1.z) / ellipsoid1.radius_z);

	return t;
}