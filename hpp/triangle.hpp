#pragma once

struct triangle: shape
{
	float x0;
	float y0;
	float z0;

	float x1;
	float y1;
	float z1;

	float x2;
	float y2;
	float z2;

	float norm_x;
	float norm_y;
	float norm_z;

	triangle
	(
		float x0,
		float y0,
		float z0,

		float x1,
		float y1,
		float z1,

		float x2,
		float y2,
		float z2,

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
		this->primitive = shape_type::st_triangle;

		this->x0 = x0;
		this->y0 = y0;
		this->z0 = z0;

		this->x1 = x1;
		this->y1 = y1;
		this->z1 = z1;

		this->x2 = x2;
		this->y2 = y2;
		this->z2 = z2;

		this->r = r;
		this->g = g;
		this->b = b;

		this->material1 = material1;
		this->material2 = material2;
		this->material3 = material3;
		this->material4 = material4;
		this->material5 = material5;

		// Surface normal.

		float v1v0_x = x1 - x0;
		float v1v0_y = y1 - y0;
		float v1v0_z = z1 - z0;

		float v2v0_x = x2 - x0;
		float v2v0_y = y2 - y0;
		float v2v0_z = z2 - z0;

		float cross_x = v1v0_y * v2v0_z - v2v0_y * v1v0_z;
		float cross_y = v1v0_x * v2v0_z - v2v0_x * v1v0_z;
		float cross_z = v1v0_x * v2v0_y - v2v0_x * v1v0_y;

		float cross_len = sqrtf
		(
			cross_x * cross_x +
			cross_y * cross_y +
			cross_z * cross_z
		);

		norm_x = cross_x / cross_len;
		norm_y = cross_y / cross_len;
		norm_z = cross_z / cross_len;
	}
};

inline triangle TO_TRIANGLE(shape* __victim)
{
	return *((triangle*)__victim);
}

inline float triangle_intersect
(
	triangle triangle1,

	float ray_ox, float ray_oy, float ray_oz,
	float ray_dx, float ray_dy, float ray_dz,

	float* norm_x,
	float* norm_y,
	float* norm_z,

	float* texture_u,
	float* texture_v
)
{
	#define v0_x (triangle1.x0)
	#define v0_y (triangle1.y0)
	#define v0_z (triangle1.z0)

	#define v1_x (triangle1.x1)
	#define v1_y (triangle1.y1)
	#define v1_z (triangle1.z1)

	#define v2_x (triangle1.x2)
	#define v2_y (triangle1.y2)
	#define v2_z (triangle1.z2)

	float v0v1_x = v1_x - v0_x;
	float v0v1_y = v1_y - v0_y;
	float v0v1_z = v1_z - v0_z;

	float v0v2_x = v2_x - v0_x;
	float v0v2_y = v2_y - v0_y;
	float v0v2_z = v2_z - v0_z;

	float pvec_x = ray_dy * v0v2_z - ray_dz * v0v2_y;
	float pvec_y = ray_dz * v0v2_x - ray_dx * v0v2_z;
	float pvec_z = ray_dx * v0v2_y - ray_dy * v0v2_x;

	float inv_det = 1.0f /
	(
		v0v1_x * pvec_x +
		v0v1_y * pvec_y +
		v0v1_z * pvec_z
	);

	float tvec_x = ray_ox - v0_x;
	float tvec_y = ray_oy - v0_y;
	float tvec_z = ray_oz - v0_z;

	float u = inv_det *
	(
		tvec_x * pvec_x +
		tvec_y * pvec_y +
		tvec_z * pvec_z
	);

	if (u < 0.0f || u > 1.0f)
	{
		return -1.0f;
	}

	float qvec_x = tvec_y * v0v1_z - tvec_z * v0v1_y;
	float qvec_y = tvec_z * v0v1_x - tvec_x * v0v1_z;
	float qvec_z = tvec_x * v0v1_y - tvec_y * v0v1_x;

	float v = inv_det *
	(
		ray_dx * qvec_x +
		ray_dy * qvec_y +
		ray_dz * qvec_z
	);

	if (v < 0.0f || u + v > 1.0f)
	{
		return -1.0f;
	}

	float t = inv_det *
	(
		v0v2_x * qvec_x +
		v0v2_y * qvec_y +
		v0v2_z * qvec_z
	);

	set_ptr(norm_x, triangle1.norm_x);
	set_ptr(norm_y, triangle1.norm_y);
	set_ptr(norm_z, triangle1.norm_z);

	set_ptr(texture_u, u);
	set_ptr(texture_v, v);

	return t;
}