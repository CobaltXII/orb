#pragma once

struct capsule: shape
{
	float a_x;
	float a_y;
	float a_z;

	float b_x;
	float b_y;
	float b_z;

	float radius;

	capsule
	(
		float a_x,
		float a_y,
		float a_z,

		float b_x,
		float b_y,
		float b_z,

		float r,
		float g,
		float b,

		float radius,

		float material1 = 0.0f,
		float material2 = 0.0f,
		float material3 = 0.0f,
		float material4 = 0.0f,
		float material5 = 0.0f
	)
	{
		this->primitive = shape_type::st_capsule;

		this->a_x = a_x;
		this->a_y = a_y;
		this->a_z = a_z;

		this->b_x = b_x;
		this->b_y = b_y;
		this->b_z = b_z;

		this->r = r;
		this->g = g;
		this->b = b;

		this->radius = radius;

		this->material1 = material1;
		this->material2 = material2;
		this->material3 = material3;
		this->material4 = material4;
		this->material5 = material5;
	}
};

inline capsule TO_CAPSULE(shape* __victim)
{
	return *((capsule*)__victim);
}

inline float capsule_intersect
(
	capsule capsule1,

	float ray_ox, float ray_oy, float ray_oz,
	float ray_dx, float ray_dy, float ray_dz,

	float* norm_x,
	float* norm_y,
	float* norm_z
)
{
	float ba_x = capsule1.b_x - capsule1.a_x;
	float ba_y = capsule1.b_y - capsule1.a_y;
	float ba_z = capsule1.b_z - capsule1.a_z;

	float oa_x = ray_ox - capsule1.a_x;
	float oa_y = ray_oy - capsule1.a_y;
	float oa_z = ray_oz - capsule1.a_z;

	float baba =
	(
		ba_x * ba_x +
		ba_y * ba_y +
		ba_z * ba_z
	);

	float bard =
	(
		ba_x * ray_dx +
		ba_y * ray_dy +
		ba_z * ray_dz
	);

	float baoa =
	(
		ba_x * oa_x +
		ba_y * oa_y +
		ba_z * oa_z
	);

	float rdoa =
	(
		ray_dx * oa_x +
		ray_dy * oa_y +
		ray_dz * oa_z
	);

	float oaoa =
	(
		oa_x * oa_x +
		oa_y * oa_y +
		oa_z * oa_z
	);

	float a = baba - bard * bard;

	float b = baba * rdoa - baoa * bard;

	float c =
	(
		baba * oaoa -
		baoa * baoa -

		capsule1.radius *
		capsule1.radius

		* baba
	);

	float h = b * b - a * c;

	if (h >= 0.0f)
	{
		float t = (-b - sqrtf(h)) / a;

		float y = baoa + t * bard;

		if (y > 0.0f && y < baba)
		{
			float hit_x = ray_ox + ray_dx * t;
			float hit_y = ray_oy + ray_dy * t;
			float hit_z = ray_oz + ray_dz * t;

			float pa_x = hit_x - capsule1.a_x;
			float pa_y = hit_y - capsule1.a_y;
			float pa_z = hit_z - capsule1.a_z;

			float paba =
			(
				pa_x * ba_x +
				pa_y * ba_y +
				pa_z * ba_z
			);

			float q = fmax(0.0f, fmin(1.0f, paba / baba));

			set_ptr(norm_x, (pa_x - q * ba_x) / capsule1.radius);
			set_ptr(norm_y, (pa_y - q * ba_y) / capsule1.radius);
			set_ptr(norm_z, (pa_z - q * ba_z) / capsule1.radius);

			return t;
		}

		float oc_x;
		float oc_y;
		float oc_z;

		if (y <= 0.0f)
		{
			oc_x = oa_x;
			oc_y = oa_y;
			oc_z = oa_z;
		}
		else
		{
			oc_x = ray_ox - capsule1.b_x;
			oc_y = ray_oy - capsule1.b_y;
			oc_z = ray_oz - capsule1.b_z;
		}

		b =
		(
			ray_dx * oc_x +
			ray_dy * oc_y +
			ray_dz * oc_z
		);

		c =
		(
			oc_x * oc_x +
			oc_y * oc_y +
			oc_z * oc_z

			- capsule1.radius * capsule1.radius
		);

		h = b * b - c;

		if (h > 0.0f)
		{
			t = -b - sqrtf(h);

			float hit_x = ray_ox + ray_dx * t;
			float hit_y = ray_oy + ray_dy * t;
			float hit_z = ray_oz + ray_dz * t;

			float pa_x = hit_x - capsule1.a_x;
			float pa_y = hit_y - capsule1.a_y;
			float pa_z = hit_z - capsule1.a_z;

			float paba =
			(
				pa_x * ba_x +
				pa_y * ba_y +
				pa_z * ba_z
			);

			float q = fmax(0.0f, fmin(1.0f, paba / baba));

			set_ptr(norm_x, (pa_x - q * ba_x) / capsule1.radius);
			set_ptr(norm_y, (pa_y - q * ba_y) / capsule1.radius);
			set_ptr(norm_z, (pa_z - q * ba_z) / capsule1.radius);

			return t;
		}
	}

	return -1.0f;
}