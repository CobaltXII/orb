#pragma once

struct cone: shape
{
	float a_x;
	float a_y;
	float a_z;

	float b_x;
	float b_y;
	float b_z;

	float radius_a;
	float radius_b;

	cone
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

		float radius_a,
		float radius_b,

		float material1 = 0.0f,
		float material2 = 0.0f,
		float material3 = 0.0f,
		float material4 = 0.0f,
		float material5 = 0.0f
	)
	{
		this->primitive = shape_type::st_cone;

		this->a_x = a_x;
		this->a_y = a_y;
		this->a_z = a_z;

		this->b_x = b_x;
		this->b_y = b_y;
		this->b_z = b_z;

		this->r = r;
		this->g = g;
		this->b = b;

		this->radius_a = radius_a;
		this->radius_b = radius_b;

		this->material1 = material1;
		this->material2 = material2;
		this->material3 = material3;
		this->material4 = material4;
		this->material5 = material5;
	}
};

inline cone TO_CONE(shape* __victim)
{
	return *((cone*)__victim);
}

inline float cone_intersect
(
	cone cone1,

	float ray_ox, float ray_oy, float ray_oz,
	float ray_dx, float ray_dy, float ray_dz,

	float* norm_x,
	float* norm_y,
	float* norm_z
)
{
	float ba_x = cone1.b_x - cone1.a_x;
	float ba_y = cone1.b_y - cone1.a_y;
	float ba_z = cone1.b_z - cone1.a_z;

	float oa_x = ray_ox - cone1.a_x;
	float oa_y = ray_oy - cone1.a_y;
	float oa_z = ray_oz - cone1.a_z;

	float ob_x = ray_ox - cone1.b_x;
	float ob_y = ray_oy - cone1.b_y;
	float ob_z = ray_oz - cone1.b_z;

	float baba =
	(
		ba_x * ba_x +
		ba_y * ba_y +
		ba_z * ba_z
	);

	float rdba =
	(
		ray_dx * ba_x +
		ray_dy * ba_y +
		ray_dz * ba_z
	);

	float oaba =
	(
		oa_x * ba_x +
		oa_y * ba_y +
		oa_z * ba_z
	);

	float obba =
	(
		ob_x * ba_x +
		ob_y * ba_y +
		ob_z * ba_z
	);

	if (oaba < 0.0f)
	{
		float c1_x = oa_x * rdba - ray_dx * oaba;
		float c1_y = oa_y * rdba - ray_dy * oaba;
		float c1_z = oa_z * rdba - ray_dz * oaba;

		float c1 =
		(
			c1_x * c1_x +
			c1_y * c1_y +
			c1_z * c1_z
		);

		if (c1 < cone1.radius_a * cone1.radius_a * rdba * rdba)
		{
			float isqr_baba = 1.0f / sqrtf(baba);

			set_ptr(norm_x, -ba_x * isqr_baba);
			set_ptr(norm_y, -ba_y * isqr_baba);
			set_ptr(norm_z, -ba_z * isqr_baba);

			return -oaba / rdba;
		}
	}
	else if (obba > 0.0f)
	{
		float t = -obba / rdba;

		float c1_x = ob_x + ray_dx * t;
		float c1_y = ob_y + ray_dy * t;
		float c1_z = ob_z + ray_dz * t;

		float c1 =
		(
			c1_x * c1_x +
			c1_y * c1_y +
			c1_z * c1_z
		);

		if (c1 < cone1.radius_b * cone1.radius_b)
		{
			float isqr_baba = 1.0f / sqrtf(baba);

			set_ptr(norm_x, ba_x * isqr_baba);
			set_ptr(norm_y, ba_y * isqr_baba);
			set_ptr(norm_z, ba_z * isqr_baba);

			return t;
		}
	}

	float rbra = cone1.radius_b - cone1.radius_a;

	float oc_x = oa_x * cone1.radius_b - ob_x * cone1.radius_a;
	float oc_y = oa_y * cone1.radius_b - ob_y * cone1.radius_a;
	float oc_z = oa_z * cone1.radius_b - ob_z * cone1.radius_a;

	float hyhy = baba + rbra * rbra;

	float ocba =
	(
		oc_x * ba_x +
		oc_y * ba_y +
		oc_z * ba_z
	);

	float ocrd =
	(
		oc_x * ray_dx +
		oc_y * ray_dy +
		oc_z * ray_dz
	);

	float ococ =
	(
		oc_x * oc_x +
		oc_y * oc_y +
		oc_z * oc_z
	);

	float k2 = baba * baba - hyhy * rdba * rdba;

	float k1 = baba * baba * ocrd - hyhy * rdba * ocba;
	float k0 = baba * baba * ococ - hyhy * ocba * ocba;

	float h = k1 * k1 - k2 * k0;

	if (h < 0.0f)
	{
		return -1.0f;
	}

	float t = (-k1 - sign(rbra) * sqrtf(h)) / (k2 * rbra);

	if (t < 0.0f)
	{
		return -1.0f;
	}

	float y = oaba + t * rdba;

	if (y > 0.0f && y < baba)
	{
		set_ptr(norm_x, baba * (baba * (oa_x + t * ray_dx) - rbra * ba_x * cone1.radius_a) - ba_x * hyhy * y);
		set_ptr(norm_y, baba * (baba * (oa_y + t * ray_dy) - rbra * ba_y * cone1.radius_a) - ba_y * hyhy * y);
		set_ptr(norm_z, baba * (baba * (oa_z + t * ray_dz) - rbra * ba_z * cone1.radius_a) - ba_z * hyhy * y);

		return t;
	}

	return -1.0f;
}