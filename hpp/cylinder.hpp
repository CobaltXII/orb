#pragma once

struct cylinder: shape
{
	float a_x;
	float a_y;
	float a_z;

	float b_x;
	float b_y;
	float b_z;

	float radius;

	cylinder
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
		this->primitive = shape_type::st_cylinder;

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

inline cylinder TO_CYLINDER(shape* __victim)
{
	return *((cylinder*)__victim);
}

inline float cylinder_intersect
(
	cylinder cylinder1,

	float ray_ox, float ray_oy, float ray_oz,
	float ray_dx, float ray_dy, float ray_dz,

	float* norm_x,
	float* norm_y,
	float* norm_z
)
{
	float ca_x = cylinder1.b_x - cylinder1.a_x;
	float ca_y = cylinder1.b_y - cylinder1.a_y;
	float ca_z = cylinder1.b_z - cylinder1.a_z;

	float oc_x = ray_ox - cylinder1.a_x;
	float oc_y = ray_oy - cylinder1.a_y;
	float oc_z = ray_oz - cylinder1.a_z;

	float caca =
	(
		ca_x * ca_x +
		ca_y * ca_y +
		ca_z * ca_z
	);

	float card =
	(
		ca_x * ray_dx +
		ca_y * ray_dy +
		ca_z * ray_dz
	);

	float caoc =
	(
		ca_x * oc_x +
		ca_y * oc_y +
		ca_z * oc_z
	);

	float a = caca - card * card;

	float ocrd =
	(
		oc_x * ray_dx +
		oc_y * ray_dy +
		oc_z * ray_dz
	);

	float b = caca * ocrd - caoc * card;

	float ococ =
	(
		oc_x * oc_x +
		oc_y * oc_y +
		oc_z * oc_z
	);

	float c =
	(
		caca * ococ -
		caoc * caoc -

		cylinder1.radius *
		cylinder1.radius

		* caca
	);

	float h = b * b - a * c;

	if (h < 0.0f)
	{
		return -1.0f;
	}

	h = sqrtf(h);

	float t = (-b - h) / a;

	float y = caoc + t * card;

	if (y > 0.0f && y < caca)
	{
		set_ptr(norm_x, (oc_x + t * ray_dx - ca_x * y / caca) / cylinder1.radius);
		set_ptr(norm_y, (oc_y + t * ray_dy - ca_y * y / caca) / cylinder1.radius);
		set_ptr(norm_z, (oc_z + t * ray_dz - ca_z * y / caca) / cylinder1.radius);

		return t;
	}

	if (y < 0.0f)
	{
		t = -caoc / card;
	}
	else
	{
		t = (caca - caoc) / card;
	}

	if (fabsf(b + a * t) < h)
	{
		set_ptr(norm_x, ca_x * sign(y) / caca);
		set_ptr(norm_y, ca_y * sign(y) / caca);
		set_ptr(norm_z, ca_z * sign(y) / caca);

		return t;
	}

	return -1.0f;
}