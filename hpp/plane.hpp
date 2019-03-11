#pragma once

struct plane: shape
{
	float x;
	float y;
	float z;

	float norm_x;
	float norm_y;
	float norm_z;

	plane
	(
		float x,
		float y,
		float z,

		float norm_x,
		float norm_y,
		float norm_z,

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
		this->primitive = shape_type::st_plane;

		this->x = x;
		this->y = y;
		this->z = z;

		this->norm_x = norm_x;
		this->norm_y = norm_y;
		this->norm_z = norm_z;

		float norm_length = sqrtf
		(
			this->norm_x * this->norm_x +
			this->norm_y * this->norm_y +
			this->norm_z * this->norm_z
		);

		this->norm_x /= norm_length;
		this->norm_y /= norm_length;
		this->norm_z /= norm_length;

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

inline plane TO_PLANE(shape* __victim)
{
	return *((plane*)__victim);
}

inline void plane_uv
(
	float intersection_x,
	float intersection_y,
	float intersection_z,

	float plane_origin_x,
	float plane_origin_y,
	float plane_origin_z,

	float plane_normal_x,
	float plane_normal_y,
	float plane_normal_z,

	float& out_u,
	float& out_v
)
{
	if (plane_normal_y == 0.0f)
	{
		plane_normal_y = plane_normal_z;
	}

	float x_axis_x = plane_normal_y * 1.0f - 0.0f * plane_normal_z;
	float x_axis_y = plane_normal_x * 1.0f - 0.0f * plane_normal_z;
	float x_axis_z = plane_normal_x * 0.0f - 0.0f * plane_normal_y;

	float x_axis_length =
	(
		x_axis_x * x_axis_x +
		x_axis_y * x_axis_y +
		x_axis_z * x_axis_z
	);

	if (x_axis_length <= 0.0f)
	{
		x_axis_x = plane_normal_y * 0.0f - 0.0f * plane_normal_y;
		x_axis_y = plane_normal_z * 0.0f - 0.0f * plane_normal_z;
		x_axis_z = plane_normal_x * 1.0f - 1.0f * plane_normal_x;
	}

	float y_axis_x = plane_normal_y * x_axis_z - x_axis_y * plane_normal_z;
	float y_axis_y = plane_normal_x * x_axis_z - x_axis_x * plane_normal_z;
	float y_axis_z = plane_normal_x * x_axis_y - x_axis_x * plane_normal_y;

	float hit_vec_x = intersection_x - plane_origin_x;
	float hit_vec_y = intersection_y - plane_origin_y;
	float hit_vec_z = intersection_z - plane_origin_z;

	out_u =
	(
		hit_vec_x * x_axis_x +
		hit_vec_y * x_axis_y +
		hit_vec_z * x_axis_z
	);

	out_v =
	(
		hit_vec_x * y_axis_x +
		hit_vec_y * y_axis_y +
		hit_vec_z * y_axis_z
	);
}

inline float plane_intersect
(
	plane plane1,

	float ray_ox, float ray_oy, float ray_oz,
	float ray_dx, float ray_dy, float ray_dz,

	float* norm_x,
	float* norm_y,
	float* norm_z,

	float* texture_u,
	float* texture_v
)
{
	float denom =
	(
		-plane1.norm_x * ray_dx +
		-plane1.norm_y * ray_dy +
		-plane1.norm_z * ray_dz
	);

	if (denom > 1e-6f)
	{
		float v_x = plane1.x - ray_ox;
		float v_y = plane1.y - ray_oy;
		float v_z = plane1.z - ray_oz;

		float distance =
		(
			v_x * -plane1.norm_x +
			v_y * -plane1.norm_y +
			v_z * -plane1.norm_z
		);

		float t = distance / denom;

		set_ptr(norm_x, plane1.norm_x);
		set_ptr(norm_y, plane1.norm_y);
		set_ptr(norm_z, plane1.norm_z);

		float temp_texture_u;
		float temp_texture_v;

		plane_uv
		(
			ray_ox + ray_dx * t,
			ray_oy + ray_dy * t,
			ray_oz + ray_dz * t,

			plane1.x,
			plane1.y,
			plane1.z,

			plane1.norm_x,
			plane1.norm_y,
			plane1.norm_z,

			temp_texture_u,
			temp_texture_v
		);

		set_ptr(texture_u, temp_texture_u);
		set_ptr(texture_v, temp_texture_v);

		return t;
	}

	return -1.0f;
}