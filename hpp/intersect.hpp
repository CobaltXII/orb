#pragma once

float do_intersect
(
	float ray_ox,
	float ray_oy,
	float ray_oz,

	float ray_dx,
	float ray_dy,
	float ray_dz,

	float* norm_x,
	float* norm_y,
	float* norm_z,

	float* texture_u,
	float* texture_v,

	shape* shape1
)
{
	if (shape1->primitive == shape_type::st_sphere)
	{
		return sphere_intersect
		(
			TO_SPHERE(shape1),

			ray_ox, ray_oy, ray_oz,
			ray_dx, ray_dy, ray_dz,

			norm_x,
			norm_y,
			norm_z
		);
	}
	else if (shape1->primitive == shape_type::st_plane)
	{
		return plane_intersect
		(
			TO_PLANE(shape1),

			ray_ox, ray_oy, ray_oz,
			ray_dx, ray_dy, ray_dz,

			norm_x,
			norm_y,
			norm_z,

			texture_u,
			texture_v
		);
	}
	else if (shape1->primitive == shape_type::st_ellipsoid)
	{
		return ellipsoid_intersect
		(
			TO_ELLIPSOID(shape1),

			ray_ox, ray_oy, ray_oz,
			ray_dx, ray_dy, ray_dz,

			norm_x,
			norm_y,
			norm_z
		);
	}
	else if (shape1->primitive == shape_type::st_cone)
	{
		return cone_intersect
		(
			TO_CONE(shape1),

			ray_ox, ray_oy, ray_oz,
			ray_dx, ray_dy, ray_dz,

			norm_x,
			norm_y,
			norm_z
		);
	}
	else if (shape1->primitive == shape_type::st_capsule)
	{
		return capsule_intersect
		(
			TO_CAPSULE(shape1),

			ray_ox, ray_oy, ray_oz,
			ray_dx, ray_dy, ray_dz,

			norm_x,
			norm_y,
			norm_z
		);
	}
	else if (shape1->primitive == shape_type::st_cylinder)
	{
		return cylinder_intersect
		(
			TO_CYLINDER(shape1),

			ray_ox, ray_oy, ray_oz,
			ray_dx, ray_dy, ray_dz,

			norm_x,
			norm_y,
			norm_z
		);
	}
	else if (shape1->primitive == shape_type::st_triangle)
	{
		return triangle_intersect
		(
			TO_TRIANGLE(shape1),

			ray_ox, ray_oy, ray_oz,
			ray_dx, ray_dy, ray_dz,

			norm_x,
			norm_y,
			norm_z,

			texture_u,
			texture_v
		);
	}
	else if (shape1->primitive == shape_type::st_intersection)
	{
		float atemp_norm_x, btemp_norm_x;
		float atemp_norm_y, btemp_norm_y;
		float atemp_norm_z, btemp_norm_z;

		float atemp_texture_u, btemp_texture_u;
		float atemp_texture_v, btemp_texture_v;

		float a_t = do_intersect
		(
			ray_ox, ray_oy, ray_oz,
			ray_dx, ray_dy, ray_dz,

			&atemp_norm_x,
			&atemp_norm_y,
			&atemp_norm_z,

			&atemp_texture_u,
			&atemp_texture_v,

			TO_INTERSECTION(shape1).sa
		);

		float b_t = do_intersect
		(
			ray_ox, ray_oy, ray_oz,
			ray_dx, ray_dy, ray_dz,

			&btemp_norm_x,
			&btemp_norm_y,
			&btemp_norm_z,

			&btemp_texture_u,
			&btemp_texture_v,

			TO_INTERSECTION(shape1).sb
		);

		if (a_t > 0.0f && b_t > 0.0f)
		{
			if (a_t > b_t)
			{
				set_ptr(norm_x, atemp_norm_x);
				set_ptr(norm_y, atemp_norm_y);
				set_ptr(norm_z, atemp_norm_z);

				set_ptr(texture_u, atemp_texture_u);
				set_ptr(texture_v, atemp_texture_v);

				shape1->r = TO_INTERSECTION(shape1).sa->r;
				shape1->g = TO_INTERSECTION(shape1).sa->g;
				shape1->b = TO_INTERSECTION(shape1).sa->b;

				shape1->material1 = TO_INTERSECTION(shape1).sa->material1;
				shape1->material2 = TO_INTERSECTION(shape1).sa->material2;
				shape1->material3 = TO_INTERSECTION(shape1).sa->material3;
				shape1->material4 = TO_INTERSECTION(shape1).sa->material4;
				shape1->material5 = TO_INTERSECTION(shape1).sa->material5;

				return a_t;
			}
			else
			{
				set_ptr(norm_x, btemp_norm_x);
				set_ptr(norm_y, btemp_norm_y);
				set_ptr(norm_z, btemp_norm_z);

				set_ptr(texture_u, btemp_texture_u);
				set_ptr(texture_v, btemp_texture_v);

				shape1->r = TO_INTERSECTION(shape1).sb->r;
				shape1->g = TO_INTERSECTION(shape1).sb->g;
				shape1->b = TO_INTERSECTION(shape1).sb->b;

				shape1->material1 = TO_INTERSECTION(shape1).sb->material1;
				shape1->material2 = TO_INTERSECTION(shape1).sb->material2;
				shape1->material3 = TO_INTERSECTION(shape1).sb->material3;
				shape1->material4 = TO_INTERSECTION(shape1).sb->material4;
				shape1->material5 = TO_INTERSECTION(shape1).sb->material5;

				return b_t;
			}
		}
		else
		{
			return -1.0f;
		}
	}

	return -1.0f;
}