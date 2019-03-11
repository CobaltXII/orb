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
			norm_z
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

	return -1.0f;
}