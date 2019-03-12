float signed_distance_field
(
	float x,
	float y,
	float z
)
{
	float distance = std::numeric_limits<float>::max();

	for (int i = 0; i < shapes.size(); i++)
	{
		shape* shape1 = shapes[i];

		// Use the appropriate signed distance function.

		if (shape1->primitive == shape_type::st_sphere)
		{
			sphere sphere1 = TO_SPHERE(shape1);

			float p_x = x - sphere1.x;
			float p_y = y - sphere1.y;
			float p_z = z - sphere1.z;

			float p_length = sqrtf
			(
				p_x * p_x +
				p_y * p_y +
				p_z * p_z
			);

			distance = fmin(p_length - sphere1.radius, distance);
		}
		else if (shape1->primitive == shape_type::st_ellipsoid)
		{
			ellipsoid ellipsoid1 = TO_ELLIPSOID(shape1);

			float k0_x = (x - ellipsoid1.x) / ellipsoid1.radius_x;
			float k0_y = (y - ellipsoid1.y) / ellipsoid1.radius_y;
			float k0_z = (z - ellipsoid1.z) / ellipsoid1.radius_z;

			float k1_x = (x - ellipsoid1.x) / (ellipsoid1.radius_x * ellipsoid1.radius_x);
			float k1_y = (y - ellipsoid1.y) / (ellipsoid1.radius_y * ellipsoid1.radius_y);
			float k1_z = (z - ellipsoid1.z) / (ellipsoid1.radius_z * ellipsoid1.radius_z);

			float k0 = sqrtf
			(
				k0_x * k0_x +
				k0_y * k0_y +
				k0_z * k0_z
			);

			float k1 = sqrtf
			(
				k1_x * k1_x +
				k1_y * k1_y +
				k1_z * k1_z
			);

			distance = fmin(k0 * (k0 - 1.0f) / k1, distance);
		}
	}

	return distance;
}