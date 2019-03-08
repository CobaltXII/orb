#define GAMMA

#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include <memory>
#include <limits>
#include <random>

#define fmin(a, b) ((a) < (b) ? (a) : (b))
#define fmax(a, b) ((a) > (b) ? (a) : (b))

const int R = 0;
const int G = 1;
const int B = 2;

#include "stb_image.h"

#include "stb_image_write.h"

#include "stb_perlin.h"

template <typename T>

inline void set_ptr(T* __ptr, T __value)
{
	if (__ptr)
	{
		*__ptr = __value;
	}
}

enum shape_type
{
	st_sphere, st_plane, st_ellipsoid
};

struct shape
{
	float material1;
	float material2;
	float material3;
	float material4;

	float r;
	float g;
	float b;

	shape_type primitive;
};

struct sphere: shape
{
	float x;
	float y;
	float z;

	float radius;

	sphere
	(
		float x,
		float y,
		float z,

		float r,
		float g,
		float b,

		float radius,

		float material1 = 0.0f,
		float material2 = 0.0f,
		float material3 = 0.0f,
		float material4 = 0.0f
	)
	{
		this->primitive = shape_type::st_sphere;

		this->x = x;
		this->y = y;
		this->z = z;

		this->r = r;
		this->g = g;
		this->b = b;

		this->radius = radius;

		this->material1 = material1;
		this->material2 = material2;
		this->material3 = material3;
		this->material4 = material4;
	}
};

inline sphere TO_SPHERE(shape* __victim)
{
	return *((sphere*)__victim);
}

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
		float material4 = 0.0f
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
	}
};

inline ellipsoid TO_ELLIPSOID(shape* __victim)
{
	return *((ellipsoid*)__victim);
}

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
		float material4 = 0.0f
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
	}
};

inline plane TO_PLANE(shape* __victim)
{
	return *((plane*)__victim);
}

inline float rand11()
{
	return float(rand()) / float(RAND_MAX) * 2.0f - 1.0f;
}

inline float rand01()
{
	return float(rand()) / float(RAND_MAX);
}

void nuke(std::string note)
{
	std::cout << note << std::endl;

	exit(EXIT_FAILURE);
}

struct sampler
{
	int x_res;
	int y_res;

	unsigned char* data;

	sampler(std::string path)
	{
		data = stbi_load(path.c_str(), &x_res, &y_res, NULL, 3);

		if (!data)
		{
			nuke("Could not load image (stbi_load).");
		}
	}

	inline void sample
	(
		float u, 
		float v,

		float& out_color_r,
		float& out_color_g,
		float& out_color_b
	)
	{
		int x = fmin(fmax(0, u * x_res), x_res - 1);
		int y = fmin(fmax(0, v * y_res), y_res - 1);

		unsigned char* pixel = data + (y * x_res + x) * 3;

		out_color_r = float(pixel[R]) / 255.0f;
		out_color_g = float(pixel[G]) / 255.0f;
		out_color_b = float(pixel[B]) / 255.0f;
	}
};

std::vector<shape*> shapes;

inline void sphere_uv
(
	float intersection_x,
	float intersection_y,
	float intersection_z,

	float sphere_center_x,
	float sphere_center_y,
	float sphere_center_z,

	float sphere_radius,

	float& out_u,
	float& out_v
)
{
	float hit_vec_x = intersection_x - sphere_center_x;
	float hit_vec_y = intersection_y - sphere_center_y;
	float hit_vec_z = intersection_z - sphere_center_z;

	out_u = (1.0f + atan2f(hit_vec_z, hit_vec_x) / M_PI) * 0.5f;

	out_v = acosf(hit_vec_y / sphere_radius) / M_PI;
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

std::vector<light> lights;

float cast
(
	float ray_ox,
	float ray_oy,
	float ray_oz,

	float ray_dx,
	float ray_dy,
	float ray_dz,

	float* hit_shape_material1 = NULL,
	float* hit_shape_material2 = NULL,
	float* hit_shape_material3 = NULL,
	float* hit_shape_material4 = NULL,

	float* hit_shape_r = NULL,
	float* hit_shape_g = NULL,
	float* hit_shape_b = NULL,

	shape** hit_shape = NULL
)
{
	float min_dist = std::numeric_limits<float>::max();

	for (int i = 0; i < shapes.size(); i++)
	{
		shape* shape1 = shapes[i];

		// Use the appropriate intersector.

		if (shape1->primitive == shape_type::st_sphere)
		{
			sphere sphere1 = TO_SPHERE(shape1);

			float i_lx = sphere1.x - ray_ox;
			float i_ly = sphere1.y - ray_oy;
			float i_lz = sphere1.z - ray_oz;

			float i_adj2 =
			(
				i_lx * ray_dx +
				i_ly * ray_dy +
				i_lz * ray_dz
			);

			float i_d2 =
			(
				i_lx * i_lx +
				i_ly * i_ly +
				i_lz * i_lz
			);

			i_d2 -= i_adj2 * i_adj2;

			float i_r2 =
			(
				sphere1.radius *
				sphere1.radius
			);

			if (i_d2 > i_r2)
			{
				continue;
			}

			float i_thc = sqrtf(i_r2 - i_d2);

			float i_t0 = i_adj2 - i_thc;
			float i_t1 = i_adj2 + i_thc;

			if (i_t0 < 0.0f && i_t1 < 0.0f)
			{
				continue;
			}

			float i_d = std::numeric_limits<float>::max();

			if (i_t0 < 0.0f)
			{
				i_d = i_t1;
			}
			else if (i_t1 < 0.0f)
			{
				i_d = i_t0;
			}
			else
			{
				i_d = fmin(i_t0, i_t1);
			}

			if (i_d < min_dist)
			{
				min_dist = i_d;

				set_ptr(hit_shape_material1, sphere1.material1);
				set_ptr(hit_shape_material2, sphere1.material2);
				set_ptr(hit_shape_material3, sphere1.material3);
				set_ptr(hit_shape_material4, sphere1.material4);

				set_ptr(hit_shape_r, sphere1.r);
				set_ptr(hit_shape_g, sphere1.g);
				set_ptr(hit_shape_b, sphere1.b);

				set_ptr(hit_shape, shape1);
			}
		}
		else if (shape1->primitive == shape_type::st_plane)
		{
			plane plane1 = TO_PLANE(shape1);

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

				if (distance >= 0.0f)
				{
					float i_d = distance / denom;
					
					if (i_d < min_dist)
					{
						min_dist = i_d;

						set_ptr(hit_shape_material1, plane1.material1);
						set_ptr(hit_shape_material2, plane1.material2);
						set_ptr(hit_shape_material3, plane1.material3);
						set_ptr(hit_shape_material4, plane1.material4);

						set_ptr(hit_shape_r, plane1.r);
						set_ptr(hit_shape_g, plane1.g);
						set_ptr(hit_shape_b, plane1.b);

						set_ptr(hit_shape, shape1);
					}
				}
			}
		}
		else if (shape1->primitive == shape_type::st_ellipsoid)
		{
			ellipsoid ellipsoid1 = TO_ELLIPSOID(shape1);

			float oc_x = ray_ox - ellipsoid1.x;
			float oc_y = ray_oy - ellipsoid1.y;
			float oc_z = ray_oz - ellipsoid1.z;

			float ocn_x = oc_x / ellipsoid1.radius_x;
			float ocn_y = oc_y / ellipsoid1.radius_y;
			float ocn_z = oc_z / ellipsoid1.radius_z;

			float rdn_x = ray_dx / ellipsoid1.radius_x;
			float rdn_y = ray_dy / ellipsoid1.radius_y;
			float rdn_z = ray_dz / ellipsoid1.radius_z;

			float a =
			(
				rdn_x * rdn_x +
				rdn_y * rdn_y +
				rdn_z * rdn_z
			);

			float b =
			(
				ocn_x * rdn_x +
				ocn_y * rdn_y +
				ocn_z * rdn_z
			);

			float c =
			(
				ocn_x * ocn_x +
				ocn_y * ocn_y +
				ocn_z * ocn_z
			);

			float h = b * b - a * (c - 1.0f);

			if (h < 0.0f)
			{
				continue;
			}

			float i_d = (-b - sqrtf(h)) / a;

			// This line is kind of buggy, it detects a lot of false
			// collisions but if we set the epsilon too high then it misses a
			// lot of true collisions.
			//
			// Shadows look really sketchy with a big constant like 1e-2f, but
			// you will see artifacts on the tips of tall and reflective
			// ellipsoids if you use something like 1e-3f.
			//
			// I still don't know how to fix this, so I use 0.0f (good
			// shadows) and just don't put large ellipsoids.

			if (i_d < 0.0f)
			{
				continue;
			}

			if (i_d < min_dist)
			{
				min_dist = i_d;

				set_ptr(hit_shape_material1, ellipsoid1.material1);
				set_ptr(hit_shape_material2, ellipsoid1.material2);
				set_ptr(hit_shape_material3, ellipsoid1.material3);
				set_ptr(hit_shape_material4, ellipsoid1.material4);

				set_ptr(hit_shape_r, ellipsoid1.r);
				set_ptr(hit_shape_g, ellipsoid1.g);
				set_ptr(hit_shape_b, ellipsoid1.b);

				set_ptr(hit_shape, shape1);
			}
		}
	}

	return min_dist;
}

void trace
(
	float ray_ox,
	float ray_oy,
	float ray_oz,

	float ray_dx,
	float ray_dy,
	float ray_dz,

	float& out_r,
	float& out_g,
	float& out_b,

	int depth = 0
)
{
	// Cast primary ray.

	float hit_shape_material1;
	float hit_shape_material2;
	float hit_shape_material3;
	float hit_shape_material4;

	float hit_shape_r;
	float hit_shape_g;
	float hit_shape_b;

	shape* hit_shape;

	float min_dist = cast
	(
		ray_ox,
		ray_oy,
		ray_oz,

		ray_dx,
		ray_dy,
		ray_dz,

		&hit_shape_material1,
		&hit_shape_material2,
		&hit_shape_material3,
		&hit_shape_material4,

		&hit_shape_r,
		&hit_shape_g,
		&hit_shape_b,

		&hit_shape
	);

	out_r = 0.0f;
	out_g = 0.0f;
	out_b = 0.0f;

	if (min_dist == std::numeric_limits<float>::max())
	{
		// Missed everything, use background shader.

		int mode = 1;

		if (mode == 0)
		{
			// Checkerboard.

			float ray_intersection_u;
			float ray_intersection_v;

			sphere_uv
			(
				ray_dx,
				ray_dy,
				ray_dz,

				0.0f, 0.0f, 0.0f,

				1.0f,

				ray_intersection_u,
				ray_intersection_v
			);

			float domain_u = 0.025f;
			float domain_v = 0.025f;

			bool check_u = fmod(ray_intersection_u + domain_u * 100.0f, domain_u * 2.0f) >= domain_u;
			bool check_v = fmod(ray_intersection_v + domain_v * 100.0f, domain_v * 2.0f) >= domain_v;

			if (check_u ^ check_v)
			{
				out_r = 0.0f;
				out_g = 0.0f;
				out_b = 0.0f;
			}
			else
			{
				out_r = 0.8f;
				out_g = 0.8f;
				out_b = 0.8f;
			}
		}
		else if (mode == 1)
		{
			// Blood sky.

			float frequency = 8.0f;

			float noise = stb_perlin_fbm_noise3
			(
				ray_dx * frequency, 
				ray_dy * frequency, 
				ray_dz * frequency, 

				2.0f, 0.5f, 6
			);

			out_r = fmax(0.0f, fmin(1.0f, (noise + 1.0f) / 2.0f * 0.032f));
		}

		return;
	}

	// Hit something. Find the color using contributions from light, shadow,
	// reflection rays, etc.

	float hit_x = ray_ox + ray_dx * min_dist;
	float hit_y = ray_oy + ray_dy * min_dist;
	float hit_z = ray_oz + ray_dz * min_dist;

	float norm_x;
	float norm_y;
	float norm_z;

	// Compute surface normal.

	if (hit_shape->primitive == shape_type::st_sphere)
	{
		norm_x = hit_x - TO_SPHERE(hit_shape).x;
		norm_y = hit_y - TO_SPHERE(hit_shape).y;
		norm_z = hit_z - TO_SPHERE(hit_shape).z;
	}
	else if (hit_shape->primitive == shape_type::st_plane)
	{
		norm_x = TO_PLANE(hit_shape).norm_x;
		norm_y = TO_PLANE(hit_shape).norm_y;
		norm_z = TO_PLANE(hit_shape).norm_z;
	}
	else if (hit_shape->primitive == shape_type::st_ellipsoid)
	{
		norm_x = (hit_x - TO_ELLIPSOID(hit_shape).x) / TO_ELLIPSOID(hit_shape).radius_x;
		norm_y = (hit_y - TO_ELLIPSOID(hit_shape).y) / TO_ELLIPSOID(hit_shape).radius_y;
		norm_z = (hit_z - TO_ELLIPSOID(hit_shape).z) / TO_ELLIPSOID(hit_shape).radius_z;
	}

	float norm_length = sqrtf
	(
		norm_x * norm_x +
		norm_y * norm_y +
		norm_z * norm_z
	);

	norm_x /= norm_length;
	norm_y /= norm_length;
	norm_z /= norm_length;

	// Small spheres should be procedurally textured.

	if (hit_shape->primitive == shape_type::st_sphere || hit_shape->primitive == shape_type::st_ellipsoid)
	{
		float frequency = 0.48f;

		float noise = stb_perlin_ridge_noise3
		(
			hit_x * frequency,
			hit_y * frequency,
			hit_z * frequency,

			2.0f, 0.5f, 1.0f, 6
		);

		if (hit_shape_r > 0.5f)
		{
			hit_shape_r = (noise + 1.0f) / 2.0f;
		}
		else if (hit_shape_g > 0.5f)
		{
			hit_shape_g = (noise + 1.0f) / 2.0f;
		}
		else if (hit_shape_b > 0.5f)
		{
			hit_shape_b = (noise + 1.0f) / 2.0f;
		}
	}

	// Check for planes, which should be checkered.

	if (hit_shape->primitive == shape_type::st_plane)
	{
		float plane_u;
		float plane_v;

		plane_uv
		(
			hit_x,
			hit_y,
			hit_z,

			TO_PLANE(hit_shape).x,
			TO_PLANE(hit_shape).y,
			TO_PLANE(hit_shape).z,

			norm_x,
			norm_y,
			norm_z,

			plane_u,
			plane_v
		);

		float domain_u = 10.0f;
		float domain_v = 10.0f;

		float expand = 16384.0f;

		int check_u = fmod(plane_u + domain_u * expand, domain_u * 2.0f) >= domain_u;
		int check_v = fmod(plane_v + domain_v * expand, domain_v * 2.0f) >= domain_v;

		if (check_u ^ check_v)
		{
			hit_shape_r = 0.0f * 2.0f;
			hit_shape_g = 0.0f * 2.0f;
			hit_shape_b = 0.0f * 2.0f;
		}
		else
		{
			hit_shape_r = 0.85f * 2.0f;
			hit_shape_g = 0.85f * 2.0f;
			hit_shape_b = 0.85f * 2.0f;
		}
	}

	// Calculate lighting.

	for (int i = 0; i < lights.size(); i++)
	{
		light light1 = lights[i];

		// Direction to shape.

		float dtl_x = light1.x - hit_x;
		float dtl_y = light1.y - hit_y;
		float dtl_z = light1.z - hit_z;

		float dr = sqrtf
		(
			dtl_x * dtl_x +
			dtl_y * dtl_y +
			dtl_z * dtl_z
		);

		dtl_x /= dr;
		dtl_y /= dr;
		dtl_z /= dr;

		// Send shadow ray.

		float shadow_dist = cast
		(
			hit_x + 1e-1f * dtl_x,
			hit_y + 1e-1f * dtl_y,
			hit_z + 1e-1f * dtl_z,

			dtl_x,
			dtl_y,
			dtl_z
		);

		if (shadow_dist < dr)
		{
			#ifdef VERY_HARD_SHADOWS

			return;

			#else

			continue;

			#endif
		}

		// Diffuse.

		float l_pow =
		(
			norm_x * dtl_x +
			norm_y * dtl_y +
			norm_z * dtl_z
		);

		l_pow /= dr * dr;

		out_r += fmax(0.0f, hit_shape_r * light1.r * l_pow);
		out_g += fmax(0.0f, hit_shape_g * light1.g * l_pow);
		out_b += fmax(0.0f, hit_shape_b * light1.b * l_pow);

		// Specular.

		float specular_constant = 0.5f;

		float view_hit_x = ray_ox - hit_x;
		float view_hit_y = ray_oy - hit_y;
		float view_hit_z = ray_oz - hit_z;

		float view_hit_length = sqrtf
		(
			view_hit_x * view_hit_x +
			view_hit_y * view_hit_y +
			view_hit_z * view_hit_z
		);

		view_hit_x /= view_hit_length;
		view_hit_y /= view_hit_length;
		view_hit_z /= view_hit_length;

		float dot_n_i =
		(
			-dtl_x * norm_x +
			-dtl_y * norm_y +
			-dtl_z * norm_z
		);

		float specular_x = -dtl_x - 2.0f * dot_n_i * norm_x;
		float specular_y = -dtl_y - 2.0f * dot_n_i * norm_y;
		float specular_z = -dtl_z - 2.0f * dot_n_i * norm_z;

		float dot_view_specular =
		(
			view_hit_x * specular_x +
			view_hit_y * specular_y +
			view_hit_z * specular_z
		);

		float specular_coefficient = powf(fmax(dot_view_specular, 0.0f), hit_shape_material4);

		out_r += fmax(0.0f, light1.r * specular_coefficient);
		out_g += fmax(0.0f, light1.g * specular_coefficient);
		out_b += fmax(0.0f, light1.b * specular_coefficient);
	}

	// Prevent infinite recursion.

	if (depth > 5)
	{
		return;
	}

	// Add reflections.

	if (hit_shape_material1 > 0.0f)
	{
		float eps = 1e-3f;

		float reflect_ox = hit_x + norm_x * eps;
		float reflect_oy = hit_y + norm_y * eps;
		float reflect_oz = hit_z + norm_z * eps;

		float incident_dot_norm =
		(
			ray_dx * norm_x +
			ray_dy * norm_y +
			ray_dz * norm_z
		);

		float reflect_dx = ray_dx - (2.0f * incident_dot_norm * norm_x);
		float reflect_dy = ray_dy - (2.0f * incident_dot_norm * norm_y);
		float reflect_dz = ray_dz - (2.0f * incident_dot_norm * norm_z);

		float reflect_r = 0.0f;
		float reflect_g = 0.0f;
		float reflect_b = 0.0f;

		if (false)
		{
			// Glossy reflections.

			float reflection_rays = 1.0f;

			for (int i = 0; i < reflection_rays; i++)
			{
				float gloss = 0.0f;

				float additive_reflect_r;
				float additive_reflect_g;
				float additive_reflect_b;

				trace
				(
					reflect_ox,
					reflect_oy,
					reflect_oz,

					reflect_dx + rand11() * gloss,
					reflect_dy + rand11() * gloss,
					reflect_dz + rand11() * gloss,

					additive_reflect_r,
					additive_reflect_g,
					additive_reflect_b,

					depth + 1
				);

				reflect_r += additive_reflect_r;
				reflect_g += additive_reflect_g;
				reflect_b += additive_reflect_b;
			}

			reflect_r /= reflection_rays;
			reflect_g /= reflection_rays;
			reflect_b /= reflection_rays;
		}
		else
		{
			// Standard reflections.

			trace
			(
				reflect_ox,
				reflect_oy,
				reflect_oz,

				reflect_dx,
				reflect_dy,
				reflect_dz,

				reflect_r,
				reflect_g,
				reflect_b,

				depth + 1
			);
		}

		out_r = out_r * (1.0f - hit_shape_material1) + reflect_r * hit_shape_material1;
		out_g = out_g * (1.0f - hit_shape_material1) + reflect_g * hit_shape_material1;
		out_b = out_b * (1.0f - hit_shape_material1) + reflect_b * hit_shape_material1;
	}

	// Add refractions.

	if (hit_shape_material2 > 0.0f)
	{
		float ref_n_x = norm_x;
		float ref_n_y = norm_y;
		float ref_n_z = norm_z;

		float eta_i = 1.0f;

		float i_dot_n =
		(
			ray_dx * norm_x +
			ray_dy * norm_y +
			ray_dz * norm_z
		);

		float eta_t = hit_shape_material3;

		if (i_dot_n < 0.0f)
		{
			i_dot_n = -i_dot_n;
		}
		else
		{
			eta_i = hit_shape_material3;

			ref_n_x = -norm_x;
			ref_n_y = -norm_y;
			ref_n_z = -norm_z;

			eta_t = 1.0f;
		}

		float eta = eta_i / eta_t;

		float refract_r = 0.0f;
		float refract_g = 0.0f;
		float refract_b = 0.0f;

		float k = sqrtf(1.0f - (eta * eta) * (1.0f - i_dot_n * i_dot_n));

		if (k * k >= 0.0f)
		{
			float eps = 1e-3f;

			float refract_ox = hit_x + ref_n_x * -eps;
			float refract_oy = hit_y + ref_n_y * -eps;
			float refract_oz = hit_z + ref_n_z * -eps;

			float refract_dx = (ray_dx + i_dot_n * ref_n_x) * eta - ref_n_x * k;
			float refract_dy = (ray_dy + i_dot_n * ref_n_y) * eta - ref_n_y * k;
			float refract_dz = (ray_dz + i_dot_n * ref_n_z) * eta - ref_n_z * k;

			trace
			(
				refract_ox,
				refract_oy,
				refract_oz,

				refract_dx,
				refract_dy,
				refract_dz,

				refract_r,
				refract_g,
				refract_b,

				depth + 1
			);
		}	

		out_r = out_r * (1.0f - hit_shape_material2) + refract_r * hit_shape_material2;
		out_g = out_g * (1.0f - hit_shape_material2) + refract_g * hit_shape_material2;
		out_b = out_b * (1.0f - hit_shape_material2) + refract_b * hit_shape_material2;
	}

	#ifndef REALLY_BRIGHT_LIGHTS

	out_r = fmax(0.0f, fmin(1.0f, out_r));
	out_g = fmax(0.0f, fmin(1.0f, out_g));
	out_b = fmax(0.0f, fmin(1.0f, out_b));

	#endif
}

int main(int argc, char** argv)
{
	shapes.push_back(new plane(0.0f, 0.0f - 24.0f, 0.0f, 0.0f, 0.0f + 1.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.000f, 0.0f, 0.0f, 2048.0f));
	shapes.push_back(new plane(0.0f, 0.0f + 64.0f, 0.0f, 0.0f, 0.0f - 1.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.000f, 0.0f, 0.0f, 2048.0f));

	// shapes.push_back(new sphere(0.0f - 32.0f * 1.0f, -8.0f, -56.0f, 1.000f, 0.000f, 0.000f, 16.0f, 0.200f, 0.0f, 0.0f, 2048.0f));
	// shapes.push_back(new sphere(0.0f + 32.0f * 0.0f, -8.0f, -56.0f, 0.000f, 1.000f, 0.000f, 16.0f, 0.200f, 0.0f, 0.0f, 2048.0f));
	// shapes.push_back(new sphere(0.0f + 32.0f * 1.0f, -8.0f, -56.0f, 0.000f, 0.000f, 1.000f, 16.0f, 0.200f, 0.0f, 0.0f, 2048.0f));

	shapes.push_back(new ellipsoid(0.0f - 32.0f * 1.0f, -8.0f, -56.0f, 1.000f, 0.000f, 0.000f, 16.0f, 24.0f, 16.0f, 0.200f, 0.0f, 0.0f, 2048.0f));
	shapes.push_back(new ellipsoid(0.0f + 32.0f * 0.0f, -8.0f, -56.0f, 0.000f, 1.000f, 0.000f, 16.0f, 24.0f, 16.0f, 0.200f, 0.0f, 0.0f, 2048.0f));
	shapes.push_back(new ellipsoid(0.0f + 32.0f * 1.0f, -8.0f, -56.0f, 0.000f, 0.000f, 1.000f, 16.0f, 24.0f, 16.0f, 0.200f, 0.0f, 0.0f, 2048.0f));

	// shapes.push_back(new ellipsoid(0.0f + 32.0f * 0.0f,32.0f, -128.0f, 0.000f, 1.000f, 0.000f, 320000.0f, 48.0f, 48.0f, 0.200f, 0.0f, 0.0f, 2048.0f));

	lights.push_back(light(25.0f, 50.0f, 0.0f, 1e+3f * 5.6f, 1e+3f * 5.6f, 1e+3f * 5.6f, 50.0f));

	int supersample = 3;

	int x_res = (128 * 8) * supersample;
	int y_res = (128 * 8) * supersample;

	float x_resf = x_res;
	float y_resf = y_res;

	float fov = 90.0f;

	float aspect =
	(
		x_resf /
		y_resf
	);

	float fov_adjust = tanf(fov * M_PI / 360.0f);

	unsigned char* frame_buffer = (unsigned char*)malloc(x_res * y_res * 3 * sizeof(unsigned char));

	if (!frame_buffer)
	{
		nuke("Could not allocate a frame buffer (malloc).");
	}

	for (int i = 0; i < x_res * y_res; i++)
	{
		unsigned char* pixel = frame_buffer + i * 3;

		pixel[R] = 0.0f * 255.0f;
		pixel[G] = 0.0f * 255.0f;
		pixel[B] = 0.0f * 255.0f;
	}

	float gamma = 1.0f / 2.2f;

	for (int j = 0; j < y_res; j++)
	{
		std::cout << "Row " << j + 1 << "/" << y_res << "\r" << std::flush;

		for (int i = 0; i < x_res; i++)
		{
			unsigned char* pixel = frame_buffer + (j * x_res + i) * 3;

			// Generate prime ray.

			float ray_ox = 0.0f;
			float ray_oy = 0.0f;
			float ray_oz = 0.0f;

			float ray_dx = (0.0f + (float(i + 0.5f) / x_resf) * 2.0f - 1.0f) * fov_adjust;
			float ray_dy = (1.0f - (float(j + 0.5f) / y_resf) * 2.0f + 0.0f) * fov_adjust;

			float ray_dz = -1.0f;

			ray_dx *= aspect;

			float ray_length = sqrtf
			(
				ray_dx * ray_dx +
				ray_dy * ray_dy +
				ray_dz * ray_dz
			);

			ray_dx /= ray_length;
			ray_dy /= ray_length;
			ray_dz /= ray_length;

			float color_r = 0.0f;
			float color_g = 0.0f;
			float color_b = 0.0f;

			trace
			(
				ray_ox,
				ray_oy,
				ray_oz,

				ray_dx,
				ray_dy,
				ray_dz,

				color_r,
				color_g,
				color_b
			);

			#ifdef GAMMA

			// Gamma correction (slightly broken).
			
			pixel[R] = fmax(0.0f, fmin(255.0f, powf(color_r, gamma) * 255.0f));
			pixel[G] = fmax(0.0f, fmin(255.0f, powf(color_g, gamma) * 255.0f));
			pixel[B] = fmax(0.0f, fmin(255.0f, powf(color_b, gamma) * 255.0f));

			#else

			pixel[R] = fmax(0.0f, fmin(255.0f, color_r * 255.0f));
			pixel[G] = fmax(0.0f, fmin(255.0f, color_g * 255.0f));
			pixel[B] = fmax(0.0f, fmin(255.0f, color_b * 255.0f));

			#endif
		}
	}

	if (supersample > 1)
	{
		unsigned char* supersampled = (unsigned char*)malloc
		(
			x_res / supersample *
			y_res / supersample *

			3 * sizeof(unsigned char)
		);

		if (!supersample)
		{
			nuke("Could not allocate a frame buffer for supersampling (malloc).");
		}

		x_res = x_res / supersample;
		y_res = y_res / supersample;

		for (int j = 0; j < y_res; j++)
		for (int i = 0; i < x_res; i++)
		{
			unsigned char* pixel = supersampled + (j * x_res + i) * 3;

			int ssi = i * supersample;
			int ssj = j * supersample;

			float sum_r = 0.0f;
			float sum_g = 0.0f;
			float sum_b = 0.0f;

			for (int v = 0; v < supersample; v++)
			for (int u = 0; u < supersample; u++)
			{
				unsigned char* sspixel = frame_buffer + ((ssj + v) * (x_res * supersample) + (ssi + u)) * 3;

				sum_r += sspixel[R];
				sum_g += sspixel[G];
				sum_b += sspixel[B];
			}

			sum_r /= supersample * supersample;
			sum_g /= supersample * supersample;
			sum_b /= supersample * supersample;

			pixel[R] = fmin(255, fmax(0, int(sum_r)));
			pixel[G] = fmin(255, fmax(0, int(sum_g)));
			pixel[B] = fmin(255, fmax(0, int(sum_b)));
		}

		frame_buffer = supersampled;
	}

	if (!stbi_write_png("orb.png", x_res, y_res, 3, frame_buffer, x_res * 3))
	{
		nuke("Could not save a frame buffer (stbi_write_png).");
	}

	return EXIT_SUCCESS;
}