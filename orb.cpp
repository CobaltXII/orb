#define GAMMA

#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include <memory>
#include <limits>
#include <random>
#include <map>

#include "inih/ini.h"

#define sign(x) ((x > 0.0f) - (x < 0.0f))

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
	st_sphere, st_plane, st_ellipsoid, st_cone, st_capsule, st_cylinder
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
		float material4 = 0.0f
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
	}
};

inline cone TO_CONE(shape* __victim)
{
	return *((cone*)__victim);
}

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
		float material4 = 0.0f
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
	}
};

inline capsule TO_CAPSULE(shape* __victim)
{
	return *((capsule*)__victim);
}

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
		float material4 = 0.0f
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
	}
};

inline cylinder TO_CYLINDER(shape* __victim)
{
	return *((cylinder*)__victim);
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

	float* norm_x = NULL,
	float* norm_y = NULL,
	float* norm_z = NULL,

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

				// Surface normal.

				set_ptr(norm_x, (ray_ox + ray_dx * i_d) - sphere1.x);
				set_ptr(norm_y, (ray_oy + ray_dy * i_d) - sphere1.y);
				set_ptr(norm_z, (ray_oz + ray_dz * i_d) - sphere1.z);
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

						// Surface normal.

						set_ptr(norm_x, plane1.norm_x);
						set_ptr(norm_y, plane1.norm_y);
						set_ptr(norm_z, plane1.norm_z);
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
				continue;
			}

			d = sqrtf(d);

			float t1 = (-b + d) / (2.0f * a);
			float t2 = (-b - d) / (2.0f * a);

			float eps = 0.0f;

			if (t1 <= eps && t2 <= eps)
			{
				continue;
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
				continue;
			}

			float i_d = t;

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

				// Surface normal.

				set_ptr(norm_x, ((ray_ox + ray_dx * i_d) - ellipsoid1.x) / ellipsoid1.radius_x);
				set_ptr(norm_y, ((ray_oy + ray_dy * i_d) - ellipsoid1.y) / ellipsoid1.radius_y);
				set_ptr(norm_z, ((ray_oz + ray_dz * i_d) - ellipsoid1.z) / ellipsoid1.radius_z);
			}
		}
		else if (shape1->primitive == shape_type::st_cone)
		{
			cone cone1 = TO_CONE(shape1);

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
					float t = -oaba / rdba;

					if (t < min_dist)
					{
						min_dist = t;

						set_ptr(hit_shape_material1, cone1.material1);
						set_ptr(hit_shape_material2, cone1.material2);
						set_ptr(hit_shape_material3, cone1.material3);
						set_ptr(hit_shape_material4, cone1.material4);

						set_ptr(hit_shape_r, cone1.r);
						set_ptr(hit_shape_g, cone1.g);
						set_ptr(hit_shape_b, cone1.b);

						set_ptr(hit_shape, shape1);

						// Surface normal.

						float isqr_baba = 1.0f / sqrtf(baba);

						set_ptr(norm_x, -ba_x * isqr_baba);
						set_ptr(norm_y, -ba_y * isqr_baba);
						set_ptr(norm_z, -ba_z * isqr_baba);
					}
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
					if (t < min_dist)
					{
						min_dist = t;

						set_ptr(hit_shape_material1, cone1.material1);
						set_ptr(hit_shape_material2, cone1.material2);
						set_ptr(hit_shape_material3, cone1.material3);
						set_ptr(hit_shape_material4, cone1.material4);

						set_ptr(hit_shape_r, cone1.r);
						set_ptr(hit_shape_g, cone1.g);
						set_ptr(hit_shape_b, cone1.b);

						set_ptr(hit_shape, shape1);

						// Surface normal.

						float isqr_baba = 1.0f / sqrtf(baba);

						set_ptr(norm_x, ba_x * isqr_baba);
						set_ptr(norm_y, ba_y * isqr_baba);
						set_ptr(norm_z, ba_z * isqr_baba);
					}
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
				continue;
			}

			float t = (-k1 - sign(rbra) * sqrtf(h)) / (k2 * rbra);

			if (t < 0.0f)
			{
				continue;
			}

			float y = oaba + t * rdba;

			if (y > 0.0f && y < baba)
			{
				if (t < min_dist)
				{
					min_dist = t;

					set_ptr(hit_shape_material1, cone1.material1);
					set_ptr(hit_shape_material2, cone1.material2);
					set_ptr(hit_shape_material3, cone1.material3);
					set_ptr(hit_shape_material4, cone1.material4);

					set_ptr(hit_shape_r, cone1.r);
					set_ptr(hit_shape_g, cone1.g);
					set_ptr(hit_shape_b, cone1.b);

					set_ptr(hit_shape, shape1);

					// Surface normal.

					float isqr_baba = 1.0f / sqrtf(baba);

					set_ptr(norm_x, baba * (baba * (oa_x + t * ray_dx) - rbra * ba_x * cone1.radius_a) - ba_x * hyhy * y);
					set_ptr(norm_y, baba * (baba * (oa_y + t * ray_dy) - rbra * ba_y * cone1.radius_a) - ba_y * hyhy * y);
					set_ptr(norm_z, baba * (baba * (oa_z + t * ray_dz) - rbra * ba_z * cone1.radius_a) - ba_z * hyhy * y);
				}
			}
		}
		else if (shape1->primitive == shape_type::st_capsule)
		{
			capsule capsule1 = TO_CAPSULE(shape1);

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
					if (t < 0.0f)
					{
						continue;
					}

					if (t < min_dist)
					{
						min_dist = t;

						set_ptr(hit_shape_material1, capsule1.material1);
						set_ptr(hit_shape_material2, capsule1.material2);
						set_ptr(hit_shape_material3, capsule1.material3);
						set_ptr(hit_shape_material4, capsule1.material4);

						set_ptr(hit_shape_r, capsule1.r);
						set_ptr(hit_shape_g, capsule1.g);
						set_ptr(hit_shape_b, capsule1.b);

						set_ptr(hit_shape, shape1);

						// Surface normal.

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
					}
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

					if (t < 0.0f)
					{
						continue;
					}

					if (t < min_dist)
					{
						min_dist = t;

						set_ptr(hit_shape_material1, capsule1.material1);
						set_ptr(hit_shape_material2, capsule1.material2);
						set_ptr(hit_shape_material3, capsule1.material3);
						set_ptr(hit_shape_material4, capsule1.material4);

						set_ptr(hit_shape_r, capsule1.r);
						set_ptr(hit_shape_g, capsule1.g);
						set_ptr(hit_shape_b, capsule1.b);

						set_ptr(hit_shape, shape1);

						// Surface normal.

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
					}
				}
			}
		}
		else if (shape1->primitive == shape_type::st_cylinder)
		{
			cylinder cylinder1 = TO_CYLINDER(shape1);

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
				continue;
			}

			h = sqrtf(h);

			float t = (-b - h) / a;

			float y = caoc + t * card;

			if (y > 0.0f && y < caca)
			{
				if (t < 0.0f)
				{
					continue;
				}

				if (t < min_dist)
				{
					min_dist = t;

					set_ptr(hit_shape_material1, cylinder1.material1);
					set_ptr(hit_shape_material2, cylinder1.material2);
					set_ptr(hit_shape_material3, cylinder1.material3);
					set_ptr(hit_shape_material4, cylinder1.material4);

					set_ptr(hit_shape_r, cylinder1.r);
					set_ptr(hit_shape_g, cylinder1.g);
					set_ptr(hit_shape_b, cylinder1.b);

					set_ptr(hit_shape, shape1);

					// Surface normal.

					set_ptr(norm_x, (oc_x + t * ray_dx - ca_x * y / caca) / cylinder1.radius);
					set_ptr(norm_y, (oc_y + t * ray_dy - ca_y * y / caca) / cylinder1.radius);
					set_ptr(norm_z, (oc_z + t * ray_dz - ca_z * y / caca) / cylinder1.radius);
				}
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
				if (t < 0.0f)
				{
					continue;
				}

				if (t < min_dist)
				{
					min_dist = t;

					set_ptr(hit_shape_material1, cylinder1.material1);
					set_ptr(hit_shape_material2, cylinder1.material2);
					set_ptr(hit_shape_material3, cylinder1.material3);
					set_ptr(hit_shape_material4, cylinder1.material4);

					set_ptr(hit_shape_r, cylinder1.r);
					set_ptr(hit_shape_g, cylinder1.g);
					set_ptr(hit_shape_b, cylinder1.b);

					set_ptr(hit_shape, shape1);

					// Surface normal.

					set_ptr(norm_x, ca_x * sign(y) / caca);
					set_ptr(norm_y, ca_y * sign(y) / caca);
					set_ptr(norm_z, ca_z * sign(y) / caca);
				}
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

	float norm_x;
	float norm_y;
	float norm_z;

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

		&norm_x,
		&norm_y,
		&norm_z,

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

	// Normalize the surface normal.

	float norm_length = sqrtf
	(
		norm_x * norm_x +
		norm_y * norm_y +
		norm_z * norm_z
	);

	norm_x /= norm_length;
	norm_y /= norm_length;
	norm_z /= norm_length;

	// Hit something. Find the color using contributions from light, shadow,
	// reflection rays, etc.

	float hit_x = ray_ox + ray_dx * min_dist;
	float hit_y = ray_oy + ray_dy * min_dist;
	float hit_z = ray_oz + ray_dz * min_dist;

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
	else
	{
		// Everything else gets a procedural texture.

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

std::map<std::string, std::map<std::string, std::string>> ini_file;

int inii(std::string section, std::string name)
{
	return std::stoi(ini_file.at(section).at(name));
}

float inif(std::string section, std::string name)
{
	return std::stof(ini_file.at(section).at(name));
}

std::string inis(std::string section, std::string name)
{
	return ini_file.at(section).at(name);
}

int ini_parser(void* user, const char* section, const char* name, const char* value)
{
	ini_file[std::string(section)][std::string(name)] = value;

	return 1;
};

int main(int argc, char** argv)
{
	if (argc != 2)
	{
		std::cout << "Usage: " << argv[0] << " <scene>" << std::endl;

		exit(EXIT_FAILURE);
	}

	if (ini_parse(argv[1], &ini_parser, NULL) < 0)
	{
		nuke("Could not load scene (ini_parse).");
	}

	for (auto section: ini_file)
	{
		if (section.first != "export" && section.first != "camera")
		{
			std::string type = inis(section.first, "type");

			if (type == "sphere")
			{
				shapes.push_back
				(
					new sphere
					(
						inif(section.first, "x"),
						inif(section.first, "y"),
						inif(section.first, "z"),

						inif(section.first, "r"),
						inif(section.first, "g"),
						inif(section.first, "b"),

						inif(section.first, "radius"),

						inif(section.first, "material1"),
						inif(section.first, "material2"),
						inif(section.first, "material3"),
						inif(section.first, "material4")
					)
				);
			}
			else if (type == "ellipsoid")
			{
				shapes.push_back
				(
					new ellipsoid
					(
						inif(section.first, "x"),
						inif(section.first, "y"),
						inif(section.first, "z"),

						inif(section.first, "r"),
						inif(section.first, "g"),
						inif(section.first, "b"),

						inif(section.first, "radius_x"),
						inif(section.first, "radius_y"),
						inif(section.first, "radius_z"),

						inif(section.first, "material1"),
						inif(section.first, "material2"),
						inif(section.first, "material3"),
						inif(section.first, "material4")
					)
				);
			}
			else if (type == "cone")
			{
				shapes.push_back
				(
					new cone
					(
						inif(section.first, "a_x"),
						inif(section.first, "a_y"),
						inif(section.first, "a_z"),

						inif(section.first, "b_x"),
						inif(section.first, "b_y"),
						inif(section.first, "b_z"),

						inif(section.first, "r"),
						inif(section.first, "g"),
						inif(section.first, "b"),

						inif(section.first, "radius_a"),
						inif(section.first, "radius_b"),

						inif(section.first, "material1"),
						inif(section.first, "material2"),
						inif(section.first, "material3"),
						inif(section.first, "material4")
					)
				);
			}
			else if (type == "capsule")
			{
				shapes.push_back
				(
					new capsule
					(
						inif(section.first, "a_x"),
						inif(section.first, "a_y"),
						inif(section.first, "a_z"),

						inif(section.first, "b_x"),
						inif(section.first, "b_y"),
						inif(section.first, "b_z"),

						inif(section.first, "r"),
						inif(section.first, "g"),
						inif(section.first, "b"),

						inif(section.first, "radius"),

						inif(section.first, "material1"),
						inif(section.first, "material2"),
						inif(section.first, "material3"),
						inif(section.first, "material4")
					)
				);
			}
			else if (type == "cylinder")
			{
				shapes.push_back
				(
					new cylinder
					(
						inif(section.first, "a_x"),
						inif(section.first, "a_y"),
						inif(section.first, "a_z"),

						inif(section.first, "b_x"),
						inif(section.first, "b_y"),
						inif(section.first, "b_z"),

						inif(section.first, "r"),
						inif(section.first, "g"),
						inif(section.first, "b"),

						inif(section.first, "radius"),

						inif(section.first, "material1"),
						inif(section.first, "material2"),
						inif(section.first, "material3"),
						inif(section.first, "material4")
					)
				);
			}
			else if (type == "plane")
			{
				shapes.push_back
				(
					new plane
					(
						inif(section.first, "x"),
						inif(section.first, "y"),
						inif(section.first, "z"),

						inif(section.first, "norm_x"),
						inif(section.first, "norm_y"),
						inif(section.first, "norm_z"),

						inif(section.first, "r"),
						inif(section.first, "g"),
						inif(section.first, "b"),

						inif(section.first, "material1"),
						inif(section.first, "material2"),
						inif(section.first, "material3"),
						inif(section.first, "material4")
					)
				);
			}
			else if (type == "light")
			{
				lights.push_back
				(
					light
					(
						inif(section.first, "x"),
						inif(section.first, "y"),
						inif(section.first, "z"),

						inif(section.first, "r"),
						inif(section.first, "g"),
						inif(section.first, "b"),

						inif(section.first, "radius")
					)
				);
			}
		}
	}

	int supersample = inii("export", "supersample");

	int x_res = inii("export", "x_res") * supersample;
	int y_res = inii("export", "y_res") * supersample;

	float ray_ox = inif("camera", "ray_ox");
	float ray_oy = inif("camera", "ray_oy");
	float ray_oz = inif("camera", "ray_oz");

	float fov = inif("camera", "fov");

	float x_resf = x_res;
	float y_resf = y_res;

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