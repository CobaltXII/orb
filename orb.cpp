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

#define STB_IMAGE_WRITE_IMPLEMENTATION

#include "stb_image_write.h"

struct orb
{
	float x;
	float y;
	float z;

	float r;
	float g;
	float b;

	float radius;
};

typedef orb lamp;

float rand11()
{
	return float(rand()) / float(RAND_MAX) * 2.0f - 1.0f;
}

float rand01()
{
	return float(rand()) / float(RAND_MAX);
}

void nuke(std::string note)
{
	std::cout << note << std::endl;

	exit(EXIT_FAILURE);
}

std::vector<orb> orbs;

std::vector<lamp> lamps;

float cast
(
	float ray_ox,
	float ray_oy,
	float ray_oz,

	float ray_dx,
	float ray_dy,
	float ray_dz,

	float& hit_orb_x,
	float& hit_orb_y,
	float& hit_orb_z,

	float& hit_orb_r,
	float& hit_orb_g,
	float& hit_orb_b,

	float& hit_orb_radius
)
{
	float min_dist = std::numeric_limits<float>::max();

	for (int i = 0; i < orbs.size(); i++)
	{
		orb orb = orbs[i];

		float i_lx = orb.x - ray_ox;
		float i_ly = orb.y - ray_oy;
		float i_lz = orb.z - ray_oz;

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
			orb.radius *
			orb.radius
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

		float i_d = std::min(i_t0, i_t1);

		if (i_d < min_dist)
		{
			min_dist = i_d;

			hit_orb_x = orb.x;
			hit_orb_y = orb.y;
			hit_orb_z = orb.z;

			hit_orb_r = orb.r;
			hit_orb_g = orb.g;
			hit_orb_b = orb.b;

			hit_orb_radius = orb.radius;
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

	float hit_orb_x;
	float hit_orb_y;
	float hit_orb_z;

	float hit_orb_r;
	float hit_orb_g;
	float hit_orb_b;

	float hit_orb_radius;

	float min_dist = cast
	(
		ray_ox,
		ray_oy,
		ray_oz,

		ray_dx,
		ray_dy,
		ray_dz,

		hit_orb_x,
		hit_orb_y,
		hit_orb_z,

		hit_orb_r,
		hit_orb_g,
		hit_orb_b,

		hit_orb_radius
	);

	out_r = 0.0f;
	out_g = 0.0f;
	out_b = 0.0f;

	if (min_dist == std::numeric_limits<float>::max())
	{
		// Missed everything.

		return;
	}

	// Hit something. Find the color using contributions from light, shadow,
	// reflection rays, etc.

	float hit_x = ray_ox + ray_dx * min_dist;
	float hit_y = ray_oy + ray_dy * min_dist;
	float hit_z = ray_oz + ray_dz * min_dist;

	float norm_x = hit_x - hit_orb_x;
	float norm_y = hit_y - hit_orb_y;
	float norm_z = hit_z - hit_orb_z;

	float norm_length = sqrtf
	(
		norm_x * norm_x +
		norm_y * norm_y +
		norm_z * norm_z
	);

	norm_x /= norm_length;
	norm_y /= norm_length;
	norm_z /= norm_length;

	// Check for large spheres, which should be checkered.

	if (hit_orb_radius > 200.0f)
	{
		float domain_x = 10.0f;
		float domain_y = 10.0f;
		float domain_z = 10.0f;

		float expand = 512.0f;

		int check_x = fmod(hit_x + domain_x * expand, domain_x * 2.0f) >= domain_x;
		int check_y = fmod(hit_y + domain_y * expand, domain_y * 2.0f) >= domain_y;
		int check_z = fmod(hit_z + domain_z * expand, domain_z * 2.0f) >= domain_z;

		if (check_x ^ check_y ^ check_z)
		{
			hit_orb_r = 0.0f * 2.0f;
			hit_orb_g = 0.0f * 2.0f;
			hit_orb_b = 0.0f * 2.0f;
		}
		else
		{
			hit_orb_r = 0.85f * 2.0f;
			hit_orb_g = 0.85f * 2.0f;
			hit_orb_b = 0.85f * 2.0f;
		}
	}

	// Calculate lighting.

	float shadow_rays = 16.0f;

	for (int k = 0; k < int(shadow_rays); k++)
	{
		for (int i = 0; i < lamps.size(); i++)
		{
			lamp lamp = lamps[i];

			// Perturb lamp for smooth shadows.

			float r1 = rand11();
			float r2 = rand11();
			float r3 = rand11();

			float deviation = shadow_rays - 1.0f;

			lamp.x += r1 * deviation;
			lamp.y += r2 * deviation;
			lamp.z += r3 * deviation;

			// Direction to lamp.

			float dtl_x = lamp.x - hit_x;
			float dtl_y = lamp.y - hit_y;
			float dtl_z = lamp.z - hit_z;

			lamp.x -= r1 * deviation;
			lamp.y -= r2 * deviation;
			lamp.z -= r3 * deviation;

			float dr = sqrtf
			(
				dtl_x * dtl_x +
				dtl_y * dtl_y +
				dtl_z * dtl_z
			);

			dtl_x /= dr;
			dtl_y /= dr;
			dtl_z /= dr;

			// Send shadow rays.

			float shadow = 1.0f;

			for (int j = 0; j < orbs.size(); j++)
			{
				orb orb = orbs[j];

				// Prevent shadow acne.

				float eps = 1e-3f;

				// Send shadow ray.

				float shadow_ray_dummy;

				float shadow_dist = cast
				(
					hit_x + eps * dtl_x,
					hit_y + eps * dtl_y,
					hit_z + eps * dtl_z,

					dtl_x,
					dtl_y,
					dtl_z,

					shadow_ray_dummy,
					shadow_ray_dummy,
					shadow_ray_dummy,

					shadow_ray_dummy,
					shadow_ray_dummy,
					shadow_ray_dummy,

					shadow_ray_dummy
				);

				if (shadow_dist < dr)
				{
					shadow -= 0.2f;
				}
			}

			float l_pow =
			(
				norm_x * dtl_x +
				norm_y * dtl_y +
				norm_z * dtl_z
			);

			l_pow /= dr * dr;

			out_r += fmax(0.0f, hit_orb_r * lamp.r * l_pow * shadow);
			out_g += fmax(0.0f, hit_orb_g * lamp.g * l_pow * shadow);
			out_b += fmax(0.0f, hit_orb_b * lamp.b * l_pow * shadow);

			l_pow = powf(l_pow * dr * dr, 128.0f) / dr / dr * lamp.radius / shadow_rays;

			out_r += fmax(0.0f, lamp.r * l_pow * shadow);
			out_g += fmax(0.0f, lamp.g * l_pow * shadow);
			out_b += fmax(0.0f, lamp.b * l_pow * shadow);
		}
	}

	out_r /= shadow_rays;
	out_g /= shadow_rays;
	out_b /= shadow_rays;

	// Add reflections.

	if (depth > 5)
	{
		return;
	}

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

	float reflection_rays = 1.0f;

	if (depth != 0)
	{
		reflection_rays = 1.0f;
	}

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

	float reflectivity = 0.8f;

	out_r = out_r * (1.0f - reflectivity) + reflect_r * reflectivity;
	out_g = out_g * (1.0f - reflectivity) + reflect_g * reflectivity;
	out_b = out_b * (1.0f - reflectivity) + reflect_b * reflectivity;
}

int main(int argc, char** argv)
{
	srand(2048);

	#include "scene/three.cpp"

	int x_res = 128 * 8;
	int y_res = 128 * 8;

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

			float ray_dx = ((float(i + 0.5f) / float(x_res)) * 2.0f - 1.0f) * fov_adjust * aspect;

			float ray_dy = (1.0f - (float(j + 0.5f) / float(y_res)) * 2.0f) * fov_adjust;

			float ray_dz = -1.0f;

			float ray_length = sqrtf
			(
				ray_dx * ray_dx +
				ray_dy * ray_dy +
				ray_dz * ray_dz
			);

			ray_dx /= ray_length;
			ray_dy /= ray_length;
			ray_dz /= ray_length;

			float color_r;
			float color_g;
			float color_b;

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

			pixel[R] = fmax(0.0f, fmin(255.0f, color_r * 255.0f));
			pixel[G] = fmax(0.0f, fmin(255.0f, color_g * 255.0f));
			pixel[B] = fmax(0.0f, fmin(255.0f, color_b * 255.0f));
		}
	}

	if (!stbi_write_png("orb.png", x_res, y_res, 3, frame_buffer, x_res * 3))
	{
		nuke("Could not save a frame buffer (stbi_write_png).");
	}

	return EXIT_SUCCESS;
}