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

#define STB_PERLIN_IMPLEMENTATION

#include "stb_perlin.h"

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

std::vector<orb> orbs1;
std::vector<orb> orbs2;

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

	for (int i = 0; i < orbs1.size(); i++)
	{
		orb orb1 = orbs1[i];

		float i_lx = orb1.x - ray_ox;
		float i_ly = orb1.y - ray_oy;
		float i_lz = orb1.z - ray_oz;

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
			orb1.radius *
			orb1.radius
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

			hit_orb_x = orb1.x;
			hit_orb_y = orb1.y;
			hit_orb_z = orb1.z;

			hit_orb_r = orb1.r;
			hit_orb_g = orb1.g;
			hit_orb_b = orb1.b;

			hit_orb_radius = orb1.radius;
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
		// Missed everything, use background shader.

		int mode = 1;

		if (mode == -1)
		{
			return;
		}
		else if (mode == 0)
		{
			// Checkerboard.

			float domain_u = 0.1f;
			float domain_v = 0.1f;

			bool check_u = fmod(ray_dx + domain_u * 100.0f, domain_u * 2.0f) >= domain_u;
			bool check_v = fmod(ray_dy + domain_v * 100.0f, domain_v * 2.0f) >= domain_v;

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
			// Sky.

			float frequency = 8.0f;

			float noise = stb_perlin_fbm_noise3
			(
				ray_dx * frequency, 
				ray_dy * frequency, 
				ray_dz * frequency, 

				2.0f, 0.5f, 6
			);

			out_r = (noise + 1.0f) / 2.0f * 0.3f;

			// out_r = 186.0f / 255.0f * 0.64f + (noise + 1.0f) / 2.0f * 0.3f;
			// out_g = 214.0f / 255.0f * 0.64f + (noise + 1.0f) / 2.0f * 0.3f;
			// out_b = 254.0f / 255.0f * 0.64f + (noise + 1.0f) / 2.0f * 0.3f;
		}

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

	// Middle sphere should be textured.

	if (hit_orb_r > 0.5f || hit_orb_g > 0.5f || hit_orb_b > 0.5f)
	{
		float frequency = 0.48f;

		float noise = stb_perlin_ridge_noise3
		(
			hit_x * frequency,
			hit_y * frequency,
			hit_z * frequency,

			2.0f, 0.5f, 1.0f, 6
		);

		if (hit_orb_r > 0.5f)
		{
			hit_orb_r = (noise + 1.0f) / 2.0f;
		}
		else if (hit_orb_g > 0.5f)
		{
			hit_orb_g = (noise + 1.0f) / 2.0f;
		}
		else if (hit_orb_b > 0.5f)
		{
			hit_orb_b = (noise + 1.0f) / 2.0f;
		}
	}

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

	float shadow_rays = 1.0f;

	for (int k = 0; k < int(shadow_rays); k++)
	{
		for (int i = 0; i < orbs2.size(); i++)
		{
			orb orb2 = orbs2[i];

			// Perturb orb for smooth shadows.

			float r1 = rand11();
			float r2 = rand11();
			float r3 = rand11();

			float deviation = shadow_rays - 1.0f;

			// Direction to orb1.

			float dtl_x = (orb2.x + r1 * deviation) - hit_x;
			float dtl_y = (orb2.y + r2 * deviation) - hit_y;
			float dtl_z = (orb2.z + r3 * deviation) - hit_z;

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

			float shadow_ray_dummy;

			float shadow_dist = cast
			(
				hit_x + 1e-3f * dtl_x,
				hit_y + 1e-3f * dtl_y,
				hit_z + 1e-3f * dtl_z,

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

			float shadow = 1.0f;

			if (shadow_dist < dr)
			{
				shadow = 0.0f;
			}

			// Diffuse.

			float l_pow =
			(
				norm_x * dtl_x +
				norm_y * dtl_y +
				norm_z * dtl_z
			);

			l_pow /= dr * dr;

			out_r += fmax(0.0f, hit_orb_r * orb2.r * l_pow * shadow);
			out_g += fmax(0.0f, hit_orb_g * orb2.g * l_pow * shadow);
			out_b += fmax(0.0f, hit_orb_b * orb2.b * l_pow * shadow);

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

			float specular_coefficient = powf(fmax(dot_view_specular, 0.0f), 2048.0f);

			out_r += fmax(0.0f, orb2.r * specular_coefficient * shadow);
			out_g += fmax(0.0f, orb2.g * specular_coefficient * shadow);
			out_b += fmax(0.0f, orb2.b * specular_coefficient * shadow);
		}
	}

	out_r /= shadow_rays;
	out_g /= shadow_rays;
	out_b /= shadow_rays;

	// Add reflections.

	if (depth > 2)
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

	if (hit_orb_radius < 200.0f)
	{
		reflectivity = 0.6f;
	}

	out_r = out_r * (1.0f - reflectivity) + reflect_r * reflectivity;
	out_g = out_g * (1.0f - reflectivity) + reflect_g * reflectivity;
	out_b = out_b * (1.0f - reflectivity) + reflect_b * reflectivity;
}

int main(int argc, char** argv)
{
	orbs1.push_back({0.0f, -3024.0f, 0.0f, 0.0f, 0.0f, 0.0f, 3000.0f});

	orbs1.push_back({0.0f, -8.0f, -56.0f, 0.0f, 1.0f, 0.0f, 16.0f});

	orbs1.push_back({0.0f - 32.0f, -8.0f, -56.0f, 1.0f, 0.0f, 0.0f, 16.0f});
	orbs1.push_back({0.0f + 32.0f, -8.0f, -56.0f, 0.0f, 0.0f, 1.0f, 16.0f});

	orbs2.push_back({25.0f, 50.0f, 0.0f, 1e+4f * 1.4f, 1e+4f * 1.4f, 1e+4f * 1.4f, 50.0f});

	int x_res = 128 * 16;
	int y_res = 128 * 16;

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

	int supersample = 2;

	if (supersample)
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