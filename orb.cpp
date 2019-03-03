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

int main(int argc, char** argv)
{
	srand(2048);

	std::vector<orb> orbs;

	std::vector<lamp> lamps;

	orbs.push_back({0.0f, 0.0f, -3064.0f, 0.5f, 0.75f, 0.25f, 3000.0f});

	for (int i = 0; i < 128; i++)
	{
		float r1 = rand01();
		float r2 = rand01();
		float r3 = rand01();

		float radius = 32.0f;

		if (true)
		{
			r2 = sqrtf(r2);
		}

		orbs.push_back
		(
			{
				sinf(r1 * 2.0f * M_PI) * r2 * radius,
				cosf(r1 * 2.0f * M_PI) * r2 * radius,

				-24.0f + cosf(r3 * 2.0f * M_PI) * r2 * radius,

				(sinf(r1 * 2.0f * M_PI) * r2 * radius + radius) / radius,
				(cosf(r1 * 2.0f * M_PI) * r2 * radius + radius) / radius,

				(-32.0f + cosf(r3 * 2.0f * M_PI) * r2 * radius + radius) / radius,

				3.0f
			}
		);
	}

	lamps.push_back({0.0f, 0.0f, 104.0f, 1e+4f * 1.4f, 1e+4f * 1.4f, 1e+4f * 1.4f, 50.0f});

	int x_res = 128 * 32;
	int y_res = 128 * 32;

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

			float ray_x = ((float(i + 0.5f) / float(x_res)) * 2.0f - 1.0f) * fov_adjust * aspect;

			float ray_y = (1.0f - (float(j + 0.5f) / float(y_res)) * 2.0f) * fov_adjust;

			float ray_z = -1.0f;

			float ray_length = sqrtf
			(
				ray_x * ray_x +
				ray_y * ray_y +
				ray_z * ray_z
			);

			ray_x /= ray_length;
			ray_y /= ray_length;
			ray_z /= ray_length;

			// Find nearest geometry intersection.

			float min_dist = std::numeric_limits<float>::max();

			float min_ox;
			float min_oy;
			float min_oz;

			float min_or;
			float min_og;
			float min_ob;

			float min_radius;

			for (int i = 0; i < orbs.size(); i++)
			{
				orb orb = orbs[i];

				float i_lx = orb.x - ray_ox;
				float i_ly = orb.y - ray_oy;
				float i_lz = orb.z - ray_oz;

				float i_adj2 =
				(
					i_lx * ray_x +
					i_ly * ray_y +
					i_lz * ray_z
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

					min_ox = orb.x;
					min_oy = orb.y;
					min_oz = orb.z;

					min_or = orb.r;
					min_og = orb.g;
					min_ob = orb.b;

					min_radius = orb.radius;
				}
			}

			// Generate color.

			if (min_dist != std::numeric_limits<float>::max())
			{
				float r = 0.0f;
				float g = 0.0f;
				float b = 0.0f;

				float hit_x = ray_ox + ray_x * min_dist;
				float hit_y = ray_oy + ray_y * min_dist;
				float hit_z = ray_oz + ray_z * min_dist;

				if (min_radius > 200.0f)
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
						min_or = 0.15f * 2.0f;
						min_og = 0.15f * 2.0f;
						min_ob = 0.15f * 2.0f;
					}
					else
					{
						min_or = 0.85f * 2.0f;
						min_og = 0.85f * 2.0f;
						min_ob = 0.85f * 2.0f;
					}
				}

				float norm_x = hit_x - min_ox;
				float norm_y = hit_y - min_oy;
				float norm_z = hit_z - min_oz;

				float norm_length = sqrtf
				(
					norm_x * norm_x +
					norm_y * norm_y +
					norm_z * norm_z
				);

				norm_x /= norm_length;
				norm_y /= norm_length;
				norm_z /= norm_length;

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

						bool shadow = false;

						for (int j = 0; j < orbs.size(); j++)
						{
							orb orb = orbs[j];

							// Prevent shadow acne.

							float eps = 1e-3f;

							float i_lx = orb.x - hit_x - dtl_x * eps;
							float i_ly = orb.y - hit_y - dtl_y * eps;
							float i_lz = orb.z - hit_z - dtl_z * eps;

							float i_adj2 =
							(
								i_lx * dtl_x +
								i_ly * dtl_y +
								i_lz * dtl_z
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

							float i_d = fmin(i_t0, i_t1);

							if (i_d < dr)
							{
								shadow = true;
							}
						}

						if (shadow)
						{
							continue;
						}

						float l_pow =
						(
							norm_x * dtl_x +
							norm_y * dtl_y +
							norm_z * dtl_z
						);

						l_pow /= dr * dr;

						r += fmax(0.0f, min_or * lamp.r * l_pow);
						g += fmax(0.0f, min_og * lamp.g * l_pow);
						b += fmax(0.0f, min_ob * lamp.b * l_pow);

						// Fake specularity.

						l_pow = powf(l_pow * dr * dr, 256.0f) / dr / dr * shadow_rays;

						r += fmax(0.0f, min_or * lamp.r * l_pow);
						g += fmax(0.0f, min_og * lamp.g * l_pow);
						b += fmax(0.0f, min_ob * lamp.b * l_pow);
					}
				}

				r /= shadow_rays;
				g /= shadow_rays;
				b /= shadow_rays;

				pixel[R] = fmax(0.0f, fmin(255.0f, r * 255.0f));
				pixel[G] = fmax(0.0f, fmin(255.0f, g * 255.0f));
				pixel[B] = fmax(0.0f, fmin(255.0f, b * 255.0f));
			}
		}
	}

	if (!stbi_write_png("orb.png", x_res, y_res, 3, frame_buffer, x_res * 3))
	{
		nuke("Could not save a frame buffer (stbi_write_png).");
	}

	return EXIT_SUCCESS;
}