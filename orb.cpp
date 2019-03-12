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

#include "tinyobjloader/tiny_obj_loader.h"

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

#include "hpp/shape.hpp"

#include "hpp/intersection.hpp"

#include "hpp/ellipsoid.hpp"
#include "hpp/cylinder.hpp"
#include "hpp/triangle.hpp"
#include "hpp/capsule.hpp"
#include "hpp/sphere.hpp"
#include "hpp/plane.hpp"
#include "hpp/cone.hpp"

std::vector<shape*> shapes;

#include "hpp/estimator.hpp"
#include "hpp/intersect.hpp"

#include "hpp/light.hpp"

std::vector<light> lights;

#include "hpp/cast.hpp"

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
	float hit_shape_material1;
	float hit_shape_material2;
	float hit_shape_material3;
	float hit_shape_material4;
	float hit_shape_material5;

	float hit_shape_r;
	float hit_shape_g;
	float hit_shape_b;

	float norm_x;
	float norm_y;
	float norm_z;

	float texture_u;
	float texture_v;

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
		&hit_shape_material5,

		&hit_shape_r,
		&hit_shape_g,
		&hit_shape_b,

		&norm_x,
		&norm_y,
		&norm_z,

		&texture_u,
		&texture_v,

		&hit_shape
	);

	out_r = 0.0f;
	out_g = 0.0f;
	out_b = 0.0f;

	if (min_dist == std::numeric_limits<float>::max())
	{
		// Missed everything, use background shader.

		float frequency = 8.0f;

		float noise = stb_perlin_fbm_noise3
		(
			ray_dx * frequency, 
			ray_dy * frequency, 
			ray_dz * frequency, 

			2.0f, 0.5f, 6
		);

		out_r = fmax(0.0f, fmin(1.0f, (noise + 1.0f) / 2.0f * 0.032f));

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

	// Check for planes, which should be checkered. Everything else gets a
	// procedural noise texture.

	if (hit_shape->primitive == shape_type::st_plane)
	{
		float domain_u = 10.0f;
		float domain_v = 10.0f;

		float expand = 16384.0f;

		int check_u = fmod(texture_u + domain_u * expand, domain_u * 2.0f) >= domain_u;
		int check_v = fmod(texture_v + domain_v * expand, domain_v * 2.0f) >= domain_v;

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
		float frequency = 0.48f;

		float noise = stb_perlin_ridge_noise3
		(
			hit_x * frequency,
			hit_y * frequency,
			hit_z * frequency,

			2.0f, 0.5f, 1.0f, 6
		);

		hit_shape_r = hit_shape_r * ((noise + 1.0f) / 2.0f);
		hit_shape_g = hit_shape_g * ((noise + 1.0f) / 2.0f);
		hit_shape_b = hit_shape_b * ((noise + 1.0f) / 2.0f);
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

		float shadow = shadow_ray
		(
			dr,

			hit_x + 1e-1f * dtl_x,
			hit_y + 1e-1f * dtl_y,
			hit_z + 1e-1f * dtl_z,

			dtl_x,
			dtl_y,
			dtl_z,

			light1.radius
		);

		// Diffuse.

		float l_pow =
		(
			norm_x * dtl_x +
			norm_y * dtl_y +
			norm_z * dtl_z
		);

		l_pow /= dr * dr;

		out_r += fmax(0.0f, hit_shape_r * light1.r * l_pow * shadow);
		out_g += fmax(0.0f, hit_shape_g * light1.g * l_pow * shadow);
		out_b += fmax(0.0f, hit_shape_b * light1.b * l_pow * shadow);

		// Specular.

		float specular_constant = hit_shape_material5;

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

		float specular_coefficient = powf(fmax(dot_view_specular, 0.0f), hit_shape_material4) / dr / dr;

		out_r += fmax(0.0f, light1.r * specular_coefficient * specular_constant * shadow);
		out_g += fmax(0.0f, light1.g * specular_coefficient * specular_constant * shadow);
		out_b += fmax(0.0f, light1.b * specular_coefficient * specular_constant * shadow);
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

	#ifdef CLAMP_TRACE

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
						inif(section.first, "material4"),
						inif(section.first, "material5")
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
						inif(section.first, "material4"),
						inif(section.first, "material5")
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
						inif(section.first, "material4"),
						inif(section.first, "material5")
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
						inif(section.first, "material4"),
						inif(section.first, "material5")
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
						inif(section.first, "material4"),
						inif(section.first, "material5")
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
						inif(section.first, "material4"),
						inif(section.first, "material5")
					)
				);
			}
			else if (type == "triangle")
			{
				shapes.push_back
				(
					new triangle
					(
						inif(section.first, "x0"),
						inif(section.first, "y0"),
						inif(section.first, "z0"),

						inif(section.first, "x1"),
						inif(section.first, "y1"),
						inif(section.first, "z1"),

						inif(section.first, "x2"),
						inif(section.first, "y2"),
						inif(section.first, "z2"),

						inif(section.first, "r"),
						inif(section.first, "g"),
						inif(section.first, "b"),

						inif(section.first, "material1"),
						inif(section.first, "material2"),
						inif(section.first, "material3"),
						inif(section.first, "material4"),
						inif(section.first, "material5")
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
			else if (type == "model")
			{
				std::string path = inis(section.first, "path");

				float scale_x = inif(section.first, "scale_x");
				float scale_y = inif(section.first, "scale_y");
				float scale_z = inif(section.first, "scale_z");

				float off_x = inif(section.first, "off_x");
				float off_y = inif(section.first, "off_y");
				float off_z = inif(section.first, "off_z");

				float r = inif(section.first, "r");
				float g = inif(section.first, "g");
				float b = inif(section.first, "b");

				float material1 = inif(section.first, "material1");
				float material2 = inif(section.first, "material2");
				float material3 = inif(section.first, "material3");
				float material4 = inif(section.first, "material4");
				float material5 = inif(section.first, "material5");

				tinyobj::attrib_t obj_attrib;

				std::vector<tinyobj::shape_t> obj_shapes;

				std::vector<tinyobj::material_t> obj_materials;

				std::string obj_warning;

				std::string obj_error;

				if (!tinyobj::LoadObj(&obj_attrib, &obj_shapes, &obj_materials, &obj_warning, &obj_error, path.c_str()))
				{
					if (!obj_error.empty())
					{
						std::cout << obj_error << std::endl;
					}

					return EXIT_FAILURE;
				}

				if (!obj_warning.empty())
				{
					std::cout << obj_warning << std::endl;
				}

				for (size_t s = 0; s < obj_shapes.size(); s++)
				{
					size_t index_offset = 0;

					for (size_t f = 0; f < obj_shapes[s].mesh.num_face_vertices.size(); f++)
					{
						int fv = obj_shapes[s].mesh.num_face_vertices[f];

						if (fv == 3)
						{
							tinyobj::real_t avx[3];
							tinyobj::real_t avy[3];
							tinyobj::real_t avz[3];

							for (size_t v = 0; v < fv; v++)
							{
								tinyobj::index_t idx = obj_shapes[s].mesh.indices[index_offset + v];

								avx[v] = obj_attrib.vertices[3 * idx.vertex_index + 0];
								avy[v] = obj_attrib.vertices[3 * idx.vertex_index + 1];
								avz[v] = obj_attrib.vertices[3 * idx.vertex_index + 2];
							}

							for (int i = 0; i < 3; i++)
							{
								avx[i] *= scale_x;
								avy[i] *= scale_y;
								avz[i] *= scale_z;

								avx[i] += off_x;
								avy[i] += off_y;
								avz[i] += off_z;
							}

							shapes.push_back
							(
								new triangle
								(
									avx[0], avy[0], avz[0],
									avx[1], avy[1], avz[1],
									avx[2], avy[2], avz[2],

									r, g, b,

									material1,
									material2,
									material3,
									material4,
									material5
								)
							);
						}

						index_offset += fv;
			    	}
				}
			}
			else
			{
				nuke("Unknown type name (main).");
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