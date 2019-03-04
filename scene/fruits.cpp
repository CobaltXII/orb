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