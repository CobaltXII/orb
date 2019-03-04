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

	orbs.push_back({sinf(r1 * 2.0f * M_PI) * r2 * radius, cosf(r1 * 2.0f * M_PI) * r2 * radius, -32.0f + cosf(r3 * 2.0f * M_PI) * r2 * radius, rand01(), rand01(), rand01(), 3.0f});
}