for (float i = -32.0f; i <= 32.0f; i += 8.0f)
for (float j = -32.0f; j <= 32.0f; j += 8.0f)
{
	orbs.push_back({i, j, -48.0f, (i + 32.0f) / 128.0f, (j + 32.0f) / 128.0f, 0.16f, 3.84f});
}

lamps.push_back({0.0f, 0.0f, 10004.0f, 1e+8f * 2.4f, 1e+8f * 2.4f, 1e+8f * 2.4f, 10000.0f});