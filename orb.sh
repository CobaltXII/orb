clang++ orb.cpp inih/ini.o lib/stb_image.o lib/stb_image_write.o lib/stb_perlin.o -o orb -std=c++11 -Ofast -march=native && ./orb ini/spheres.ini