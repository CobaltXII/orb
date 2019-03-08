# Orb

![Alt text](https://github.com/CobaltXII/orb/blob/master/img/orb_capped_cones.png?raw=true)

## Summary

Orb started out as a weekend raytracing project, but got a little bigger. It currently supports diffuse lighting, specular lighting, recursion, reflection, refraction, procedural environments, procedural texturing, supersampling, and gamma correction. Orb was named Orb originally because it could only raytrace spheres, but it can now raytrace ellipsoids, planes and cones as well.

## Libraries

I make heavy Sean T. Barrett's libraries. They add a significant significant amount to compilation time, so I compile them once and then link them for evermore. To compile them, use the following commands.

### stb_image

```bash
clang -D STB_IMAGE_IMPLEMENTATION -c stb/stb_image.c -o lib/stb_image.o -Ofast
```

### stb_image_write

```bash
clang -D STB_IMAGE_WRITE_IMPLEMENTATION -c stb/stb_image_write.c -o lib/stb_image_write.o -Ofast
```

### stb_perlin

```bash
clang -D STB_PERLIN_IMPLEMENTATION -c stb/stb_perlin.c -o lib/stb_perlin.o -Ofast
```

## Compiling

Just link the libraries and compile as default.

```bash
clang++ orb.cpp lib/stb_image.o lib/stb_image_write.o lib/stb_perlin.o -o orb -std=c++11 -Ofast && ./orb
```

## Credits

Thanks to Inigo Quilez for his great articles and intersector functions. Thanks to Cyrille Favreau for his ellipsoid intersection function. Thanks to Brook Heisler for a great tutorial on getting started with raytracing. Without you guys, this project would never have happened, so thanks!