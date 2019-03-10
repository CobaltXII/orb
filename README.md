# Orb

![Alt text](https://github.com/CobaltXII/orb/blob/master/img/orb_soft_shadow.png?raw=true)

## Summary

Orb started out as a weekend raytracing project, but got a little bigger. It currently supports scene files, diffuse lighting, specular lighting, recursion, reflection, refraction, procedural environments, procedural texturing, supersampling, and gamma correction. Orb was originally named Orb because it could only raytrace spheres, but it can now raytrace ellipsoids, planes, capsules, cylinders and cones as well.

## Libraries

I make heavy use of Sean T. Barrett's libraries. They add a significant amount to compilation time, so I compile them once and then link them for evermore. To compile them, use the following commands.

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

I also use Ben Hoyt's INI parsing library for loading portable scene files. The library adds a small delay to compilation time, so I compile it once and then link it. To compile it, use the following command.

### inih

```bash
clang -c inih/ini.c -o inih/ini.o -Ofast
```

I also use Syoyo Fujita's Wavefront OBJ parsing library for loading portable model files. The library adds a large delay to compilation time, so I compile it once and then link it. To compile it, use the following command.

### tinyobjloader

```bash
clang -c tinyobjloader/tiny_obj_loader.cc -o tinyobjloader/tiny_obj_loader.o -Ofast
```

## Compiling

Just link the libraries and compile as default. I was too lazy to write a Makefile.

```bash
./orb.sh
```

## Scenes

There are a bunch of scene files included in this repository that you can look at to see how they work. Every primitive is mentioned at least once, so it should be easy to tweak a few parameters and see what they do. There is no formal documentation, so you will have to take a peek at the code if you get stumped.

## Credits

Thanks to Inigo Quilez for his great articles and intersector functions. Thanks to Cyrille Favreau for his ellipsoid intersection function. Thanks to Brook Heisler for a great tutorial on getting started with raytracing. Thanks to Ben Hoyt for his great INI library. Without you guys, this project would never have happened, so thanks again!