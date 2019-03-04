# Orb

![Alt text](https://github.com/CobaltXII/orb/blob/master/img/orb_2.png?raw=true)

A tiny little raytracer that I made in a weekend. Does not support recursion, reflection, refraction or any of those cool effects. Hell, it doesn't even support proper specular lighting. However, it still makes cool images if you use it properly.

## Compiling

I use this script

```bash
clang++ orb.cpp -o orb -std=c++11 -Ofast && ./orb
```