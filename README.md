# Recursive Ray Tracer

**Jake Mingolla**

**April 2015**

**Collaborators: [Casey Gowrie](https://github.com/ctgowrie)**

**Category: Computer Graphics**

**Language(s): C++**

### About

This project implements a recursive ray tracer for rendering 3D scenes composed of primitive objects.

For each pixel on the near panel of the viewing frustum, a ray is cast out into the scene to calculate the color intensity of the pixel. Each object within the scene is tested for intersection with each ray - if it intersects, a color intensity calculation is performed using the [Phong Lighting Model](https://en.wikipedia.org/wiki/Phong_reflection_model) to generate ambient, diffuse, and specular contributions from each of the lights in the scene. This process continues with the reflected ray to accumulate contributions from recursive rays until the recursion limit is reached. After 6-7 recursive calls the contributions become negligible.

These calculations are handled entirely by the CPU and only use OpenGL functionality to set pixel colors. While this is orders of magnitude slower than using a GPU for scene rendering, understanding ray tracing is an important step in graphics development.

Scenes can be loaded in XML format and can be downloaded from the Comp 175 website [here](http://www.cs.tufts.edu/comp/175/assignments/a5/data-ray.tar.gz). The primitives supported by the project are Cone, Cube, Cylinder, Hourglass (Special), and Sphere.

The scene also supports "isect-only" mode to speed up calculations to ignore color intensity and only show white for intersection and black for rays that do not intersect with any objects in the scene.

Finally, this project also includes texture mapping across all of the primitive shapes. Any texture in a .ppm file format can be "stretched" across a primitive shape such as a cube, sphere, cone, etc.

### Screenshots

![Image](http://i.imgur.com/CAsZRof.jpg)

Due to some compression artifacts a little bit of the clarity is lost in this screenshot. However, you can clearly see the recursive aspect of the ray tracer in the reflections within the spheres. In addition, this screenshot shows the texture mapping capabilities as a low quality image of the earth is mapped over a cube.


### Note

In order to simplify the submission process in Comp 175 - Computer Graphics, all files are in the same root directory. In an ideal system all files should be separated into different folders based on type.

In addition, the code has relatively sparse comments since a working submission was heavily prioritized over documentation. Since most of the code is relatively straight-forward I don't see this being a readability issue - in future versions I will continue to add documentation and comments.

### Dependencies
- OpenGL > 1.x
  - Standard graphics library
  - Any version past 1.x has support for primitive object loading through the standard OpenGL TRIANGLE_STRIP pipeline.
- GLUI
  - Provides user interface support for OpenGL
- [tinyxml](www.sourceforge.net/projects/tinyxml) parser

### Collaborators

- **Mike Shah** for providing support code necessary for loading .ppm images.

- **Remco Chang** for introducing this project as a part of Comp 175 - Computer Graphics at Tufts Spring 2015. In addition, Remco provided support code listed in Algebra.h as well as all of the modules necessary for parsing the scene files.
