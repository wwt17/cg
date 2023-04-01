Assignment #5
========================

Use the same way to compile the code as suggested:

```bash
cd Assignment_5
mkdir build
cd build
cmake ..
make
```

Then you will get the compiled executable `assignment5`.

I provided multiple command line arguments for the executable:

* `--filename [filename]`: The filename of the output image. Default to `triangle.png`.
* `--mesh_filename [mesh_filename]`: The filename of the mesh file. Default to `{DATA_DIR}/bunny.off`.
* `--no-refit-mesh`: Without this flag, I will refit the mesh into the $[-1, +1]^3$ cube.
* `-n [n]`: The distance from the camera to the near plane in the camera space. Default to 1.
* `--field_of_view [field_of_view]`:The field of view in degree. Default to 45.
* `--shading [shading]`: The type of shading (one of `silhouette, wireframe, flat, per-vertex`). Default to `wireframe`.
* `--alpha [alpha]` or `-a [alpha]`: Transparency. Default is set by shading.
* `--transformation` or `-t`: With this flag, will output a gif with transformation in Ex.3.
* `--timesteps [timesteps]`: How many timesteps to sample for the gif. Default to 50.
* `--delay [delay]`: Delay of each timestep. Default to 10.
* `--projection_type [projection_type]` or `-p [projection_type]`: The projection type, either `orthographic` (`o`) or `perspective` (`p`). Default to `perspective`.
* `--obj_specular_exponent [obj_specular_exponent]`: The `obj_specular_exponent` in the original code.
* `-w [w]`: Width of the image in pixel. Default to 500.
* `-h [h]`: Height of the image in pixel. Default to 500.

## Ex.1: Load and Render a 3D model

Silhouette of the bunny:

`./assignment5 --shading silhouette -po --filename ../img/bunny_silhouette_orthographic.png`

![](img/bunny_silhouette_orthographic.png?raw=true)

## Ex.2: Shading

Orthographic bunny:

`./assignment5 --shading wireframe -po --filename ../img/bunny_wireframe_orthographic.png`
![](img/bunny_wireframe_orthographic.png?raw=true)

`./assignment5 --shading flat -po --filename ../img/bunny_flat_orthographic.png`
![](img/bunny_flat_orthographic.png?raw=true)

`./assignment5 --shading per-vertex -po --filename ../img/bunny_per-vertex_orthographic.png`
![](img/bunny_per-vertex_orthographic.png?raw=true)

## Ex.3: Object Transformation

Orthographic bunny:

`./assignment5 --shading wireframe -po -t --filename ../img/bunny_wireframe_orthographic.gif`
![](img/bunny_wireframe_orthographic.gif?raw=true)

`./assignment5 --shading flat -po -t --filename ../img/bunny_flat_orthographic.gif`
![](img/bunny_flat_orthographic.gif?raw=true)

`./assignment5 --shading per-vertex -po -t --filename ../img/bunny_per-vertex_orthographic.gif`
![](img/bunny_per-vertex_orthographic.gif?raw=true)

## Ex.4: Camera

Perspective bunny:

`./assignment5 --shading wireframe --filename ../img/bunny_wireframe.png`
![](img/bunny_wireframe.png?raw=true)

`./assignment5 --shading flat --filename ../img/bunny_flat.png`
![](img/bunny_flat.png?raw=true)

`./assignment5 --shading per-vertex --filename ../img/bunny_per-vertex.png`
![](img/bunny_per-vertex.png?raw=true)

`./assignment5 --shading wireframe -t --filename ../img/bunny_wireframe.gif`
![](img/bunny_wireframe.gif?raw=true)

`./assignment5 --shading flat -t --filename ../img/bunny_flat.gif`
![](img/bunny_flat.gif?raw=true)

`./assignment5 --shading per-vertex -t --filename ../img/bunny_per-vertex.gif`
![](img/bunny_per-vertex.gif?raw=true)

Adapting the aspect ratio:

`./assignment5 --shading per-vertex -w 1000 -h 500 --filename ../img/bunny_per-vertex_wide.png`

![](img/bunny_per-vertex_wide.png?raw=true)

`./assignment5 --shading per-vertex -w 500 -h 1000 --field_of_view 80 --filename ../img/bunny_per-vertex_tall.png` (I use a larger `field_of_view` so that the whole object can be seen)

![](img/bunny_per-vertex_tall.png?raw=true)

Other objects:

bumpy cube:

`./assignment5 --mesh_filename ../data/bumpy_cube.off --shading per-vertex --filename ../img/bumpy_cube.png`

![](img/bumpy_cube.png?raw=true)
