Assignment #4
========================

Use the same way to compile the code as suggested:

```bash
cd Assignment_4
mkdir build
cd build
cmake ..
make
```

Then you will get the compiled executable `assignment4`.

I provided multiple command line arguments for the executable:

* `--filename [filename]`: The filename of the output image. Default to `raytrace.png`.
* `--mesh_filename [mesh_filename]`: The filename of the mesh file. Default to `{DATA_DIR}/dodeca.off`.
* `--brute_force`: If set, use brute force without BVH. This is used to check if the result using the BVH is the same as brute force.
* `--centroid_sorting_split`: If set, in top-down construction of the tree, split the set of triangles by sorting them by their centroids along the direction with maximum span of the bounding box of the centroids. Default to unset because this is too slow for large number of meshes like `dragon.off`.
* `--construction [construction_type]`: The type of construction, either `top-down` (default) or `bottom-up`.
* `--criterion [criterion]`: The criterion (cost function) used to decide which pair of node (with the least value on this function) should be merged next. Can be `centroid_distance` (default) or `volume_increase`.
* `--focal_length [focal_length]` or `-f [focal_length]`: The focal length. Default to 2.
* `--field_of_view [field_of_view]`:The field of view in degree. Default to 45.
* `--projection_type [projection_type]` or `-p [projection_type]`: The projection type, either `orthographic` (`o`) or `perspective` (`p`). Default to `perspective`.
* `--obj_specular_exponent [obj_specular_exponent]`: The `obj_specular_exponent` in the original code.
* `--grid_size [grid_size]` or `-g [grid_size]`: The grid size in the original code used in Perlin noise. Default to 20.
* `-w [w]`: Width of the image in pixel. Default to 640.
* `-h [h]`: Height of the image in pixel. Default to 480.

## Ex.1: Triangle Mesh Ray Tracing

Output with orthographic camera (by running `./assignment4 --brute_force --projection_type orthographic --filename ../img/dodeca_bf_orthographic.png`):
![](img/dodeca_bf_orthographic.png?raw=true)

Output with perspective camera (by running `./assignment4 --brute_force --filename ../img/dodeca_bf.png`):
![](img/dodeca_bf.png?raw=true)

## Ex.2: AABB Trees

### Top-Down Construction

Note that splitting the set of triangles by sorting them by their centroids along the direction with maximum span of the bounding box of the centroids (`--centroid_sorting_split`) is very slow on `dragon.off`, so I do not use this in default.

Note I added a little speed-up optimization when finding the nearest object on the tree: I terminate the recursion on the subtree if the distance to the bbox of this subtree is farther than the already-known nearest object.

Output of dodeca: (by running `./assignment4 --filename ../img/dodeca_top-down.png`):
![](img/dodeca_top-down.png?raw=true)

Output of cube: (by running `./assignment4 --mesh_filename ../data/cube.off --filename ../img/cube_top-down.png`):
![](img/cube_top-down.png?raw=true)

Output of bunny: (by running `./assignment4 --mesh_filename ../data/bunny.off --filename ../img/bunny_top-down.png`):
![](img/bunny_top-down.png?raw=true)

Output of dragon: (by running `./assignment4 --mesh_filename ../data/dragon.off --filename ../img/dragon_top-down.png`):

![](img/dragon_top-down.png?raw=true)

### Bottom-Up Construction

Since the time complexity of our simple implementation is $O(n^3)$, it takes prohibitively long to build the tree for the dragon. So I just show examples on the bunny. Results on other meshes are the same.

Output of bunny using `centroid_distance` criterion (by running `./assignment4 --mesh_filename ../data/bunny.off --construction bottom-up --criterion centroid_distance --filename ../img/bunny_bottom-up_centroid_distance.png`):

![](img/bunny_bottom-up_centroid_distance.png?raw=true)

Output of bunny using `centroid_distance` criterion (by running `./assignment4 --mesh_filename ../data/bunny.off --construction bottom-up --criterion volume_increase --filename ../img/bunny_bottom-up_volume_increase.png`):

![](img/bunny_bottom-up_volume_increase.png?raw=true)
