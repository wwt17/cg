Assignment #6
========================

**The video recording is [illustration.mkv](illustration.mkv). Watch the video for illustration of all operations.**

Use the same way to compile the code as suggested:

```bash
cd Assignment_6
mkdir build
cd build
cmake ..
make
```

Then you will get the compiled executable `RasterViewer`.

I provided multiple command line arguments for the executable:

* `--redraw_interval/-i [redraw_interval]`: The `redraw_interval` passed to the SDL viewer. Default to 30.
* `--width/-w [width]`: Width of the viewer in pixel. Default to 500.
* `--height/-h [height]`: Height of the viewer in pixel. Default to 500.
* `--line_thickness/-t [line_thickness]`: the thickness of the line when drawing triangles. Default to 1.
* `--n_inter_frames/-f [n_inter_frames]`: interpolating two adjacent keyframes by splitting the time span into `n_inter_frames` segments with equal span, so we have `n_inter_frames - 1` frames in-between. Default to 5.

Ex.1: Triangle Soup Editor
--------------------------

I implemented the insertion (enabled by the key <kbd>i</kbd>), translation (enabled by the key <kbd>o</kbd>), and delete (enabled by the key <kbd>p</kbd>) mode as described. Watch the video for the results.


Ex.2: Rotation/Scale
--------------------

I implemented the translations when pressing the keys <kbd>h</kbd>, <kbd>j</kbd>, <kbd>k</kbd>, <kbd>l</kbd> as described. Watch the video for the results.

Ex.3: Colors
------------

I implemented the color mode (enabled by the key <kbd>c</kbd>) as described. Colors are picked by a key from <kbd>1</kbd> to <kbd>9</kbd>. Watch the video for the results.


Ex.4: Shader Translation/Scaling/Rotation
----------------------------------------------------

I did all spatial transformations inside the vertex shader. I maintained a transform matrix for each triangle. In the vertex shader, I first look up the triangle of the vertex, and then apply the transform matrix. See the code for details.

Ex.5: View Control
------------------

I implemented the view control as described. Note that for unknown reason my SDL can only recognize the key <kbd>+</kbd>, which is obtained by the combination <kbd>Shift + =</kbd> as <kbd>Shift</kbd> + <kbd>=</kbd> separately. So I actually used <kbd>=</kbd> for the zoom-in operation. Watch the video for the results.


Ex.6: Add Keyframing
--------------------

I implemented the keyframing and interpolation. Assume the initial state is the first keyframe. Each time I present a (interpolated) frame. By pressing the keys described below, you can navigate the frames by advancing the next frame. Assume we have `n_keyframes` keyframes now. There will be `n_inter_frames - 1` interpolated frames between each two adjacent keyframes. (So there are `(n_keyframes - 1) * n_inter_frames + 1` frames in total.) We refer to the keyframe before the current interpolated frame as the current keyframe. When trying to advance to the next frame at the last frame, it will instead advance to the first frame. You are able to do all the operations in all modes described above on the current keyframe. (Note that translations on interpolated frames are undefined.)

Here is the set of keys used:

* <kbd>;</kbd>: insert a new keyframe after the current keyframe. The new keyframe is initialized to be the same as the current keyframe.
* <kbd>'</kbd>: insert a new keyframe after the last keyframe.
* <kbd>/</kbd>: delete the current keyframe.
* <kbd>\\</kbd>: clear keyframes (except the first one).
* <kbd>Space</kbd>: advance to the next frame by linear interpolation.
* <kbd>Tab</kbd>: advance to the next frame by Bézier interpolation.

Watch the video for the illustration of these operations.
