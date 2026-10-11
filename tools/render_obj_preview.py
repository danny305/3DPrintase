"""Headless Blender preview for comparing printable OBJ proportions."""

import math
import sys
from pathlib import Path

import bpy
from mathutils import Vector


def args_after_separator():
    return sys.argv[sys.argv.index("--") + 1 :] if "--" in sys.argv else []


source, destination = map(Path, args_after_separator()[:2])
bpy.ops.wm.read_factory_settings(use_empty=True)
bpy.ops.wm.obj_import(filepath=str(source))
objects = [obj for obj in bpy.context.scene.objects if obj.type == "MESH"]

# Use the reference print's single teal material so only geometry is compared.
material = bpy.data.materials.new("reference_teal")
material.diffuse_color = (0.02, 0.62, 0.48, 1.0)
material.metallic = 0.0
material.roughness = 0.30
for obj in objects:
    obj.data.materials.clear()
    obj.data.materials.append(material)

world_points = [obj.matrix_world @ Vector(corner) for obj in objects for corner in obj.bound_box]
minimum = Vector(tuple(min(point[i] for point in world_points) for i in range(3)))
maximum = Vector(tuple(max(point[i] for point in world_points) for i in range(3)))
center = (minimum + maximum) / 2
extent = max(maximum - minimum)

# Ground plane and neutral studio lighting approximate the photographed model.
bpy.ops.mesh.primitive_plane_add(size=extent * 4, location=(center.x, center.y, minimum.z - 0.03 * extent))
plane = bpy.context.object
ground = bpy.data.materials.new("ground")
ground.diffuse_color = (0.08, 0.08, 0.08, 1.0)
ground.roughness = 0.7
plane.data.materials.append(ground)
plane.hide_render = True

for location, energy, size in [
    ((center.x - extent, center.y - extent, center.z + 2 * extent), 900, extent),
    ((center.x + 2 * extent, center.y + extent, center.z + extent), 500, extent),
]:
    bpy.ops.object.light_add(type="AREA", location=location)
    light = bpy.context.object
    light.data.energy = energy
    light.data.shape = "DISK"
    light.data.size = size
    direction = center - light.location
    light.rotation_euler = direction.to_track_quat("-Z", "Y").to_euler()

bpy.ops.object.camera_add(location=(center.x + 1.35 * extent, center.y - 1.55 * extent, center.z + 1.05 * extent))
camera = bpy.context.object
camera.rotation_euler = (center - camera.location).to_track_quat("-Z", "Y").to_euler()
camera.data.type = "ORTHO"
camera.data.ortho_scale = extent * 1.42
bpy.context.scene.camera = camera

scene = bpy.context.scene
scene.render.engine = "BLENDER_WORKBENCH"
scene.display.shading.light = "STUDIO"
scene.display.shading.color_type = "MATERIAL"
scene.display.shading.show_shadows = True
scene.display.shading.show_cavity = True
scene.display.shading.cavity_type = "WORLD"
scene.render.resolution_x = 900
scene.render.resolution_y = 900
scene.render.resolution_percentage = 100
scene.render.image_settings.file_format = "PNG"
scene.render.filepath = str(destination)
scene.render.film_transparent = False
if scene.world is None:
    scene.world = bpy.data.worlds.new("preview_world")
scene.world.color = (0.08, 0.08, 0.08)
scene.view_settings.look = "AgX - Medium High Contrast"
bpy.ops.render.render(write_still=True)
