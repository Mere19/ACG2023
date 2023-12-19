import mitsuba as mi

fname = "cbox_disney_metallic_mitsuba"

mi.set_variant('scalar_rgb')
scene = mi.load_file(fname + ".xml")
img = mi.render(scene)
bmp_small = mi.Bitmap(img).convert(mi.Bitmap.PixelFormat.RGB, mi.Struct.Type.UInt8, True)
bmp_small.write(fname + ".png")
mi.Bitmap(img).write(fname + ".exr")