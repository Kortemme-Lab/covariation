import os
import subprocess
import time
import colorsys
from fsio import write_temp_file, read_file


def run_r_script_contents(r_script):
    r_temp_filename = write_temp_file(".", r_script)
    try:
        run_r_script(r_temp_filename)
        os.remove(r_temp_filename)
    except:
        os.remove(r_temp_filename)
        raise


def run_r_script(r_script_filename, cwd = '.'):
    p = subprocess.Popen(["R", "CMD", "BATCH", r_script_filename], cwd = cwd)
    while True:
        time.sleep(0.3)
        errcode = p.poll()
        if errcode != None:
            break
    rout = "{0}out".format(r_script_filename)
    rout_contents = None
    if os.path.exists(rout):
        rout_contents = read_file(rout)
        os.remove(rout)
    rdata_file = os.path.join(os.path.split(r_script_filename)[0], '.RData')
    if os.path.exists(rdata_file):
        os.remove(rdata_file)
    if errcode != 0:
        print(rout_contents)
        raise Exception("The R script failed with error code %d." % errcode)
    return rout_contents



#Change saturation of hexadecimal color
#input: hexadecimal color, saturation adjustment
#output: saturation-adjusted hexadecimal color string
def saturate_hex_color(hexcolor, adjustment = 1.0):
    assert(adjustment >= 0 and len(hexcolor) >= 1)
    prefix = ""
    if hexcolor[0] == '#':
        hexcolor = hexcolor[1:]
        prefix = "#"
    #-
    assert(len(hexcolor) == 6)
    if adjustment == 1.0:
        return "%s%s" % (prefix, hexcolor)
    else:
        hsvColor = list(colorsys.rgb_to_hsv(int(hexcolor[0:2], 16)/255.0, int(hexcolor[2:4], 16)/255.0, int(hexcolor[4:6], 16)/255.0))
        hsvColor[1] = min(1.0, hsvColor[1] * adjustment)
        rgbColor = [min(255, 255 * v) for v in colorsys.hsv_to_rgb(hsvColor[0], hsvColor[1], hsvColor[2])]
        return "%s%.2x%.2x%.2x" % (prefix, rgbColor[0], rgbColor[1], rgbColor[2])


#Get list of hexadecimal colors
#input: number of colors needed (int), start value (optional), saturation adjustment (float; optional)
#output: list of hexadecimal color strings
# Added by Shane
def color_wheel(n, start = 15, saturation_adjustment = None):
    hues = range(start, start + 360, 360/n)
    rgbcolors = ['%x%x%x' % (255 * hlscol[0], 255 * hlscol[1], 255 * hlscol[2]) for hlscol in [colorsys.hls_to_rgb(float(h % 360) / 360.0, 0.65, 1.00) for h in hues]]
    if saturation_adjustment:
        return [saturate_hex_color(rgbcol, saturation_adjustment) for rgbcol in rgbcolors]
    else:
        return rgbcolors
