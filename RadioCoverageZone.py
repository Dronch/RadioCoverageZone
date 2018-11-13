#!/usr/bin/env python
"""
The MIT License
Copyright (c) 2018 Zybin Andrey
Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:
The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.
THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.
"""

from PIL import Image, TiffImagePlugin
import math


def calc_coverage_area(input_filename, output_filename, elevation, wavelength, cutoff=2.7):
    """Calc coverage area of radio transmitter with Fresnel zone method

    Args:
        input_filename: filepath of geotiff with elevation matrix of target area
        output_filename: filepath of output geotiff with signal reception zone
		elevation: radio transmitter elevation above ground
		wavelength: transmitter wavelength
		cutoff: intensity cutoff
    """
	
    TiffImagePlugin.WRITE_LIBTIFF = False

    image_src = Image.open(input_filename)
    w = image_src.size[0]
    h = image_src.size[1]
    pix = image_src.load()

    col = 0
    row = 1
    z = 2
    center = (int(w / 2), int(h / 2), pix[int(w / 2), int(h / 2)])

	# Building elevation profile
    def get_points_on_radius(x, y):
        delta_x = math.fabs(x - center[0])
        delta_y = math.fabs(y - center[1])

        step_x = 1.0 if delta_x >= delta_y else delta_x / delta_y
        step_y = 1.0 if delta_x <= delta_y else delta_y / delta_x

        if delta_x > delta_y:
            if x < center[col]:
                start_x = x
                start_y = y
                finish_x = center[col]
                finish_y = center[row]
            else:
                start_x = center[col]
                start_y = center[row]
                finish_x = x
                finish_y = y
            step_y *= 1 if finish_y > start_y else -1
        else:
            if y < center[row]:
                start_x = x
                start_y = y
                finish_x = center[col]
                finish_y = center[row]
            else:
                start_x = center[col]
                start_y = center[row]
                finish_x = x
                finish_y = y
            step_x *= 1 if finish_x > start_x else -1

        result = []
        tmp_x = start_x
        tmp_y = start_y
        result.append((tmp_x, tmp_y, pix[tmp_x, tmp_y]))
        done = not tmp_x < finish_x if delta_x > delta_y else not tmp_y < finish_y
        while not done:
            tmp_x += step_x
            tmp_y += step_y
            index = (int(tmp_x), int(tmp_y), pix[int(tmp_x), int(tmp_y)])
            result.append(index)
            done = not tmp_x < finish_x if delta_x > delta_y else not tmp_y < finish_y
        return result

	# Get index of max element
    def max_value(indexes):
        if len(indexes) == 0:
            return None
        result = indexes[0]
        for index in indexes:
            if index[z] > result[z]:
                result = index
        return result

	# Get distance from center
    def get_length(x, y):
        return math.sqrt(math.pow(x - center[col], 2) + math.pow(y - center[row], 2))

    image = Image.open(input_filename)
    pix_result = image.load()

	# Iterate over all points
	max_value = 200
    for i in range(w):
        for j in range(h):
            if pix[i, j] == 255:
                pix_result[i, j] = 255
                continue
            pix_result[i, j] = max_value

            points = get_points_on_radius(i, j)
            try:
			    # We don't need center point and target point
                points.pop(len(points) - 1)
                points.pop(0)
            except:
                continue
            hill = max_value(points)
            if not hill:
                pix_result[i, j] = max_value
                continue

            hill_z = hill[z]
            hill_l = get_length(hill[col], hill[row])
            value_l = get_length(i, j)
            value_z = pix[i, j]

            x1 = 0
            y1 = center[2] + elevation
            x0 = hill_l
            y0 = hill_z
            x2 = value_l
            y2 = value_z

            a = y2-y1
            b = x1-x2
            c = x2*y1-x1*y2
            sign = lambda x: (1, -1)[x < 0]
            d = -sign((y2-y1) * (x0 - x1) / (x2 - x1) + y1 - y0) * math.fabs(a * x0 + b * y0 + c) / math.sqrt(a * a + b * b)

            xd = (b*(b*x0-a*y0)-a*c)/(a*a+b*b)
            yd = (a*(a*y0-b*x0)-b*c)/(a*a+b*b)

            r1 = math.sqrt(math.pow(xd-x1,2)+math.pow(yd-y1,2))
            r2 = math.sqrt(math.pow(x2-xd,2)+math.pow(y2-yd,2))

            fresnel_zone_radius = math.sqrt(wavelength * r1 * r2 / (r1 + r2))
            if fresnel_zone_radius == 0:
                pix_result[i, j] = 255
                continue
            if d < -1 * fresnel_zone_radius:
                d = -1 * fresnel_zone_radius
            if d > fresnel_zone_radius:
                d = fresnel_zone_radius
            s = fresnel_zone_radius - d
            pix_result[i, j] = int(max_value * 0.5 * s / fresnel_zone_radius)
			if max_value / cutoff > pix_result[i, j]:
			    pix_result[i, j] = 255

    image.save(output_filename, tiffinfo=image.tag)
