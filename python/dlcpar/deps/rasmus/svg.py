"""
svg.py
SVG drawing
"""

# python libraries
import os


#=============================================================================
# colors

def _color2string(color):
    """Convert color to string"""
    return "rgb(%d,%d,%d)" % (int(255 * color[0]),
                              int(255 * color[1]),
                              int(255 * color[2]))


def _color_fields(stroke_color, fill_color):
    """Return string to set color fields"""

    txt = ""

    if stroke_color:
        if isinstance(stroke_color, str):
            stroke_str = stroke_color
            stroke_op = None
        else:
            stroke_str = _color2string(stroke_color)
            if len(stroke_color) > 3:
                stroke_op = stroke_color[3]
            else:
                stroke_op = None

        txt += "stroke='%s' " % stroke_str
        if stroke_op is not None:
            txt += "stroke-opacity='%f' " % stroke_op

    if fill_color:
        if isinstance(fill_color, str):
            fill_str = fill_color
            fill_op = None
        else:
            fill_str = _color2string(fill_color)
            if len(fill_color) > 3:
                fill_op = fill_color[3]
            else:
                fill_op = None

        txt += "fill='%s' " % fill_str
        if fill_op is not None:
            txt += "fill-opacity='%f' " % fill_op

    return txt

# common colors
#         r   g   b   a
RED    = (1., 0., 0., 1)
ORANGE = (1., .5, 0., 1)
YELLOW = (1., 1., 0., 1)
GREEN  = (0., 1., 0., 1)
BLUE   = (0., 0., 1., 1)
PURPLE = (1., 0., 1., 1)
BLACK  = (0., 0., 0., 1)
GREY   = (.5, .5, .5, 1)
WHITE  = (1., 1., 1., 1)
NULL   = (0., 0., 0., 0)


#=============================================================================

class Svg(object):
    """Class for writing SVG files"""

    def __init__(self, stream):
        self.out = stream

    def close(self):
        """Close stream"""
        self.out.close()


    def begin_svg(self, width, height):
        """Begin svg"""
        self.out.write(
            """<?xml version='1.0' encoding='UTF-8'?>
            <!DOCTYPE svg PUBLIC '-//W3C//DTD SVG 1.1//EN'
            'http://www.w3.org/Graphics/SVG/1.1/DTD/svg11.dtd'>\n""")
        self.out.write(
            """<svg width='%d' height='%d'
            xmlns='http://www.w3.org/2000/svg' version='1.1'>\n""" % \
            (width, height))

        # default style
        self.out.write("<g>")


    def end_svg(self, close=True):
        """End svg"""
        self.out.write("</g></svg>")
        if close:
            self.close()


    def linear_gradient(self, name, x1, y1, x2, y2,
                        stops=None, colors=None):
        """Draw linear gradient"""
        if stops is None:
            stops = []
        if colors is None:
            colors = []

        self.out.write('''<defs>
<linearGradient id="%s" x1="%f%%" y1="%d%%" x2="%f%%" y2="%f%%">''' %
                       (name, x1, y1, x2, y2))

        for stop, color in zip(stops, colors):
            if len(color) == 3:
                color = list(color) + [1.0]

            self.out.write('<stop offset="%f%%" style="stop-color:%s; '
                           'stop-opacity:%f"/>' % (stop,
                                                   _color2string(color[:3]),
                                                   color[3]))

        self.out.write('</linearGradient></defs>')


    def radial_gradient(self, name, cx=50, cy=50, r=50, fx=50, fy=50,
                        stops=None, colors=None):
        """Draw radial gradient"""
        if stops is None:
            stops = []
        if colors is None:
            colors = []

        self.out.write('''<defs>
<radialGradient id="%s" cx="%f%%" cy="%d%%" r="%f%%" fx="%f%%" fy="%f%%">''' %
                       (name, cx, cy, r, fx, fy))


        for stop, color in zip(stops, colors):
            if len(color) == 3:
                color = list(color) + [1.0]
            self.out.write('<stop offset="%f%%" style="stop-color:%s; '
                           'stop-opacity:%f"/>' % (stop,
                                                   _color2string(color[:3]),
                                                   color[3]))

        self.out.write('</radialGradient></defs>')


    def write_attr_options(self, color=None, **options):
        """Specify attribute options"""
        if color:
            if len(color) > 3:
                self.out.write("stroke-opacity='%f' stroke='%s' " %
                               (color[3], _color2string(color)))
            else:
                self.out.write("stroke='%s' " % (_color2string(color)))

        for key, val in options.iteritems():
            self.out.write("%s='%s' " % (key, val))


    def line(self, x1, y1, x2, y2, color=None, **options):
        """Draw line"""
        self.out.write(
            """<line x1='%f' y1='%f' x2='%f' y2='%f' """ %
            (x1, y1, x2, y2))

        if color != None:
            options["color"] = color
        self.write_attr_options(**options)

        self.out.write(" />\n")


    def polygon(self, verts, stroke_color=BLACK, fill_color=BLACK, stroke_width=1):
        """Draw polygon"""
        self.out.write(
            "<polygon %s points='" % _color_fields(stroke_color, fill_color))

        for i in xrange(0, len(verts), 2):
            self.out.write("%f,%f " % (verts[i], verts[i+1]))
        self.out.write("' stroke-width='%f'/>\n" % stroke_width)


    def rect(self, x, y, width, height, stroke_color=BLACK, fill_color=BLACK, stroke_width=1):
        """Draw rectangle"""
        self.out.write(
            """<rect x='%f' y='%f' width='%f' height='%f' %s stroke-width='%f'/>\n""" % \
            (x, y, width, height, _color_fields(stroke_color, fill_color), stroke_width))


    def circle(self, x, y, radius, stroke_color=BLACK, fill_color=BLACK, stroke_width=1):
        """Draw circle"""
        self.out.write("<circle cx='%f' cy='%f' r='%f' %s stroke-width='%f'/>\n" % \
            (x, y, radius, _color_fields(stroke_color, fill_color), stroke_width))


    def ellipse(self, x, y, xradius, yradius, stroke_color=BLACK, fill_color=BLACK, stroke_width=1):
        """Draw ellipse"""
        self.out.write("<ellipse  cx='%f' cy='%f' rx='%f' ry='%f' %s stroke-width='%f'/>\n" %\
            (x, y, xradius, yradius, _color_fields(stroke_color, fill_color), stroke_width))


    def text(self, msg, x, y, size, stroke_color=NULL, fill_color=BLACK,
             anchor="start", baseline="auto", angle=0, stroke_width=1):
        """Write text"""

        # http://www.w3.org/TR/SVG/text.html
        # anchor
        #     start | middle | end
        # dominant-baseline
        #     auto | alphabetic | ideographics |
        #     middle | central |
        #     mathematical | hanging |
        #     text-bottom | text-top

        anglestr = "transform='translate(%f,%f) rotate(%f)'" % \
                    (x, y, angle)

        self.out.write(
            "<g %s><text x='0' y='0' font-size='%f' %s stroke-width='%f' text-anchor='%s' style='dominant-baseline: %s'>%s</text></g>\n" % \
            (anglestr, size, _color_fields(stroke_color, fill_color), stroke_width,
             anchor, baseline, msg))


    def text2(self, msg, x, y, size, stroke_color=NULL, fill_color=BLACK,
              anchor="start", baseline="auto", angle=0, stroke_width=1):
        """Write text"""

        if angle != 0:
            anglestr = "" #transform='rotate(%f,0,0)'" % angle
        else:
            anglestr = ""


        self.out.write(
            "<text x='%f' y='%f' font-size='%f' %s stroke-width='%f' text-anchor='%s' style='dominant-baseline: %s' %s>%s</text>\n" % \
            (x, y, size, _color_fields(stroke_color, fill_color), stroke_width,
             anchor, baseline, anglestr, msg))


    def begin_transform(self, *options):
        """Begin transform"""
        self.out.write("<g transform='")

        for option in options:
            key = option[0]
            value = option[1:]

            if key == "scale":
                self.out.write("scale(%f, %f) " % value)

            elif key == "translate":
                self.out.write("translate(%f, %f) " % value)

            elif key == "rotate":
                self.out.write("rotate(%f, %f, %f) " % value)

            else:
                raise Exception("unknown transform option '%s'" % key)

        self.out.write("' >\n")


    def end_transform(self):
        """End transform"""
        self.out.write("</g>\n")


    def begin_style(self, style):
        """Begin style"""
        self.out.write("<g style='%s'>\n" % style)


    def end_style(self):
        """End style"""
        self.out.write("</g>\n")


    def write(self, text):
        """Write text"""
        self.out.write(text)


    def comment(self, msg):
        """Write comment"""
        self.out.write("\n<!-- %s -->\n\n" % msg)


def convert(filename, outfilename=None):
    """Convert from svg to png"""
    if outfilename is None:
        outfilename = filename.replace(".svg", ".png")
    os.system("convert " +filename+ " " +outfilename)
    #os.system("rm " + filename)

    return outfilename


# testing
if __name__ == "__main__":
    svg = Svg(file("out.svg", "w"))

    svg.begin_svg(300, 500)

    svg.comment("MY COMMENT")

    svg.begin_transform(('scale', .5, .5))

    svg.line(0, 0, 100, 100, RED)
    svg.rect(10, 10, 80, 100, BLACK, (0, 1, 1, .5))
    svg.polygon([80, 90,
                 100, 100,
                 60, 100],
                (0, 0, 0, 1),
                (0, 0, 1, .3))

    svg.end_transform()


    svg.begin_style("stroke-width:3")
    svg.begin_transform(('translate', 200, 0))

    svg.line(0, 0, 100, 100, RED)
    svg.rect(10, 10, 80, 100, BLACK, (0, 1, 1, .5))
    svg.polygon([80, 90,
                 100, 100,
                 60, 100],
                (0, 0, 0, 1),
                (0, 0, 1, .3))

    svg.end_transform()
    svg.end_style()

    svg.ellipse(150, 250, 70, 50, BLACK, (.5, .5, .9, 1))
    svg.circle(150, 250, 50, BLACK, RED)
    svg.circle(150, 250, 30, WHITE, BLUE)



    svg.begin_style("font-family: arial")
    svg.begin_transform(('translate', 0, -200))
    svg.begin_transform(('translate', 0, 400))
    svg.text("A", 0, 0, 40, BLUE, (1, 0, 0, .1))
    svg.end_transform()

    svg.begin_transform(('translate', 0, 440))
    svg.text("C", 0, 0, 40, BLUE, (1, 0, 0, .1))
    svg.end_transform()
    svg.end_transform()
    svg.end_style()

    svg.begin_style("font-family: helvetica")
    svg.begin_transform(('translate', 0, -200))
    svg.begin_transform(('translate', 0, 480), ('scale', 1, 1))
    svg.text("T", 3, 0, 40, BLUE, (1, 0, 0, .1))
    svg.end_transform()

    svg.begin_transform(('translate', 0, 520), ('scale', 1, .1))
    svg.text("G", 0, 0, 40, BLUE, (1, 0, 0, .1))
    svg.end_transform()
    svg.end_transform()
    svg.end_style()

    svg.line(35, 200, 35, 400, RED)


    svg.begin_style("font-family: courier")
    svg.text("* FIXED WIDTH    *", 100, 400, 10)
    svg.text("* IS THE DEFAULT *", 100, 410, 10)
    svg.end_style()

    for i in range(0, 300, 10):
        color = (i / 300.0, 0, 1 - i/300.0, 1)
        svg.rect(i, 450, 10, 50, color, color)

    svg.end_svg()

    convert("out.svg")
