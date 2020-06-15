"""
colors.py
Color classes and functions
"""

#=============================================================================
# Color maps


# common colors
RED    = ( 1,  0,  0,  1)
ORANGE = ( 1, .5,  0,  1)
YELLOW = ( 1,  1,  0,  1)
GREEN  = ( 0,  1,  0,  1)
BLUE   = ( 0,  0,  1,  1)
PURPLE = ( 1,  0,  1,  1)
BLACK  = ( 0,  0,  0,  1)
GREY   = (.5, .5, .5,  1)
WHITE  = ( 1,  1,  1,  1)


class ColorMap(object):
    """ColorMap maps values on the real line to colors"""

    def __init__(self, table=None):
        """
        'table' should be the following format:

        [
          [val1, color1],
          [val2, color2],
          [val3, color3],
          ...etc..
        ]

        Values bewteen val1 and val2 will be assigned a blend of color1 and color2.
        value-color pairs can be specified in any order within table.
        """
        if table is None:
            table = []
        self.table = table
        self.table.sort(key=lambda x: x[0])


    def get(self, value):
        """Returns values in [0,1]"""

        assert self.table

        # determine where color falls in table
        for i in xrange(len(self.table)):
            if value <= self.table[i][0]:
                break
        if i > 0:
            i -= 1

        if value <= self.table[i][0]:
            # return lower bound color
            return self.table[i][1]
        if value >= self.table[i+1][0]:
            # return upper bound color
            return self.table[i+1][1]

        # blend two nearest colors
        part = value - self.table[i][0]
        tot = float(self.table[i+1][0] - self.table[i][0])
        weight1 = (tot - part) / tot
        weight2 = part/tot

        newcolor = []
        color1 = self.table[i][1]
        color2 = self.table[i+1][1]
        for j in xrange(len(color1)):
            newcolor.append(weight1 * color1[j]
                            + weight2 * color2[j])
        return newcolor


    def get_int(self, value):
        """Get integer value in [0,255]"""
        return [int(x*255) for x in self.get(value)]


def get_webcolor(color, maxval=1):
    """Get Web RGB hex"""
    colstr = "#"
    for i in color:
        h = hex(int(i * 255.0 / maxval))[2:]
        if len(h) == 1:
            h = "0" + h
        colstr += h
    return colstr


def rainbow_color_map(data=None, low=None, high=None, alpha=None):
    """Return rainbow color map"""
    if data is not None:
        low = min(data)
        high = max(data)
    assert low is not None and high is not None

    if not alpha:
        alpha = 1

    def get_color(color, alpha):
        return color[:3] + (alpha,)

    return ColorMap([[low,              get_color(BLUE, alpha)],
                     [.5*low+.5*high,   get_color(GREEN, alpha)],
                     [.25*low+.75*high, get_color(YELLOW, alpha)],
                     [high,             get_color(RED, alpha)]])
