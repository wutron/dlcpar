# adapted from projects/rasmus/python/compbio/vis/transsvg.py

# python libraries
from collections import defaultdict

# ramus libraries
from rasmus import treelib, svg, util, stats
from compbio import phylo


def draw_stree(canvas, stree, slayout,
               yscale=100,
               stree_width=.8, 
               stree_color=(.4, .4, 1),
               snode_color=(.2, .2, .4),
               slabels=None):

    # set defaults
    if slabels is None:
        slabels = {}
    else:
        import colorsys
        r,g,b = stree_color[:3]
        h,s,v = colorsys.rgb_to_hsv(r,g,b)
        r,g,b = colorsys.hsv_to_rgb(h, .8*s, v)
        label_color = (r, g, b)

    # draw stree branches
    w = yscale * stree_width / 2.0
    for node in stree:
        x, y = slayout[node]
        px, py = slayout[node.parent]
        
        # draw branch
        canvas.polygon([px, py-w,
                        x, y-w,
                        x, y+w,
                        px, py+w], strokeColor=(0, 0, 0, 0),
                       fillColor=stree_color)


    # draw stree nodes
    for node in stree:
        if not node.is_leaf():
            x, y = slayout[node]
            canvas.line(x, y-w, x, y+w, color=snode_color,
                        style="stroke-dasharray: 1, 1")


    # label stree nodes
    for node in stree:
        if node.name in slabels:
            x, y = slayout[node]
            px, py = slayout[node.parent]
            canvas.text(slabels[node.name],
                        (x + px)/2., (y + py)/2.,
                        .7*w,
                        fillColor=label_color,
                        anchor="middle",
                        baseline="central")
    


def draw_tree(tree, stree, extra,
              xscale=100, yscale=100,
              leaf_padding=10, 
              label_size=None,
              label_offset=None,
              font_size=12,
              stree_font_size=20,
              canvas=None, autoclose=True,
              rmargin=10, lmargin=10, tmargin=0, bmargin=0,
              stree_color=(.4, .4, 1),
              snode_color=(.2, .2, .7),
              event_size=10,
              rootlen=None,
              stree_width=.8,
              filename=sys.stdout,
              labels=None,
              slabels=None
              ):

    recon = extra["species_map"]
    loci = extra["locus_map"]
    order = extra["order"]

    # setup color map
    all_loci = sorted(set(loci.values()))
    num_loci = len(all_loci)
    colormap = util.rainbow_color_map(low=0, high=num_loci-1)
    locus_color = {}
    for ndx, locus in enumerate(all_loci):
        locus_color[locus] = colormap.get(ndx)

    # set defaults
    font_ratio = 8. / 11.
    
    if label_size is None:
        label_size = .7 * font_size

    #if label_offset is None:
    #    label_offset = -1

    if sum(x.dist for x in tree.nodes.values()) == 0:
        legend_scale = False
        minlen = xscale

    snames = dict((x, x) for x in stree.leaf_names())

    if labels is None:
        labels = {}
    if slabels is None:
        slabels = {}

    # layout stree
    slayout = treelib.layout_tree(stree, xscale, yscale)

    if rootlen is None:
        rootlen = .1 * max(l[0] for l in slayout.values())

    # setup slayout
    x, y = slayout[stree.root]
    slayout[None] =  (x - rootlen, y)
    for node, (x, y) in slayout.items():
        slayout[node] = (x + rootlen, y  - .5 * yscale)

    # layout tree
    ylists = defaultdict(lambda: [])
    yorders = {}

    # layout speciations and genes (y)
    events = phylo.label_events(tree, recon)
    for node in tree.preorder():
        snode = recon[node]
        event = events[node]
        if event == "spec" or event == "gene":
            yorders[node] = len(ylists[snode])
            ylists[snode].append(node)

    # layout internal nodes (y)
    for node in tree.postorder():
        snode = recon[node]
        event = events[node]
        if event != "spec" and event != "gene":
            v = [yorders[child]
                 for child in node.children]
            yorders[node] = stats.mean(v)

    # layout node (x)

    xorders = {}
    xmax = defaultdict(lambda: 0)
    for node in tree.postorder():
        snode = recon[node]
        event = events[node]
        if event == "spec" or event == "gene":
            xorders[node] = 0
        else:
            v = [xorders[child] for child in node.children]
            xorders[node] = max(v) + 1
        xmax[snode] = max(xmax[snode], xorders[node])

##    # initial order
##    xpreorders = {}
##    for node in tree.postorder():
##        snode = recon[node]
##        event = events[node]
##        if event == "spec" or event == "gene":
##            xpreorders[node] = 0
##        else:
##            v = [xpreorders[child] for child in node.children]
##            xpreorders[node] = max(v) + 1
####        print node.name, xpreorders[node]
##    # hack-ish approach : shift x until order is satisfied
##    def shift(node, x):
##        xpreorders[node] += x
##        for child in node.children:
##            if events[child] != "spec":
##                shift(child, x)
##    satisfied = False
##    while not satisfied:
##        satisfied = True
##        for snode, d in order.iteritems():
##            for plocus, lst in d.iteritems():
##                # test each pair
##                for m, node1 in enumerate(lst):
##                    x1 = xpreorders[node1]
##                    for node2 in lst[m+1:]:
##                        x2 = xpreorders[node2]
####                        print node1, node2, x1, x2
##                        if x2 < x1:
##                            # violation - shift all descendants in the sbranch
##                            satisfied = False
####                            print 'violation', node1, node2, x1, x2, x1-x2+1
##                            shift(node2, x1-x2+1)
##                            break
##    # finally, "normalize" xorders
##    xorders = {}
##    xmax = defaultdict(lambda: 0)
##    for node in tree.postorder():
##        snode = recon[node]
##        xorders[node] = xpreorders[node]
##        xmax[snode] = max(xmax[snode], xorders[node])
####        print node.name, xpreorders[node]

    # setup layout
    layout = {None: slayout[None]}
    for node in tree:
        snode = recon[node]
        nx, ny = slayout[snode]
        px, py = slayout[snode.parent]

        # calc x
        frac = (xorders[node]) / float(xmax[snode] + 1)
        deltax = nx - px
        x = nx - frac * deltax

        # calc y
        deltay = ny - py
        slope = deltay / float(deltax)
        deltax2 = x - px
        deltay2 = slope * deltax2
        offset = py + deltay2
        
        frac = (yorders[node] + 1) / float(max(len(ylists[snode]), 1) + 1)
        y = offset + (frac - .5) * stree_width * yscale

        
        layout[node] = (x, y)

##        if y > max(l[1] for l in slayout.values()) + 50:
##            print nx, ny
##            print px, py
##            print offset, frac
##            print ylists[snode], yorders[node]
##            print node, snode, layout[node]


    # layout label sizes
    max_label_size = max(len(x.name)
        for x in tree.leaves()) * font_ratio * font_size
    max_slabel_size = max(len(x.name)
        for x in stree.leaves()) * font_ratio * stree_font_size
    
 
    xcoords, ycoords = zip(* slayout.values())
    maxwidth = max(xcoords) + max_label_size + max_slabel_size
    maxheight = max(ycoords) + .5 * yscale
    
    
    # initialize canvas
    if canvas is None:
        canvas = svg.Svg(util.open_stream(filename, "w"))
        width = int(rmargin + maxwidth + lmargin)
        height = int(tmargin + maxheight + bmargin)
        
        canvas.beginSvg(width, height)
        canvas.beginStyle("font-family: \"Sans\";")
        
        if autoclose == None:
            autoclose = True
    else:
        if autoclose == None:
            autoclose = False

    canvas.beginTransform(("translate", lmargin, tmargin))
    
    draw_stree(canvas, stree, slayout,
               yscale=yscale,
               stree_width=stree_width, 
               stree_color=stree_color,
               snode_color=snode_color,
               slabels=slabels)

    # draw stree leaves
    for node in stree:
        x, y = slayout[node]
        if node.is_leaf():
            canvas.text(snames[node.name], 
                        x + leaf_padding + max_label_size,
                        y+stree_font_size/2., stree_font_size,
                        fillColor=snode_color)


    # draw tree
    for node in tree:
        x, y = layout[node]
        px, py = layout[node.parent]

        if node.parent:
            color = locus_color[loci[node.parent]]
        else:
            color = locus_color[loci[tree.root]]

        canvas.line(x, y, px, py, color=color)

    # draw tree names
    for node in tree:
        x, y = layout[node]
        px, py = layout[node.parent]        

        if node.is_leaf():
            canvas.text(node.name, 
                        x + leaf_padding, y+font_size/2., font_size,
                        fillColor=(0, 0, 0))
            
        if node.name in labels:
            canvas.text(labels[node.name], x, y, label_size, fillColor=(0,0,0))



    # draw events
    for node in tree:
        if node.parent:
            locus = loci[node]
            plocus = loci[node.parent]
            
            if locus != plocus:
                color = locus_color[locus]
                x, y = layout[node]
                o = event_size / 2.0

                canvas.rect(x - o, y - o, event_size, event_size,
                            fillColor=color,
                            strokeColor=color)
        
        
    canvas.endTransform()
    
    if autoclose:
        canvas.endStyle()
        canvas.endSvg()
    
    return canvas


