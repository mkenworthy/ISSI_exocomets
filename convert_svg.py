"""
--------------------------------------------------
Script to generate custom marker from an SVG image
--------------------------------------------------

Heavily referenced from this tutorial: 
https://petercbsmith.github.io/marker-tutorial.html
"""

import matplotlib.pyplot as plt
import numpy as np

import matplotlib as mpl
from svgpathtools import svg2paths
from svgpath2mpl import parse_path

def generate_marker_from_svg(svg_path):
    image_path, attributes = svg2paths(svg_path)

    image_marker = parse_path(attributes[0]['d'])

    image_marker.vertices -= image_marker.vertices.mean(axis=0)

    image_marker = image_marker.transformed(mpl.transforms.Affine2D().rotate_deg(180))
    image_marker = image_marker.transformed(mpl.transforms.Affine2D().scale(-1,1))

    return image_marker

def demo_generation(svg_path):
    marker1 = generate_marker_from_svg(svg_path)

    x = np.linspace(0,2*np.pi,10)
    y = np.sin(x)

    plt.plot(x,y,'o',marker=marker1,markersize=30)

    plt.show()

if __name__ == "__main__":
    demo_generation("comet2.svg")
