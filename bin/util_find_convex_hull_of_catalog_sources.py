#!/usr/bin/env python
# 
import click # https://click.palletsprojects.com/en/8.0.x/quickstart/
import os, sys, re, shutil, json
import numpy as np
import astropy.units as u
from astropy.table import Table
from astropy.wcs import WCS
from astropy.io import fits
from astropy.io import ascii as asciitable
from astropy.coordinates import SkyCoord, ICRS
from astropy.visualization import ImageNormalize, PercentileInterval, AsinhStretch
from shapely.geometry import Polygon
from matplotlib import rcParams
rcParams['font.family'] = 'NGC'
from matplotlib import pyplot as plt
from matplotlib import patheffects as path_effects
from scipy.spatial.qhull import ConvexHull
from scipy.spatial.distance import euclidean
import pylab as p
n = np


def _angle_to_point(point, centre):
    '''calculate angle in 2-D between points and x axis'''
    delta = point - centre
    res = n.arctan(delta[1] / delta[0])
    if delta[0] < 0:
        res += n.pi
    return res


def _draw_triangle(p1, p2, p3, **kwargs):
    tmp = n.vstack((p1,p2,p3))
    x,y = [x[0] for x in zip(tmp.transpose())]
    p.fill(x,y, **kwargs)
    #time.sleep(0.2)


def area_of_triangle(p1, p2, p3):
    '''calculate area of any triangle given co-ordinates of the corners'''
    return n.linalg.norm(n.cross((p2 - p1), (p3 - p1)))/2.


def convex_hull(points, graphic=True, smidgen=0.0075):
    '''Calculate subset of points that make a convex hull around points

    Recursively eliminates points that lie inside two neighbouring points until only convex hull is remaining.

    :Parameters:
        points : ndarray (2 x m)
            array of points for which to find hull
        graphic : bool
            use pylab to show progress?
        smidgen : float
            offset for graphic number labels - useful values depend on your data range

    :Returns:
        hull_points : ndarray (2 x n)
            convex hull surrounding points
    '''
    if graphic:
        p.clf()
        p.plot(points[0], points[1], 'ro')
    n_pts = points.shape[1]
    assert(n_pts > 5)
    centre = points.mean(1)
    if graphic: p.plot((centre[0],),(centre[1],),'bo')
    angles = n.apply_along_axis(_angle_to_point, 0, points, centre)
    pts_ord = points[:,angles.argsort()]
    if graphic:
        for i in np.arange(n_pts):
            p.text(pts_ord[0,i] + smidgen, pts_ord[1,i] + smidgen, \
                   '%d' % i)
    pts = [x[0] for x in zip(pts_ord.transpose())]
    prev_pts = len(pts) + 1
    k = 0
    while prev_pts > n_pts:
        prev_pts = n_pts
        n_pts = len(pts)
        if graphic: p.gca().patches = []
        i = -2
        while i < (n_pts - 2):
            Aij = area_of_triangle(centre, pts[i],     pts[(i + 1) % n_pts])
            Ajk = area_of_triangle(centre, pts[(i + 1) % n_pts], \
                                   pts[(i + 2) % n_pts])
            Aik = area_of_triangle(centre, pts[i],     pts[(i + 2) % n_pts])
            if graphic:
                _draw_triangle(centre, pts[i], pts[(i + 1) % n_pts], \
                               facecolor='blue', alpha = 0.2)
                _draw_triangle(centre, pts[(i + 1) % n_pts], \
                               pts[(i + 2) % n_pts], \
                               facecolor='green', alpha = 0.2)
                _draw_triangle(centre, pts[i], pts[(i + 2) % n_pts], \
                               facecolor='red', alpha = 0.2)
            if Aij + Ajk < Aik:
                if graphic: p.plot((pts[i + 1][0],),(pts[i + 1][1],),'go')
                del pts[i+1]
            i += 1
            n_pts = len(pts)
        k += 1
    return n.asarray(pts) # not working



@click.command()
@click.argument('input_catalog_file', type=click.Path(exists=True))
def main(input_catalog_file):
    click.echo('Opening catalog file %r'%(input_catalog_file))

    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1)
    
    tb = Table.read(input_catalog_file)
    # col_ra = [t for t in tb.colnames if t.upper() in ['RA']]
    # col_dec = [t for t in tb.colnames if t.upper() in ['DEC']]
    col_ra = [t for q in ['RA_DEG', 'RA', 'ALPHA_J2000'] for t in tb.colnames if t.upper() == q.upper()]
    col_dec = [t for q in ['DEC_DEG', 'DEC', 'DELTA_J2000'] for t in tb.colnames if t.upper() == q.upper()]
    assert len(col_ra) > 0
    assert len(col_dec) > 0
    x = tb[col_ra[0]].data
    y = tb[col_dec[0]].data
    if np.ma.isMaskedArray(x) and np.ma.isMaskedArray(y):
        mask = np.logical_and(x.mask == False, y.mask == False)
        x = x[mask].filled()
        y = y[mask].filled()
    points = np.column_stack([x, y])
    print('points.shape', points.shape)
    print('points.mean()', points.mean())
    print('x.max() - x.min()', x.max() - x.min())
    print('y.max() - y.min()', y.max() - y.min())
    
    ax.plot(x, y)
    
    hull = ConvexHull(points)
    vertices = hull.vertices.tolist()
    vertices.append(hull.vertices[0])
    vertices = np.array(vertices).astype(int)
    # perimeter = np.sum([euclidean(t1, t2) for t1, t2 in zip(points[vertices], points[vertices][1:])])
    # print('hull.vertices', hull.vertices)
    # print('hull.simplices', hull.simplices)
    for simplex in hull.simplices:
        ax.plot(points[simplex, 0], points[simplex, 1], 'r-')
    
    polygon_str = ''
    for vertex in vertices:
        if polygon_str != '':
            polygon_str += ' '
        polygon_str += '{} {}'.format(points[vertex][0], points[vertex][1])
    
    print('polygon_str: ')
    print(polygon_str)
    
    # hull_points = convex_hull(points.T, graphic=False)
    # ax.plot(hull_points[0], hull_points[1], 'g-') # not working

    ax.set_xlim(ax.get_xlim()[::-1])
    plt.show(block=True)






# MAIN
if __name__ == '__main__':
    main()


