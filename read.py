import math
import numpy
import shapefile

def read_poly(shapes):
    p = None
    while True:
        shape = next(shapes)
        if p is None:
            p = shape.points[0]
        last_q = p
        for q in shape.points[:-1]:
            if q != last_q:
                yield q
            last_q = q
        if shape.points[-1] == p:
            # finally closed
            break


def read_polys():
    sf = shapefile.Reader('coastlines-split-4326/lines')
    print('%d records...' % sf.numRecords)
    shapes = sf.iterShapes()
    points = []
    while True:
        yield read_poly(shapes)


def ll_to_3d(lat, lon):
    lat *= math.pi / 180
    lon *= math.pi / 180
    x = math.cos(lat) * math.cos(lon)
    z = math.cos(lat) * math.sin(lon)
    y = math.sin(lat)
    return numpy.array([x, y, z])


def mag(v):
    return numpy.dot(v, v)**0.5


def dist(u, v):
    return mag(u-v)


def spherical_angle(a, b, c):
    n = numpy.cross(b-a, c-b) / (mag(b-a) * mag(c-b))
    alpha = math.asin(numpy.dot(n, b))
    if numpy.dot(b-a, c-b) >= 0:
        return alpha
    else:
        return numpy.fmod(2*math.pi - alpha, 2*math.pi) - math.pi


def it_circular_triplets(it):
    first_a = a = next(it)
    first_b = b = next(it)
    for c in it:
        yield (a, b, c)
        a, b = b, c
    yield (a, b, first_a)
    yield (b, first_a, first_b)


all_time_best_delta = 0.0

n = 0
for points1, points2 in zip(read_polys(), read_polys()):
    # We need to read each shape twice
    # - first to calculate angular sum and distance
    # - second to find the windingness
    n += 1
    if n % 1000 == 0:
        print('%d...' % n)
    total_outer_angle, total_distance, total_count = 0, 0, 0
    for p, q, r in it_circular_triplets(points1):
        a, b, c = [ll_to_3d(lat, lon) for lon, lat in (p, q, r)]
        outer_angle = spherical_angle(a, b, c)
        total_outer_angle += outer_angle
        total_distance += dist(a, b)
        total_count += 1
    partial_outer_angle, partial_distance = 0, 0
    min_angle, max_angle = float('inf'), float('-inf')
    min_coord, max_coord = None, None
    for p, q, r in it_circular_triplets(points2):
        a, b, c = [ll_to_3d(lat, lon) for lon, lat in (p, q, r)]
        outer_angle = spherical_angle(a, b, c)
        partial_outer_angle += outer_angle
        partial_distance += dist(a, b)
        adjusted_outer_angle = partial_outer_angle - (partial_distance / total_distance) * total_outer_angle
        if adjusted_outer_angle > max_angle:
            max_angle, max_coord = adjusted_outer_angle, r
        if adjusted_outer_angle < min_angle:
            min_angle, min_coord = adjusted_outer_angle, r
    if max_angle - min_angle > all_time_best_delta:
        all_time_best_delta = max_angle - min_angle
        max_lon, max_lat = max_coord
        min_lon, min_lat = min_coord
        print('%f (%d):      %f, %f      %f, %f' % (all_time_best_delta, total_count, min_lat, min_lon, max_lat, max_lon))
