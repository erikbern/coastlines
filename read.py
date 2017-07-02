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
    n = 0
    while True:
        n += 1
        if n % 100000 == 0:
            print('%d...' % n)
        if n >= sf.numRecords:
            break # TODO: idk why this is needed
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


for points1, points2 in zip(read_polys(), read_polys()):
    # We need to read each shape twice
    # - first to calculate angular sum and distance
    # - second to find the windingness
    total_outer_angle, total_distance, total_count, sum_outer_angle = 0, 0, 0, 0
    for p, q, r in it_circular_triplets(points1):
        a, b, c = [ll_to_3d(lat, lon) for lon, lat in (p, q, r)]
        outer_angle = spherical_angle(a, b, c)
        sum_outer_angle += total_outer_angle * dist(a, b)
        total_outer_angle += outer_angle
        total_distance += dist(a, b)
        total_count += 1
        if dist(a, b) > 1e-3:
            p_lon, p_lat = p
            q_lon, q_lat = q
            print('long distance: %f, %f to %f, %f is %.2fkm' % (p_lat, p_lon, q_lat, q_lon, dist(a, b) / (math.pi/2) * 1e4))
    average_outer_angle = sum_outer_angle / total_distance - 0.5 * total_outer_angle
    partial_outer_angle, partial_distance = -average_outer_angle, 0
    for p, q, r in it_circular_triplets(points2):
        a, b, c = [ll_to_3d(lat, lon) for lon, lat in (p, q, r)]
        outer_angle = spherical_angle(a, b, c)
        adjusted_outer_angle = partial_outer_angle - (partial_distance / total_distance) * total_outer_angle
        if abs(adjusted_outer_angle) >= math.pi*3:
            lon, lat = (numpy.array(p) + numpy.array(q)) / 2
            print('%f (%d, %f):      %f, %f' % (adjusted_outer_angle, total_count, total_outer_angle, lat, lon))
        partial_distance += dist(a, b)
        partial_outer_angle += outer_angle
