import math
import heapq
import numpy
import shapefile

def read_polys():
    sf = shapefile.Reader('land-polygons-complete-4326/land_polygons')
    print('%d records...' % sf.numRecords)
    points = []
    n = 0
    for shape in sf.iterShapes():
        n += 1
        if n % 1000 == 0:
            print('%d / %d...' % (n, sf.numRecords))
        yield shape.points[:-1] # last point is same as first point


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


def it_circular_triplets(points):
    for i in range(len(points)):
        yield points[i], points[(i+1)%len(points)], points[(i+2)%len(points)]


most_twisted = []
for points in read_polys():
    # We need to read each shape twice
    # - first to calculate angular sum and distance
    # - second to find the windingness
    total_outer_angle, total_distance, total_count, sum_outer_angle = 0, 0, 0, 0
    for p, q, r in it_circular_triplets(points):
        a, b, c = [ll_to_3d(lat, lon) for lon, lat in (p, q, r)]
        outer_angle = spherical_angle(a, b, c)
        sum_outer_angle += total_outer_angle * dist(a, b)
        total_outer_angle += outer_angle
        total_distance += dist(a, b)
        total_count += 1
        if dist(a, b) > 1e-2:
            p_lon, p_lat = p
            q_lon, q_lat = q
            print('long distance: %f, %f to %f, %f is %.2fkm' % (p_lat, p_lon, q_lat, q_lon, dist(a, b) / (math.pi/2) * 1e4))
    if total_distance <= 0:
        print('total_distance %f, total_count %d' % (total_distance, total_count))
        continue
    average_outer_angle = sum_outer_angle / total_distance - 0.5 * total_outer_angle
    partial_outer_angle, partial_distance = -average_outer_angle, 0
    max_angle, max_coord = 0, None
    for p, q, r in it_circular_triplets(points):
        a, b, c = [ll_to_3d(lat, lon) for lon, lat in (p, q, r)]
        outer_angle = spherical_angle(a, b, c)
        adjusted_outer_angle = partial_outer_angle - (partial_distance / total_distance) * total_outer_angle
        if abs(adjusted_outer_angle) > abs(max_angle):
            max_angle = adjusted_outer_angle
            max_coord = (numpy.array(p) + numpy.array(q)) / 2
        partial_distance += dist(a, b)
        partial_outer_angle += outer_angle
    heapq.heappush(most_twisted, (abs(max_angle), max_angle, max_coord))
    while len(most_twisted) > 100:
        heapq.heappop(most_twisted)

for _, max_angle, max_coord in sorted(most_twisted, reverse=True):
    lon, lat = max_coord
    print('%5.2f %.4f, %.4f' % (max_angle, lat, lon))

