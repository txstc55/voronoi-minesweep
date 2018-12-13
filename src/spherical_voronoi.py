from scipy.spatial import SphericalVoronoi
import numpy as np
import random
import math
import sys
from matplotlib import colors
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
import matplotlib.pyplot as plt
from scipy.spatial import SphericalVoronoi
from mpl_toolkits.mplot3d import proj3d


def generate_points(radius):
	x = random.uniform(-radius, radius)
	y = random.uniform(-math.sqrt(radius**2 - x**2), math.sqrt(radius**2 - x**2))
	z = math.sqrt(radius**2 - x**2 - y**2) *random.choice([-1, 1])
	# print(x**2+y**2+z**2)
	return (x, y, z)


n_points = int(sys.argv[1])

centers = set([])


while len(centers)<n_points:
	random_point = generate_points(1)
	centers.add(random_point)

centers = np.array(list(centers))
# print(centers)


points = centers

center = np.array([0, 0, 0])

radius = 1

sv = SphericalVoronoi(points, radius, center)
sv.sort_vertices_of_regions()
def row_to_str(r):
    return " ".join(list(map(lambda x: str(x), list(r))))+"\n"


file = "../src/spherical_voronoi_data"
f = open(file, 'w')

f.writelines(str(len(sv.vertices))+" "+str(len(sv.regions))+"\n")
# print(sv.vertices)


for i in range(len(sv.vertices)):
	f.writelines(row_to_str(sv.vertices[i]))

for i in range(len(sv.regions)):
	f.writelines(row_to_str(sv.regions[i]))
for i in (range(len(points))):
	f.writelines(row_to_str(points[i]))

# # print(len(sv.regions))


# # generate plot
# fig = plt.figure()
# ax = fig.add_subplot(111, projection='3d')
# # plot the unit sphere for reference (optional)
# u = np.linspace(0, 2 * np.pi, 100)
# v = np.linspace(0, np.pi, 100)
# x = np.outer(np.cos(u), np.sin(v))
# y = np.outer(np.sin(u), np.sin(v))
# z = np.outer(np.ones(np.size(u)), np.cos(v))
# ax.plot_surface(x, y, z, color='y', alpha=0.1)
# # plot generator points
# ax.scatter(points[:, 0], points[:, 1], points[:, 2], c='b')
# # plot Voronoi vertices
# ax.scatter(sv.vertices[:, 0], sv.vertices[:, 1], sv.vertices[:, 2],
#                     c='g')
# # indicate Voronoi regions (as Euclidean polygons)
# for region in sv.regions:
#     random_color = colors.rgb2hex(np.random.rand(3))
#     polygon = Poly3DCollection([sv.vertices[region]], alpha=1.0)
#     polygon.set_color(random_color)
#     ax.add_collection3d(polygon)
# plt.show()