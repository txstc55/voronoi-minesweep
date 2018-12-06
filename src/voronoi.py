import matplotlib.pyplot as pl
import numpy as np
import scipy as sp
import scipy.spatial
import sys

eps = sys.float_info.epsilon



n_towers = int(sys.argv[1])
towers = np.random.rand(n_towers, 2)
bounding_box = np.array([0., 1., 0., 1.]) # [x_min, x_max, y_min, y_max]

def in_box(towers, bounding_box):
    return np.logical_and(np.logical_and(bounding_box[0] <= towers[:, 0],
                                         towers[:, 0] <= bounding_box[1]),
                          np.logical_and(bounding_box[2] <= towers[:, 1],
                                         towers[:, 1] <= bounding_box[3]))


def row_to_str(r):
    return " ".join(list(map(lambda x: str(x), list(r))))+"\n"


def voronoi(towers, bounding_box):

    # Select towers inside the bounding box
    i = in_box(towers, bounding_box)
    # Mirror points
    points_center = towers[i, :]
    points_left = np.copy(points_center)
    points_left[:, 0] = bounding_box[0] - (points_left[:, 0] - bounding_box[0])
    points_right = np.copy(points_center)
    points_right[:, 0] = bounding_box[1] + (bounding_box[1] - points_right[:, 0])
    points_down = np.copy(points_center)
    points_down[:, 1] = bounding_box[2] - (points_down[:, 1] - bounding_box[2])
    points_up = np.copy(points_center)
    points_up[:, 1] = bounding_box[3] + (bounding_box[3] - points_up[:, 1])
    points = np.append(points_center,
                       np.append(np.append(points_left,
                                           points_right,
                                           axis=0),
                                 np.append(points_down,
                                           points_up,
                                           axis=0),
                                 axis=0),
                       axis=0)
    # Compute Voronoi
    vor = sp.spatial.Voronoi(points)
    # Filter regions
    regions = []
    for region in vor.regions:
        flag = True
        for index in region:
            if index == -1:
                flag = False
                break
            else:
                x = vor.vertices[index, 0]
                y = vor.vertices[index, 1]
                if not(bounding_box[0] - eps <= x and x <= bounding_box[1] + eps and
                       bounding_box[2] - eps <= y and y <= bounding_box[3] + eps):
                    flag = False
                    break
        if region != [] and flag:
            regions.append(region)
    vor.filtered_points = points_center
    vor.filtered_regions = regions


    fig = pl.figure()
    ax = fig.gca()


    # Plot initial points
    ax.plot(vor.filtered_points[:, 0], vor.filtered_points[:, 1], 'b.')



    # Plot ridges points
    for region in vor.filtered_regions:
        vertices = vor.vertices[region, :]
        ax.plot(vertices[:, 0], vertices[:, 1], 'g.')

    # Plot ridges
    for region in vor.filtered_regions:
        vertices = vor.vertices[region + [region[0]], :]
        ax.plot(vertices[:, 0], vertices[:, 1], 'k-')


    # print("Check points")
    # for i in range(len(vor.vertices)):
    #     for j in range(len(vor.vertices)):
    #         if i!=j:
    #             if vor.vertices[i][0]==vor.vertices[j][0] and vor.vertices[i][1]==vor.vertices[j][1]:
    #                 print(i, j)
    # print("Done checking points")

    for i, p in enumerate(vor.vertices):
        ax.annotate(str(i), (p[0], p[1]))
        # print(p)



    return vor

print("Generating voronoi graph")
vor = voronoi(towers, bounding_box)

print("Checking completeness")
while (len(vor.filtered_regions)!=n_towers):
    towers = np.random.rand(n_towers, 2)
    vor = voronoi(towers, bounding_box)



file = "../src/voronoi_data"
f = open(file, 'w')

print("Writing data to "+file)


f.writelines(str(len(vor.vertices)) + " "+ str(len(vor.filtered_regions))+"\n")
for i in range(len(vor.vertices)):
    f.writelines(row_to_str(vor.vertices[i]))

for i in range(len(vor.filtered_regions)):
    f.writelines(row_to_str(vor.filtered_regions[i]))

f.close()



# print(vor.filtered_points)
# print(len(vor.filtered_points))
# print(vor.vertices)
# print(len(vor.vertices))

# print(vor.filtered_regions)
# print(len(vor.filtered_regions))
# print(max(list(map(lambda x: max(x), vor.filtered_regions))))


fig_path = "../src/voronoi_graph.jpg"
pl.savefig(fig_path)
print("Check the graph at "+fig_path)

