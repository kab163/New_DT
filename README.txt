Newest branch of working DT code. This is a combination of some of the summer 2018 code + other add-on code for data partitioning.

HRC_DT.C contains slabbing (vertical columns) of paritioned data for Delaunay. Then it invokes CLEAP to calculate the Delaunay Triangulation.

TODO: 
- remove the boundary triangle of each partition
- piece together the partial solutions to make 1 final solution
- create similar code for slices and cubes
