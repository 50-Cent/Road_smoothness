# Road_smoothness
The road smoothness analysis of a noisy time-series data is based on the smoothness computation on a graph. It is a region-based operation, not a point-based one. A graph is constructed using a couple of points that lie within a certain scale. Therefore, using the change of scale both global and local deformation information can be computed.

The road smoothness analysis takes a feature layer with point feature class as input and gives two separate output layers.

The first output layer exhibits the absolute value smoothness of a road. The smoothness values are coloured encoded in the symbology layer based on warning threshold and critical threshold.

The second output layer shows a set of polygons based on a group of points. This tool automatically save a shapefile for the polygon featureclass in the background. The polygons are colored based on another input symbology layer. As an example, a red polygon indicates that smoothness value of the entire polygon is above a critical threshold. The minimum number of points that are allowed to form a polygon within the scale is set as '3' within the script.
