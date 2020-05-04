# This script is to be used in ArcGIS for plotting targets within the analyzed grid zones
# depending on M3. The process (1) uses a sliding window (kernel) performing
# analyses on groups of cells for a considered grid zone, then (2) plot targets given the
# sum of fav. values calculated within the window and the minimum fav. value for each cell
# considered in this window.

import arcpy
from arcpy import env
from arcpy.sa import *
import numpy as np
import pandas as pd

arcpy.CheckOutExtension("spatial")
env.extent = "MAXOF"
env.workspace = "C:/Project/ArcGIS/(grids)"
Coordsystem = "GEOGCS['GCS_ETRS_1989',DATUM['D_ETRS_1989',SPHEROID['GRS_1980',8183294.64,630550.669301,674150.669301,8145094.64]],PRIMEM['Greenwich',0.0],UNIT['Degree',0.0174532925199433]]"
arcpy.env.outCoordinateSystem = Coordsystem
arcpy.env.overwriteOutput = True

# Initialize parameters
input_areas = ["A1","A2","A3","A4","A5","A6","A7","A8","A9","A10"]
cell_size = 100         # 100m on ArcGIS
threshold = 59.53       # sum of the favorability values of model M3
favThresh = 2.358845    # minimum sum of favorability for the model M3

# Plot targets
analyse_grids()


def analyse_grids():
    # Convert threshold value into string name
    if type(favThresh) is float:
        # ArcGIS does not manage . of float when converting
        start = str(favThresh).split(".",1)[0]  # takes name before point
        end = str(favThresh).split(".",1)[1]    # takes name after point
        threshName = start+"_"+end
    elif type(favThresh) is int:
        threshName = str(favThresh)
    else:
        print("Please, enter a non-string threshold value.")

    for gArea in input_areas:
        # Get a grid zone name | identified e.g. by "grid_100_A9" in ArcGIS
        area = str(cell_size) + "_" + gArea
        gridName = r"grid_" + area
        print(gArea + " analyzed...")

        # get favorability and FID values
        grid_cells = arcpy.da.SearchCursor("C:/Project/ArcGIS/(grids)/"+gridName+".dbf",["favv","FID"])
        
        fav_cells,FID_cells = [],[]
        for row in grid_cells:
            fav_cells.append(row[0])
            FID_cells.append(row[1])

        # Grid description
        pDesc = arcpy.Describe(gridName)
        pWidth = pDesc.extent.width         # get witdh
        pHeight = pDesc.extent.height       #get height

        # Add X and Y coordinates in grid table
        fieldnames = [field.name for field in arcpy.ListFields(gridName)]
        if "X_Coord" or "Y_Coord" not in fieldnames:
            print("Calculating coordinates...")
            arcpy.AddField_management(gridName+".dbf", "X_Coord", "DOUBLE",0,field_alias="", field_is_nullable="NON_NULLABLE")
            arcpy.AddField_management(gridName+".dbf", "Y_Coord", "DOUBLE",0,field_alias="", field_is_nullable="NON_NULLABLE")

            arcpy.CalculateField_management(gridName+".dbf","X_Coord","!SHAPE.CENTROID.X!","PYTHON_9.3")
            arcpy.CalculateField_management(gridName+".dbf","Y_Coord","!SHAPE.CENTROID.Y!","PYTHON_9.3")
        else:
            pass

        # Get coordinates
        grid_cells_values = arcpy.da.SearchCursor(gridName,["X_Coord","Y_Coord"])
        grid_cells_XY = []
        for row in grid_cells_values:
            grid_cells_XY.append((row[0],row[1]))    # Coordinates of cells in tuples

        # Redefinition of arrays
        print("Redefining arrays...")
        gridCol = int(pWidth / cell_size)
        gridRow = int(pHeight / cell_size)
        try:
            f_array = np.asarray(fav_cells).reshape(gridRow,gridCol)
            fid_array = np.asarray(FID_cells).reshape(gridRow,gridCol) 
        except:
            try:
                f_array = np.asarray(fav_cells).reshape(gridRow,gridCol+1)
                fid_array = np.asarray(FID_cells).reshape(gridRow,gridCol+1) 
            except:
                try:
                    f_array = np.asarray(fav_cells).reshape(gridRow+1,gridCol)
                    fid_array = np.asarray(FID_cells).reshape(gridRow+1,gridCol) 
                except:
                    f_array = np.asarray(fav_cells).reshape(gridRow+1,gridCol+1)
                    fid_array = np.asarray(FID_cells).reshape(gridRow+1,gridCol+1)

        # Plot the targets
        plot_targets(fav=f_array,xy=grid_cells_XY,fid=fid_array)
        print("Done.\n")


def plot_targets(fav,xy,fid):
    pointCount = 0     # counts number of targets generated

    print("Plotting targets...")
    x_lst_TMP,y_lst_TMP = [],[]
    for x in range(len(fav[:,0])) :
        for y in range(len(fav[0,:])):
            try:
                # sliding window extracting fav. values through the grid
                # window ~model M3 considering e.g. model 6x6 if M3 extent does not exceed 3 cells spatially
                # Notes
                #       the geometry of window is not fundamental in the analysis (1st approximation)
                #       only the number of cells N considered in M3 matters (2nd approximation)
                #       so we take a square and collect best favorability values over N cells
                fav_TMP = fav[x:x+6,y:y+6]
                fid_TMP = fid[x:x+6,y:y+6]

                # Flatten the array
                fav_lst_TMP = fav_TMP.ravel()

                # Sort the flattened array
                fav_lst_sorted = sorted(fav_lst_TMP, reverse=True)

                # Summing ober e.g. 22 first (best) M3 cells
                # if e.g. M3 consist only of 22 cells |
                sumValue = sum(fav_lst_sorted[:22])
                fav_List_TMP = list(fav_lst_sorted[:22])
            except:
                # if we reach the grid border, skip the analysis
                continue

            # Check if (1) the sum of fav. value for N cells considered
            # is above a certain threshold to plot a target at the window location and
            # (2) if a fav_threshold is reach for individual cells
            if (sumValue >= threshold) and all(i >= favThresh for i in fav_List_TMP):
                cell_pos_1 = fid_TMP[0][0]              # cell position where a target point will be created (middle of window)
                cell_pos_2 = fid_TMP[-1][-1]

                # Calculation of coordinates of window center
                XCoord_1 = xy[cell_pos_1-1][0]
                YCoord_1 = xy[cell_pos_1-1][1]
                XCoord_2 = xy[cell_pos_2-1][0]
                YCoord_2 = xy[cell_pos_2-1][1]

                XCoord = (XCoord_1 + XCoord_2) / 2
                YCoord = (YCoord_1 + YCoord_2) / 2
                
                x_lst_TMP.append(XCoord)
                y_lst_TMP.append(YCoord)
                pointCount += 1
            else:
                pass
    print("Number of targets plotted: {0}".format(pointCount))
