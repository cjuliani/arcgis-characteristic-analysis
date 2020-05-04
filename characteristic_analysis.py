# This script must be used within ArcGIS to work on feature datasets and tables.
# Conceptual method: McCammon, R.B., Botbol, J.M., Sinding-Larsen, R., Bowen, R.W., 1983. Characteristics analysis – 1981: final program and a possible discovery. Math. Geol. 15 (1), 59–83.
# Implementation for case study: Juliani, C., Ellefmo, S.L., 2018. Probabilistic estimates of permissive areas for undiscovered
# seafloor massive sulfide deposits on an Arctic Mid-Ocean Ridge. Ore Geol. Rev. 95, 917–930.

from __future__ import division
# In Python 2, 25/100 is zero when performing an integer divison. since the result is less than 1.
# __future__ import division to your script. This will always perform a float division when using the / operator,
# and use // for integer division.
# Another option would be making at least one of the operands a float, e.g. 25.0/100.
# In Python 3, 25/100 is always 0.25.
# This is a problem of integer truncation (i.e., any fractional parts of a number are discarded).

import arcpy.da, arcpy
import numpy as np
# import matplotlib.pyplot as plt
import pandas as pd
from arcpy import env
from scipy.stats import binom
from arcpy.sa import *

# Set numpy arrays to display 3 decimals instead of none after the comma
np.set_printoptions(formatter={'float': lambda x: "{0:0.6f}".format(x)})

input_areas = ["A1",
               "A2",
               "A3",
               "A4",
               "A5",
               "A6",
               "A7",
               "A8",
               "A9",
               "A10"]      # grid zones; i.e. fishnet grids with cell coordinates
# input_areas = ["LC_M3_3"]
mtxModel = "normal"         # normal or M2 | model M2 uses binomial verification (see Appendix in Juliani and Ellefmo, 2018)
tableUpdt = "yes"           # yes or no
cell_size = 100             # cell size from fishnet grids
fldName = "favv"       # field name newly created in gris to save favorability values


run_analysis()


def run_analysis():
    # counter Initialization | to save 1st eigenvector calculated from 1st area, for next areas
    cnt = 0
    for gArea in input_areas:
        # Get a grid zone name | identified e.g. by "grid_100_A9" in ArcGIS
        area = str(cell_size) + "_" + gArea
        gridName = r"grid_" + area
        print(gArea + " analyzed...")

        # References to tables from which values will be extracted
        # Notes:
        #       tables contain spatial data obtained from Focal Statistics given the grids considered.
        #       3 examples of variables given below | statistics considered are in the MEAN form.
        #       FID_: name of field identifying grid cells (ID)
        var_1 = arcpy.da.SearchCursor(gridName + "_var1.dbf", ["FID_", "MEAN"])     # e.g. fault distance calculated in ArcGIS
        var_2 = arcpy.da.SearchCursor(gridName + "_var1.dbf", ["FID_", "MEAN"])
        var_3 = arcpy.da.SearchCursor(gridName + "_var1.dbf", ["FID_", "MEAN"])

        # Lists that will contain column data from dbf files
        var_1_lst, var_2_lst, var_3_lst = [],[],[]

        # Extraction of variables and associated FID
        for row in var_1:
            var_1_lst.append((row[0],row[1]))  # (FID,value)
        for row in var_2:
            var_2_lst.append((row[0],row[1]))
        for row in var_3:
            var_3_lst.append((row[0],row[1]))

        # ----- Define spatial rules
        # Notes:
        #       Transform data in ternary  from given spatial rules determined by the investigator.
        #       For example, if var_1 is fault density, then a value > 20.0 is considered favourable for mineralization
        #           due to permeability increase. Everything > 20.0 is +1.
        #       Value intervals for which we are not sure about the favorability are set to 0.
        #       Unfavourable value are set to -1.

        # Get FID list
        FID_list = []
        for t in range(len(var_1_lst)):
            FID_list.append(var_1_lst[t][0])

        # Rule 1 | e.g. for fault distance
        var_1_TERNARY = []
        for t in range(len(var_1_lst)):
            if var_1_lst[t][1] <= 200.0:
                # favorable distance up to 200m from fault
                var_1_TERNARY.append(1)
            elif 500.0 >= var_1_lst[t][1] > 200.0:
                # uncertain between 200 and 500m
                var_1_TERNARY.append(0)
            else:
                # unfavorable beyond 500m
                var_1_TERNARY.append(-1)

        # Rule 2
        var_2_TERNARY = []
        for t in range(len(var_2_lst)):
            if var_2_lst[t][1] <= 300.0:
                var_2_TERNARY.append(1)
            elif 300.0 < var_2_lst[t][1] <= 800.0:
                var_2_TERNARY.append(0)
            else:
                var_2_TERNARY.append(-1)

        # Rule 3
        var_3_TERNARY = []
        for t in range(len(var_3_lst)):
            if var_3_lst[t][1] <= 500.0:
                var_3_TERNARY.append(1)
            else:
                var_3_TERNARY.append(0)

        # Data management
        main_DataFrame = pd.DataFrame(data={'var1': var_1_TERNARY,
                                            'var2': var_2_TERNARY,
                                            'var3': var_3_TERNARY},
                                      columns=['var1', 'var2', 'var3'])

        # Convertion to matrix
        main_Matrix = pd.DataFrame.as_matrix(main_DataFrame)
        dim = main_Matrix.shape[0]  # rows (observations)
        var = main_Matrix.shape[1]  # columns (variables)

        # ----- Calculate scores
        inputMtx = modMatrix.copy()     # avoid modifying the original matrix

        # Calculates PC1 of 1st grid zone i.e. if counter is 0.
        if cnt == 0:

            # probability matrix for models M1 or M3 considered, do binomial analysis | see Appendix in Juliani and Ellefmo, 2018
            if mtxModel == "M1M3":
                # Calculate p and q, i.e. numbers of +1 and -1 respectively, for each variable
                for i in range(len(inputMtx[0,:])):
                    temp = inputMtx[:, i]
                    temp = np.transpose(temp)
                    x = 'p{}'.format(i + 1)
                    exec("{} = {}".format(x, list(temp).count(1)))      # number of +1 for variable i
                    y = 'q{}'.format(i + 1)
                    exec("{} = {}".format(y, list(temp).count(-1)))     # number of -1 for variable i
                n = len(temp)

                # Calculate the probability of success between each variable
                for i in range(1, len(inputMtx[0, :]) + 1):
                    for j in range(1, len(inputMtx[0, :]) - (i - 1) + 1):
                        j += i - 1
                        a = 'prob{}{}'.format(i, j)
                        b = eval('p{}'.format(i))  # eval() to get the value, not the string name
                        c = eval('p{}'.format(j))
                        d = eval('q{}'.format(i))
                        e = eval('q{}'.format(j))
                        exec("{} = ({} * {} + {} * {})/n**2".format(a, b, c, d,e))
                        #print("prob" + str(i) + str(j))

                # Matrix reconstruction
                probMtx = np.zeros((var, var))  # (3x3)

                # Construction of probability matrix (triangular form)
                for i in range(0, len(inputMtx[0,:])):
                    k = 1 + (i - 1)
                    while k < len(inputMtx[0, :]):
                        probMtx[i, k] = eval('prob{}{}'.format(i, k))
                        k += 1

                # Symmetric shaping of the matrix (by adding missing values)
                for i in range(0, len(inputMtx[0, :]) - 1):
                    k = 1 + (i - 1)
                    while k < len(inputMtx[0, :]):
                        probMtx[k, i] = probMtx[i, k]
                        k += 1

                # Print the prob. matrix to analyze which cells must be kept
                # or removed in ArcGIS for M1 and M3 (see Appendix in Juliani and Ellefmo, 2018)
                print("Probability matrix :", probMtx)

            # probability matrix for model M2
            elif mtxModel == "M2":
                probMtx_pmf = np.zeros((var, var))

                for i in range(1, len(inputMtx[0,:]) + 1):
                    for j in range(1, len(inputMtx[0,:]) - (i - 1) + 1):
                        j += i - 1
                        temp = inputMtx[:, i - 1] + inputMtx[:, j - 1]
                        a = 'r{0}{1}'.format(i, j)
                        counter = list(temp).count(2) + list(temp).count(-2)

                        if counter != 0:
                            exec('{} = {}'.format(a, counter))
                            r = eval('r{0}{1}'.format(i, j)) - 1
                            #
                            end = min(eval('p{}'.format(i)),
                                      eval('p{}'.format(j))) + min(eval('q{}'.format(i)),
                                      eval('q{}'.format(j)))    # take min. because matches cannot be higher than the min. of p?+q?
                            x = range(0,end)                    # possible numbers of matching

                            # probability of x matching for n trials (binomial process)
                            binomVal = sum(binom.pmf(x[:counter], n, eval('prob{0}{1}'.format(i, j))))
                            probMtx_pmf[i - 1, j - 1] = binomVal
                        else:
                            # if counter is 0, it means probability of finding -1 or +1 is null
                            probMtx_pmf[i - 1, j - 1] = 0.0

                # Matrix reconstruction
                for i in range(0, len(inputMtx[0,:]) - 1):
                    for j in range(1, len(inputMtx[0,:]) - (i - 1)):
                        j += i - 1
                        probMtx_pmf[j, i] = probMtx_pmf[i, j]

                # Print the prob. matrix to analyze which cells must be kept
                # or removed in ArcGIS for M2 (see Appendix in Juliani and Ellefmo, 2018)
                print("Probability matrix: ", probMtx_pmf)

            else:
                # if mtxModel does not consider models but regional data, do PCA
                # Notes
                #       the target plotting will be based on the M3 model
                #       and M3 will be analyze on the regional favorability map generated by PCA

                # Method 1 - PCA on covariance matrix
                pca = PCA(n_components=modMatrix.shape[1])
                pca.fit(modMatrix)
                # Position of highest variance:
                ix = list(pca.explained_variance_ratio_).index(max(pca.explained_variance_ratio_))
                pc1_vect = pca.components_[ix]

                # Method 2 - PCA on correlation matrix
                # results = PCAm(modMatrix)
                # results.fracs	# variances
                # results.Wt[0]	# PC1

                print("PC1: " + str(pc1_vect))
                if tableUpdt == "yes":
                    update_grid(gridName,inputMtx,pc1_vect)
        else:
            # if cnt>=1
            pass
            if tableUpdt == "yes":
                update_grid(gridName,inputMtx,pc1_vect)
        cnt = cnt + 1

def update_grid(gridName,mtx,weights):
    # Grid table considered (gArea) to be updated
    tableInput = gridName + ".dbf"

    # Calculate fav. values
    probMtx_M_list = []
    for i in range(len(mtx[:, 0])):
        probMtx_M_list.append(mtx[i,:].dot(weights))

    # Add new field for storing favorability values (favv) if not existing in grid table
    fieldnames = [field.name for field in arcpy.ListFields(gridName)]
    if fldName not in fieldnames:
        arcpy.AddField_management(tableInput, fldName,
                                  "DOUBLE", 0,
                                  field_alias="",
                                  field_is_nullable="NULLABLE")
        print("Field " + str(fldName) + " created.")
    else:
        pass

    # Write favorability values in grid table
    print("Updating...")
    tableUpdate = arcpy.UpdateCursor(gridName)
    i = 0
    for row in tableUpdate:
        row.setValue(fldName, probMtx_M_list[i])
        tableUpdate.updateRow(row)
        i = i + 1
    print("Done.\n")