import numpy as np  #NumPy is a general-purpose array-processing package
from odbAccess import * # module to have access to odb output files of Abaqus simulations
import os
import sys
#This module provides a number of functions and variables that can
#be used to manipulate different parts of the Python runtime environment

# print(sys.argv)
name        = sys.argv[-2]        # sys.argv, contains the command-line arguments passed to the script
sim_step    = sys.argv[-1] 
StepName    = 'Step-1'            # Name of the step of interest
Variable    = 'U'                 # Name of the Variable 1
PartName    = 'PART-1-1'          # Name of the Part
nodeSetName = 'ALL_CELL_NODES'    # If the node set is 'all nodes', we do not need to call the nodes 

odb         = openOdb(name)
Steps       = odb.steps[StepName]
stepFrame   = Steps.frames[-1]
subset      = odb.rootAssembly.instances[PartName].nodeSets[nodeSetName].nodes
Field       = stepFrame.fieldOutputs[Variable].getSubset(position = NODAL)
fieldValues = Field.values

# NameOfFile  = 'displacements.txt'
# FileResults = open(NameOfFile,'w')
# for i in range (0,len(subset)):
	# data = fieldValues[subset[i].label-1].data
	# FileResults.write('%1.0f\t %10.7E\t %10.7E\t %10.7E\n' %(fieldValues[subset[i].label-1].nodeLabel , data[0], data[1], data[2]))
# FileResults.close()


# Open a file to write the results
file_name   = "displacements_{}.dat".format(str(sim_step).zfill(5))    # fill with zeros to reach 5 positions
FileResults = open(file_name, 'w')
# Write the header of the file
# FileResults.write('Node Label\t U1\t U2\n')

# Write the values of the field
for i in range(len(subset)):
    # print("Iteration:", i)
    # print("Subset:", subset[i])
    # print("Node Label:", fieldValues[subset[i].label-1].nodeLabel)
    data = fieldValues[i].data
    #magn = fieldValues[i].magnitude
    # print("Data:", data)
    FileResults.write('%1.0f\t %10.7E\t %10.7E\n' %(fieldValues[subset[i].label-1].nodeLabel , data[0], data[1]))

# Close the file
FileResults.close()




Variable1    = 'COORD'                               # Integration point coordinates, requested on element output
#Field1       = stepFrame.fieldOutputs[Variable1].getSubset(position=INTEGRATION_POINT)
Field1       = stepFrame.fieldOutputs[Variable1].getSubset(position=CENTROID)
fieldValues1 = Field1.values

Variable2    = 'S'                               # Stress [Sx Sy Sz Sxy]
#Field2       = stepFrame.fieldOutputs[Variable2].getSubset(position=INTEGRATION_POINT)
Field2       = stepFrame.fieldOutputs[Variable2].getSubset(position=CENTROID)
fieldValues2 = Field2.values

file_name   = "stress_{}.dat".format(str(sim_step).zfill(5))    # fill with zeros to reach 5 positions
FileResults = open(file_name, 'w')
# Write the header of the file
# FileResults.write('int_point_x\t int_point_y\t Sx\t Sy\t Sxy\t Smaxprin\t Sminprin\n')

# contador = 0;

for i in range(len(fieldValues1)):
    # contador = contador+1
    data_1 = fieldValues1[i].data   # vector [x y]
    data_2 = fieldValues2[i].data   # vector [Sx Sy Sz Sxy]
    data_3 = fieldValues2[i].maxPrincipal   # vector [Smax princ]
    data_4 = fieldValues2[i].minPrincipal   # vector [Smin princ]
    
    # print(contador)    
    # print(data_3)
    
    FileResults.write('%10.7E\t %10.7E\t %10.7E\t %10.7E\t %10.7E\t %10.7E\t %10.7E\n' %(data_1[0], data_1[1], data_2[0], data_2[1], data_2[3], data_3, data_4))

    # print('Node label:', nodeLabel, 'Stress:', stress)
    #print('coordx:', data_1[0],'coordy:', data_1[1])
    
FileResults.close()

