from abaqus import *
from abaqusConstants import *
from odbAccess import *
import os

# Define the file path for the output
folder_name = 'Output-WithoutAir-Q1'

# Ensure the folder exists without using exist_ok
if not os.path.exists(folder_name):
    os.makedirs(folder_name)


# Path to your ODB files folder
odb_path = 'T:/Abaqus-Temp/20240506MagnetoGrowth/Cylinder-WithoutAir-Q1/'

# Function to process a single ODB file
def process_odb_file(odb_file_name, output_file_name, Node_select):
    # Open the ODB file
    o1 = session.openOdb(name=odb_path + odb_file_name)
    session.viewports['Viewport: 1'].setValues(displayedObject=o1)
    session.viewports['Viewport: 1'].makeCurrent()
    odb = session.odbs[odb_path + odb_file_name]
    
    # Get the data for the selected node
    TIMEvsCOORD = session.xyDataListFromField(
        odb=odb,
        outputPosition=NODAL,
        variable=(('COORD', NODAL, ((COMPONENT, 'COOR1'), 
                                    (COMPONENT, 'COOR2'), 
                                    (COMPONENT, 'COOR3'))),),
        nodeLabels=(('PART-1-1', (Node_select,)),)
    )
    
    # Write the TIMEvsCOORD to a local file
    output_file_path = os.path.join(folder_name, output_file_name)
    exp_data = TIMEvsCOORD[2]  # Since xyDataListFromField returns a list, we take the 3rd column, which is the COORD3
    
    # Use the session's writeXYReport method to save the data to a CSV file
    session.writeXYReport(fileName=output_file_path, xyData=(exp_data,))
    
    # Close the ODB file to release resources
    odb.close()

# Process each ODB file

# Define the node number you want to track
Node_number = 1236
process_odb_file('Cylinder-g08.odb', 'TimeVSCoord-Cylinder-g08.csv', Node_number)
process_odb_file('Cylinder-g10.odb', 'TimeVSCoord-Cylinder-g10.csv', Node_number)
process_odb_file('Cylinder-g12.odb', 'TimeVSCoord-Cylinder-g12.csv', Node_number)
