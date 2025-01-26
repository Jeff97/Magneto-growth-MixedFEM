from abaqus import *
from abaqusConstants import *
from odbAccess import *
import os

# Path to your ODB file
odb_path = 'T:/Abaqus-Temp/20240506MagnetoGrowth/SurfacePattern/3D-WithoutAir/SubFilm-w60-h345-f3-g11.odb'

# Open the ODB file
o1 = session.openOdb(name=odb_path)
session.viewports['Viewport: 1'].setValues(displayedObject=o1)
session.viewports['Viewport: 1'].makeCurrent()
odb = session.odbs[odb_path]

# Define the file path for the output
folder_name = 'OutputEnergy-g11-Film/'
if not os.path.exists(folder_name):
    os.makedirs(folder_name)

######################################################################

# Extract XY data
xy_data_list = session.xyDataListFromField(
    odb=odb,
    outputPosition=ELEMENT_CENTROID,
    variable=(('UVARM17', INTEGRATION_POINT), ('UVARM18', INTEGRATION_POINT)),
    elementSets=("PART-1-1.FILM-DUMMY",)
)

# Ensure the output is a proper XYData format
if xy_data_list:
    # Iterate through the XY data list
    for xy_data in xy_data_list:
        # Extract element label (assuming it is stored in the name or metadata)
        element_label = xy_data.description.split("Element ")[-1].split(",")[0]

        # Sanitize the element label to remove illegal characters
        sanitized_label = ''.join(c for c in element_label if c.isalnum() or c in ('-', '_'))

        # Create XYData object directly
        xy_data_obj = session.XYData(
            name=sanitized_label,
            data=xy_data.data,
            sourceDescription='Extracted from ODB'
        )

        # Define output file name based on element label
        output_file_individual = os.path.join(folder_name,  sanitized_label + ".csv")

        # Write to file
        session.writeXYReport(fileName=output_file_individual, xyData=(xy_data_obj,))
else:
    print("No XY data extracted.")

######################################################################

# Define the file path for the output
folder_name = 'OutputEnergy-g11-Substrate/'
if not os.path.exists(folder_name):
    os.makedirs(folder_name)

######################################################################

# Extract XY data
xy_data_list = session.xyDataListFromField(
    odb=odb,
    outputPosition=ELEMENT_CENTROID,
    variable=(('UVARM17', INTEGRATION_POINT), ('UVARM18', INTEGRATION_POINT)),
    elementSets=("PART-1-1.SUBSTRATE-DUMMY",)
)

# Ensure the output is a proper XYData format
if xy_data_list:
    # Iterate through the XY data list
    for xy_data in xy_data_list:
        # Extract element label (assuming it is stored in the name or metadata)
        element_label = xy_data.description.split("Element ")[-1].split(",")[0]

        # Sanitize the element label to remove illegal characters
        sanitized_label = ''.join(c for c in element_label if c.isalnum() or c in ('-', '_'))

        # Create XYData object directly
        xy_data_obj = session.XYData(
            name=sanitized_label,
            data=xy_data.data,
            sourceDescription='Extracted from ODB'
        )

        # Define output file name based on element label
        output_file_individual = os.path.join(folder_name,  sanitized_label + ".csv")

        # Write to file
        session.writeXYReport(fileName=output_file_individual, xyData=(xy_data_obj,))
else:
    print("No XY data extracted.")


# Close the ODB file to release resources
odb.close()
