# CytoLove

![image text](https://github.com/kyotohamakko/CytoLove/tree/main/images/Cyto-Love.png "DBSCAN Performance Comparison")

The codes for the method called **CytoLove**, which was published in our paper [**Machine learning-guided reconstruction of cytoskeleton network from Live-cell AFM Images**](https://www.cell.com/iscience/fulltext/S2589-0042(24)02132-1) accepted by **iscience**. The codes is about the method for reconstruction of fiber structure by a connected-particle model from 2D image.

If you use the code in your work, please cite:

Ju et al., Machine learning-guided reconstruction of cytoskeleton network from live-cell AFM images, iScience (2024), https://doi.org/10.1016/j.isci.2024.110907.

# Reconstruction of fiber network of artificial data

Run the script "reconstruct_artifact.m" to reconstruct the graph representation of a hand make crossing fiber.
Image data has been uploaded in the folder "./data". 

# Reconstruction of the actin network

Run the script "reconstruct_actin.m" to reconstruct the graph representation of actin.
By setting actin_type = "cortex" of actin_type = "lamilipodia" to decide which type of actin to reconstruct.
Image data has been uploaded in the folder "./data". 

# Angle distribution of actin

After executing "reconstruct_actin.m", result will be saved in the folder "./save". Run the script "angle_distribution.m" to read the reconstruction result and calculate the angle distribution.

# For more data of our work

For more image data, please refer to cortex1.zip, cortex2.zip and lamellipodia.zip in "./data".
