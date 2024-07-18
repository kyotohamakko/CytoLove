CytoLove
The codes for the method called "CytoLove", which was published in "Machine learning-guided reconstruction of cytoskeleton network from Live-cell AFM Images" in "iscience".
The codes is about the method for reconstruction of fiber structure by a connected-particle model from 2D image.

Reconstruction of fiber network of artificial data

Run the script "reconstruct_artifact.m" to reconstruct the graph representation of a hand make crossing fiber.
Image data has been uploaded in the folder "./data". 

Reconstruction of the actin network

Run the script "reconstruct_actin.m" to reconstruct the graph representation of actin.
By setting actin_type = "cortex" of actin_type = "lamilipodia" to decide which type of actin to reconstruct.
Image data has been uploaded in the folder "./data". 

Angle distribution of actin

After executing "reconstruct_actin.m", result will be saved in the folder "./save". Run the script "angle_distribution.m" to read the reconstruction result and calculate the angle distribution.

For more image data, please refer to cortex1.zip, cortex2.zip and lamellipodia.zip.
