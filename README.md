# CytoLove

![image text](https://github.com/kyotohamakko/CytoLove/blob/main/images/CytoLove.png "Flow of CytoLove")

The codes for the method called **CytoLove**, which was published in our paper [**Machine learning-guided reconstruction of cytoskeleton network from Live-cell AFM Images**](https://www.cell.com/iscience/fulltext/S2589-0042(24)02132-1) accepted by **iscience**. The codes is about the method for recognition and reconstruction of line/fiber structure by a connected-particle model from 2D image. Except for the application for actin filament in biology, this method can be used to track and reconstruct any line/fiber structer in the 2D image.

If you use the code in your work, please cite:

**Ju et al., Machine learning-guided reconstruction of cytoskeleton network from live-cell AFM images, iScience (2024), https://doi.org/10.1016/j.isci.2024.110907.**

# Set the strength of bending energy
Before using this code, please decide the strength of bending energy according to the situation of your line/fiber structure. This method was originally proposed to track and reconstruct actin filament, which has high persistence length. That is, line/fiber structure that almostly straight or with very low curvature. If your line/fiber structure is not straight, please find the **two of** the following lines: 

`E_bend = 100.0*(1.0 - abs(cos(pi*(p.orient-e_new.point.orient)/180)) )^(2*0.5);`

in `reconstruct_artifact.m` or `reconstruct_actin.m`, command it out or decrease the value of the parameter that controls the strength of bending energy (100.0 and 0.5 in the above code) to allow the tracking and reconstruction of high-curvature line/fiber structure. You can see **equation (17)** in the STAR METHODS session in [**our paper**](https://www.cell.com/iscience/fulltext/S2589-0042(24)02132-1) for the mathematical explanation. 


# Reconstruction of fiber network of artificial data

Run the script `reconstruct_artifact.m` to reconstruct the graph representation of a hand make crossing fiber.
Image data has been uploaded in the folder `./data`. 

# Reconstruction of the actin network

Run the script `reconstruct_actin.m` to reconstruct the graph representation of actin.
By setting `actin_type = "cortex"` or `actin_type = "lamilipodia"` to decide which type of actin to reconstruct.
Image data has been uploaded in the folder `./data`. 

# Angle distribution of actin

After executing `reconstruct_actin.m`, result will be saved in the folder `./save`. Run the script `angle_distribution.m` to read the reconstruction result and calculate the angle distribution.

# For more data of our work

For more image data, please refer to `cortex1.zip`, `cortex2.zip` and `lamellipodia.zip` in folder `./data`.
