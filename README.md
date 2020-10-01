# Exact Ligand Cone Angles

## INPUT
**ComplexDataBase1.txt**

Cartesian coordinates for the complex are stored in ComplexDataBase1.txt in the following format. *The apex atom (central atom) must be listed first!*
```
(XN) NAME
1 Z x y z
2 Z x y z
0
```
where N is the index of the complex in the database (to be called upon by ConeAngleDriver.nb), NAME is a unique identifier (for you to identify, it's not used in the program), Z is the atomic symbol (see ConeAnglePackage.nb section for list of supported atoms), and x, y, and z are the Cartesian coordiantes of the corresonding atom. 0 follows each input to separate complexes.


**ComplexDataBase.nb/ComplexDataBase.m**
ComplexDataBase.nb reads in the Cartesian coordinates and formats Apex, Ligands, and XAtoms correctly for ConeAngleDriver.nb. To change the database .txt file where the complexes are stored, edit the *lingandsIn* variable. For example, for a new database NewLigands.txt, the *lingandsIn* variable should be ligandsIn=OpenRead[“NewLigands.txt”]. To solve the cone angles of multiple ligands on a single apex, edit the variable *Ligands*. For example, if atom 1 is the metal center, atoms 2-6 belong to the first ligand, and atoms 7-14 belong to the second ligand, the variable *Ligands* sould be Ligands={Range[2,6],Range[7,14]}. Once edits have been made to ComplexDataBase.nb, the file must be saved as a Mathematica Package with the .m extension to allow ConeAngleDriver.nb to call upon the package.


## RUN FILES
**ConeAngleDriver.nb**
The notebook ConeAngleDriver.nb is the only package the user needs to modify to use the code as is. This interface calls on ComplexDataBase.m to read Cartesian coordinates from ComplexDataBase1.txt and calls ConeAnglePackage.m to calculate the cone angle and visualize the results. The variable *ConeDirectory* should be defined as the directory with the ExactLigandConeAngle files. *ComplexSet* defines the complex(es) of interest. The index should correspond to N in ComplexDataBase1.txt. For example, to find the cone angles of complexes 12-17, ComplexSet = Range[12, 17] would be used. To find the cone angle of only complex 12, the range should be ComplexSet=Range[12,12].

The output can be controlled by the variable *kPrint*:
* *kPrint* = 0, No printing within package
* *kPrint* = 1, Print {ConeAngle, ConeAxis, ConeAtoms} for each ligand
* *kPrint* = 2, For each ligand print both {ConeAngle, ConeAxis, ConeAtoms} and a 3D plot showing the ligand placed inside the cone
* *kPrint* = 3, For each ligand print both {ConeAngle, ConeAxis, ConeAtoms} and a table of van der Waals radii, vertex angles, and Cartesian coordinates
* *kPrint* = 4, For each ligand print both {ConeAngle, ConeAxis, ConeAtoms} and a histogram of candidate three-atom cone angles (°) in the range of (θ<sub>2max</sub>,θ<sub>cm</sub>)

**ConeAnglePackage.nb/ConeAnglePackage.m**
ConeAnglePackage.nb includes all mathematics used to solve for and visualize the cone angle. Atomic symbols H, He, C, N, O, F, Ne, Si, P, S, Cl, Ar, As, Se, Br, Kr, Te, I, Xe, and Fe are supported. To add atoms that are not supported, add the atomic symbol and van der Waals radius in ConeAnglePackage.nb in the *RvdW* variable. *RvdW* is a list of three lists. The first list holds atomic symbols, the second list holds van der Waals radii, and the third list holds the colors the atoms appear in the visualization in *kPrint* = 2. Once edits have been made to ConeAnglePackage.nb, the file must be saved as a Mathematica Package with the .m extension to allow ConeAngleDriver.nb to call upon the package.

## OUTPUT
Program output will appear in ConeAngleDriver.nb along with the total run time. Output level is controlled by kPrint in ConeAngleDriver.nb.

## CITATION
If you find this work useful, please cite our publication: 
* [Jenna A. Bilbrey  Arianna H. Kazez  Jason Locklin  Wesley D. Allen, "Exact ligand cone angles", Journal of Computational Chemistry, 2013, 34 (14), pp. 1189-1197.](https://onlinelibrary.wiley.com/doi/full/10.1002/jcc.23217)



