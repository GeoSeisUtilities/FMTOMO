![FMTOMO.png](/run/media/donato/ARCHIVIO/GeoSisUtilities/FMTOMO_GitHub/FMTOMO/docs/logos/FMTOMO.png)

# Graphical User Interface for performing Local Earthquake Tomography using the FMTOMO code

FMTOMO is a Fortran 90 code for undertaking 3-D traveltime tomography in spherical coordinates. This repository implements a GUI (Graphical User Interface) for easily using the code to perform LET (Local Earthquake Tomography).

## Limitations

The GUI only refers to LET (FMTOMO-extras version), also allowing the fully non-linear relocation of hypocenters and the Vp/Vs inversion. The models consist of one layer within which velocity can vary in 3-D, separated by layer boundaries that have variable geometry. The forward problem is solved using the so-called Fast Marching Method, a grid-based eikonal solver, while the inverse problem is solved using an iterative non-linear approach based on a subspace inversion scheme. Unknowns that can be inverted for include velocity, interface, and hypocenter parameters.

## Documentation

In the *docs* folder, the documentation for this GUI is provided together with the FMTOMO detailed manual, and the FMTOMO-extras relatively brief user guide.

Since the GUI is a collection of shortcuts for the FMTOMO-extras add-on, it is strongly recommended to read the full documentation for understanding the code. 

## Installation

The *<u>compile_all</u>* script provides the full installation of FMTOMO, FMTOMO-extras, and the Python dependencies for the GUI. The instruction manuals contain detailed information about the installation and the dependencies.

## Citation

Since the GUI is a derived product, you are kindly asked to cite the following.

For FMTOMO:

Rawlinson, N., and M. Urvoy (2006), Simultaneous inversion of active and passive source datasets for 3-D seismic structure with application to Tasmania, Geophys. Res. Lett., 33, L24313, doi:10.1029/2006GL028105. 

de Kool, M., Rawlinson, N. and Sambridge, M. (2006), A practical grid-based method for tracking multiple refraction and reflection phases in three-dimensional heterogeneous media, Geophysical Journal International, Volume 167, 253–270, https://doi.org/10.1111/j.1365-246X.2006.03078.x

For FMTOMO-Extras:

Pilia S., N. Rawlinson, N. G. Direen, P. R. Cummins, and N. Balfour (2013), Structural controls on localized intraplate deformation and seismicity in southern Australia: insights from local earthquake tomography of the Flinders Ranges, Journal of Geophysical Research: Solid Earth, 118, 2176–2190, doi:10.1002/jgrb.50168.

For the GUI dependencies:

Beyreuther, M., Barsch, R., Krischer, L., Megies, T., Behr, Y., and Wassermann, J. (May/June 2010), ObsPy: A Python Toolbox for Seismology, Seismological Research Letters, 81 (3), 530-533. [SRL 81:3 Electronic Seismologist](http://www.seismosoc.org/publications/SRL/SRL_81/srl_81-3_es/)

The pandas development team, 2024. pandas-dev/pandas: Pandas. [pandas-dev/pandas: Pandas](https://doi.org/10.5281/ZENODO.3509134)

Harris, C.R., Millman, K.J., Van Der Walt, S.J., Gommers, R., Virtanen, P., Cournapeau, D., Wieser, E., Taylor, J., Berg, S., Smith, N.J., Kern, R., Picus, M., Hoyer, S., Van Kerkwijk, M.H., Brett, M., Haldane, A., Del Río, J.F., Wiebe, M., Peterson, P., Gérard-Marchant, P., Sheppard, K., Reddy, T., Weckesser, W., Abbasi, H., Gohlke, C., Oliphant, T.E., 2020. Array programming with NumPy. Nature 585, 357–362. [Array programming with NumPy | Nature](https://doi.org/10.1038/s41586-020-2649-2)
