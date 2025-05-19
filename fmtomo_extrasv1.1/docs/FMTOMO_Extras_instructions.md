# FMTOMO Extras

## INTRODUCTION

This file briefly describes how to install and operate several programs that allow:

- Fully non-linear source relocation together with velocity inversion

- Direct inversion for Vp/Vs ratio

It essentially operates as an add-on package to FMTOMO.

## LIMITATIONS

At this point, the code will only work for local earthquake tomography (i.e. sources lies inside model volume), and with media containing no internal interfaces. If this is too restrictive, then you are advised to use FMTOMO, which does allow for source relocation, but under a locally linear assumption.

## ASSUMPTIONS

It is assumed that you are familiar with the FMTOMO package, which you will need to be familiar with in order to use FMTOMO Extras. It is assumed that FMTOMO is installed and properly linked on your system.
NOTE: most of the scripts are in ksh language. zsh can also be used.

## INSTALLATION

To install, edit line 18 of the <u>*compile_all*</u> script so that it points to your favourite Fortran 90 compiler (it is well tested with gfortran). Then simply execute <u>*compile_all*</u> at the command line. All executables will be placed in the subdirectory bin, *which* can then be placed anywhere and accessed by editing your <u>*.bashrc*</u> file or other startup shell script.
**NOTE**: Although some programs and their actions are similar to that found in the original FMTOMO package (e.g.  <u>*grid3dgi*</u>, <u>*synthdatai*</u>, <u>*gmtslicei*</u>), they have been modified and may have different input requirements. Thus, they have been renamed by placing an "i" at the end. This means that the Extras package can be installed without interfering with the original FMTOMO package, although the Extras package DOES require FMTOMO to be installed in order to work (not vice versa).

## CODE STRUCTURE

Two examples have been provided with the distribution: observed and synthetic; in the former case, we are dealing with an observed P-wave and S-wave dataset, and in the latter case, the corresponding checkerboard resolution test. The directory structure of these examples is as follows:

- In the root directory, there are three subdirectories plus a number of files.

- *invert_p* contains all the FMTOMO files required to do an inversion of the P-wave data for Vp-structure. However, the additional element is that instead of <u>*tomo3d*</u>, a script called <u>*tomo3di*</u> is used, which calls the program <u>*gridloc*</u> to perform a non-linear inversion for source relocation between each inversion step for velocity structure. Thus, executing <u>*tomo3di*</u> will invert all P-arrival times for Vp and the source location.

- *invert_s* is exactly the same as *invert_p* except for S-waves.

- The subdirectory *vpvs* contains all the files necessary to perform the direct inversion of P and S-times for Vp/Vs velocity structure.

## RUNNING THE CODE

Some bash scripts (guided by <u>*run_PS_tomography.sh*</u>) are provided to automatically perform the joint inversion of P, S, and Vp/Vs models. These merge the following described procedure that provides a general outline of how to run the new software, once all input files are in place: 

1) Under *invert_p/mkmodel* and *invert_s/mkmodel*, make sure that <u>*obsdatai*</u> has been run and <u>*otimes.dat*</u> and <u>*src.dat*</u> have been created. Copy these files into *invert_p* and *invert_s* directories, ensuring that <u>*otimes.dat*</u> is renamed <u>*otimesref.dat*</u>. Also, <u>*receivers.in*</u> needs to be copied across too and renamed <u>*receiversref.in*</u> (these are actually the source locations that will be adjusted during the inversion process). As with FMTOMO, <u>*sources.in*</u> needs to be copied across, as well as the velocity/interface grid files and propagation grid information.

2) Once these steps have been done, it should be possible to run <u>*tomo3di*</u> in either *invert_p* or *invert_s* to separately carry out inversion and relocation for the P and S-wave data, respectively.

3) In order to invert for both Vp and Vs and jointly relocate (otherwise, if only step 2 is done, then the P-sources and S-sources will probably be in different locations when if fact that should not be the case if a receiver detected both a P and S-wave from the same source) move to the root directory, and execute <u>*jointimes*</u>; this will create <u>*src.dat*</u> in the root directory, with contains association information for P and S times. Now execute the script <u>*tomo3dj*</u>, which will jointly invert for Vp, Vs, and source locations.

4) The traveltime misfit as a function of iteration will be listed in <u>*residuals.dat*</u> in *invert_p* (for P-times) and *invert_s* (for S-times); this is exactly the same as for the standard FMTOMO package. Likewise, Vp can be visualised using <u>*gmtslicei*</u> in *invert_p/gmtplot*, and Vs can be visualised using <u>*gmtslicei*</u> in *invert_s/gmtplot*.

5) Once <u>*tomo3dj*</u> has completed, you can then invert for Vp/Vs ratio. A linear approach is used based on the output of <u>*tomo3dj*</u>. It therefore runs much more quickly than <u>*tomo3dj*</u>. To do this, enter the directory *vpvs*, run the script <u>*cp_script*</u>, and then <u>*tomo3dps*</u>. The data residual information is contained in <u>*vpvsresids.out*</u>, and the Vp/Vs ratio can be visualised in the standard way in the subdirectory *gmtplot*. Note that the starting Vp/Vs ratio is assumed to be constant throughout the model.
   
   **Note** that steps 2 and 3 are alternative in order to perform independent relocation based on P and S models, respectively, or to consider the same source locations for both.

## NOTES ABOUT VISUALISING THE RESULTS (TRADITIONAL)

FMTOMO uses <u>*gmtslice*</u> to extract vertical and horizontal slices through the model for visualisation by GMT. With the Extras package, this is done in exactly the same way, but an option for plotting the % perturbation now exists. Hence, the code has been renamed <u>*gmtslicei*</u> to avoid confusion. To plot percentage perturbation, edit line 59 of <u>*gmtslicei.in*</u>. Then, rather than using the scripts <u>*plotd*</u>, <u>*plotew*</u> and <u>*plotns*</u>, you can use <u>*plotd*</u>, <u>*plotewp*</u> and <u>*plotnsp*</u>. Note that the code can also be set to produce estimates of location uncertainty. This can be done in <u>*gridlocj.in*</u>, and a file <u>*serror.da*</u>t is produced. Note that the plotting scripts don't currently include an option for plotting uncertainty, but this could be added quite easily. For example, to plot the EW error, you could include in the GMT script something like the following (after copying <u>*serror.dat*</u> to the current directory):

```bash
sort -n -r -k 6,6 serror.dat >! temp.dat
awk '{print $3,$1,$6}' temp.dat >! serrorew.dat
psxy serrorew.dat $bounds1 $proj -Sc0.33 -Ccolor.cpt -W1 -O -K -P >> $psfile
```

For the NS error, the equivalent is

```bash
sort -n -r -k 5,5 serror.dat >! temp.dat
awk '{print $2,$1,$5}' temp.dat >! serrorns.dat
psxy serrorns.dat $bounds1 $proj -Sc0.33 -Ccolor.cpt -W1 -O -K -P >> $psfile
```

For depth error, it would be

```bash
sort -n -r -k 4,4 serror.dat >! temp.dat
awk '{print $2,$1,$4}' temp.dat >! serrord.dat
psxy serrord.dat $bounds1 $proj -Sc0.33 -Ccolor.cpt -W1 -O -K -P >> $psfile
```

## NOTES ABOUT VISUALISING THE RESULTS (All In One)

A folder *visualization* has been added containing: a Bash script, a Python script, and an "empty" *<u>gmtslicei.in</u>* file. The Python script *<u>3D_visualization.py</u>* prepares the input for the Bash script *<u>plot_sections.sh</u>* which, in turn, navigates the P and S folder for producing horizontal and vertical images of the two models. Specifically, it produces SN, EW, and Horizontal images corresponding with the velocity grid nodes (e.g., if the velocity grid is 12 points in latitude, the script will produce 12 EW vertical cross sections). Moreover, it also creates vertical sections in the models following the azimuth and the spacing set in the first lines of the Python code. 

The user can find the results in the *gmtplot* folders in the *invert_p/* and *invert_s/* directories. Both absolute and relative velocity images are produced, so check the presence of *<u>velabsolute.cpt</u>* and *<u>velrelative.cpt</u>*.

The Python script *<u>3D_visualization.py</u>* also produces a folder *output_files* containing the velocity models (P and S), the sources (original and relocated), and the receiver data in *<u>.csv</u>* format.

## SETTING UP A SYNTHETIC TEST

Some bash scripts (guided by <u><em>run_PS_synthetic.sh</em></u>) are provided to automatically perform the joint inversion of P, S, and Vp/Vs synthetic models. These merge the following described procedure (here it is assumed that you have a real data example, and now want to run a resolution test):

1. In *invert_p/mkmodel*, edit <u>*grid3dgi.in*</u> so that the synthetic structure is created (be it checkerboard, spikes, or random structure).

2. Copy the grid files to *invert_p* and run <u>*fm3d*</u>, which creates <u>*arrivals.dat*</u>

3. In *invert_p/mkmodel*, copy the newly created <u>*arrivals.dat*</u> to the current directory and execute <u>*synthdatai*</u>. Then copy <u>*src.dat*</u> and <u>*otimes.dat*</u> to *invert_p*, ensuring that in the latter case the file is renamed to <u>*otimesref.dat*</u>.

4. Repeat steps 1-3 for the S-wave data, i.e., under *invert_s*

5. Move to the root location of the synthetic test, execute <u>*jointimes*</u>, and then <u>*tomo3dj*</u>.

6. For Vp/Vs ratio, enter *vpvs* and execute <u>*cp_script*</u> followed by <u>*tomo3dps*</u>. Note that the file <u>*vgridvpvs.in*</u>, which is also created (in addition to the solution that comes from the direct inversion for Vp/Vs, which is located in <u>*vgrids.in*</u>), contains the result of simply dividing the Vp solution (from *invert_p*) with the Vs solution (from *invert_s*). This is useful for the purposes of comparison.

## SELECTION OF PAPERS THAT HAVE USED FMTOMO EXTRAS

- Pilia S., N. Rawlinson, N. G. Direen, P. R. Cummins, and N. Balfour (2013), Structural controls on localized intraplate deformation and seismicity in southern Australia: insights from local earthquake tomography of the Flinders Ranges, Journal of Geophysical Research: Solid Earth, 118, 2176–2190, https://doi.org/10.1002/jgrb.50168.

- Wilks M., Rawlinson N., Kendall J.-M., Nowacki A., Biggs J., Ayele A. and Wookey J. (2020). The Coupled Magmatic and Hydrothermal Systems of the Restless Aluto Caldera, Ethiopia. Frontiers in Earth Sciences, 8, https://doi.org/10.3389/feart.2020.579699

- Zenonos, A., DeSiena, L., Widiyantoro, S., Rawlinson, N. (2020). Direct inversion of S-P differential arrival times for  Zenonos, A., DeSiena, L., Widiyantoro, S., Rawlinson, N. (2020). Direct inversion of S-P differential arrival times for Vp/Vs ratio in SE Asia. Journal of Geophysical Research: Solid Earth, 125, e2019JB019152. https://doi.org/10.1029/2019JB019152 

- Gauntlett, M., Hudson, T., Kendall, J.-M., Rawlinson, N., Blundy, J., Lapins, S., et al. (2023). Seismic tomography of Nabro caldera, Eritrea: Insights into the magmatic and hydrothermal systems of a recently erupted volcano. Journal of Geophysical Research: Solid Earth, 128, e2022JB025742. https://doi.org/10.1029/2022JB025742 

- Han, J., Rawlinson, N., Greenfield, T., White, R. S., Brandsdóttir, B., Winder, T., & Drouin, V. (2024). Evidence of a shallow magma reservoir beneath Askja caldera, Iceland, from body wave tomography. Geophysical Research Letters, 51, e2023GL107851. https://doi.org/10.1029/2023GL107851 
