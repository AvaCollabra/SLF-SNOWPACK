/***********************************************************************************/
/*  Copyright 2009-2010 WSL Institute for Snow and Avalanche Research    SLF-DAVOS */
/***********************************************************************************/
/* This file is part of MeteoIO.
    MeteoIO is free software: you can redistribute it and/or modify
    it under the terms of the GNU Lesser General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    MeteoIO is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU Lesser General Public License for more details.

    You should have received a copy of the GNU Lesser General Public License
    along with MeteoIO.  If not, see <http://www.gnu.org/licenses/>.
*/
#ifndef __MAINPAGE_H__
#define __MAINPAGE_H__

 /**
 * @mainpage Welcome to SNOWPACK
 * @section intro_sec Introduction
 * SNOWPACK is the operational model of the Swiss avalanche warning service and is available as an integrated software package, made of a library that can be reused
 * in third party applications and a standalone model. It simulates the evolution of the snow cover based on meteorological input data, based on the physical modeling
 * of the various processes taking place. International intercomparison studies show that SNOWPACK is successfully applied to alpine, arctic, maritime and continental snow covers.
 *
 * This library is available under LGPL version 3 or above, see <a href="http://www.gnu.org/licenses/lgpl.txt">www.gnu.org</a>. The Visual C++ version uses a BSD-licensed port of getopt for Visual C++, with a \subpage getopt_copyright "BSD copyright notice".
 *
 * @section table_of_content Table of content
 * -# \subpage getting_started "Getting Started"
 * -# End User documentation
 *    -# Model principles
 *        -# \subpage general "General concepts"
 *        -# \subpage references "References"
 *    -# Inputs / Outputs
 *        -# \subpage requirements "Data requirements"
 *        -# \subpage input_formats "Input file formats"
 *        -# \subpage output_formats "Output file formats"
 *    -# Simulation tools
 *        -# \subpage configuration "Configuring a simulation"
 *        -# \subpage sngui_config "Visualization with sngui"
 * -# Programing using %Snowpack
 *        -# \subpage libsnowpack_basics "Programming with libsnowpack"
 * -# Expanding %Snowpack
 *        -# \subpage coding_style "Coding style"
 *        -# \subpage adding_extra_models "Adding extra models"
 */

/**
 * @page getting_started Getting Started
 * After you installed a binary package or compiled and installed Snowpack, you can run your first simulation.
 * Please make sure you properly set the proper environement variables for your operating system:
 *      - on osX: set \em PATH and \em DYLD_FALLBACK_LIBRARY_PATH
 *      - on Linux: set \em PATH and \em LD_LIBRARY_PATH if you install the package to a non-standard location
 *      - on Windows: set \em PATH
 * How to do this (and much more) is explained in the online documentation at https://slfsmm.indefero.net/p/snowpack/page/Getting-started/.
 *
 * Now, let's run an example simulation:
 * -# First, copy the examples as provided with this documentation (the whole \b doc/examples directory) to a directory where you have write access.
 * -# Then edit the configuration file for the chosen simulation (in cfgfiles), for example \b io_res1exp.ini and set the PLUGINPATH key to
 *    contain the path to your input plugins (such as \em Applications/Snowpack-3.0.0/lib/meteoio/plugins on osX, \em /usr/local/lib/meteoio/plugins on Linux
 *    or <i>C:\\Program files\\Snowpack-3.0.0\\lib\\meteoio\\plugins</i> on Windows)
 * -# Select one of the examples and run its start script, for example <b>./run_res1exp.sh</b> (on Windows, you have to open the script file and copy
 *    the last command line it contains into a terminal)
 * -# Once the simulation is finished, the results are available in the \b output directory
 * -# The results can be visualized using the \ref sngui_config "sngui tool" by opening the <b>.pro</b> file that was generated in \b output.
 *
 */

/**
 * @page general General concepts
 * The one-dimensional snow cover model SNOWPACK (Lehning et al., 1999; Bartelt and Lehning, 2002; Lehning et al., 2002a, b),
 * initially programmed in C, was primarily developed for the support of avalanche warning in Switzerland (it has now been ported to C++).
 * However, this physical model is also used for other applications such as permafrost investigations (Lütschg et al., 2003),
 * the assessment of snow – vegetation interactions, climate research (Rasmus and Räisänen, 2003; Bavay et al., 2009), mass- and energy balance
 * calculations for arctic areas (Meirold-Mautner and Lehning, 2003) and calculations of chemical solute transport in snow (Waldner et al., 2003).
 *
 * @section physical_processes Physical processes
 * \image html basics.png "Principal physical processes included in the SNOWPACK model"
 * \image latex basics.eps "Principal physical processes included in the SNOWPACK model" width=0.9\textwidth
 *
 * A graphical review of the physical processes described by the SNOWPACK model is given in the above figure. SNOWPACK is based on a Lagrangian
 * finite element implementation and solves the instationary heat transfer and settlement equations. Phase changes and transport of water vapor and
 * liquid water are included. Special attention is paid to the metamorphism of snow and its connection to mechanical properties such as thermal
 * conductivity and viscosity. At present, SNOWPACK is the most advanced snow cover models worldwide in terms of microstructural detail.
 * Therefore, first attempts are being made to estimate snow cover stability from SNOWPACK simulations (Lehning et al., 2003).
 *
 * @section current_op_usage Operational usage experience
 * SNOWPACK runs operationally
 * on a network of high Alpine automatic weather and snow measurement stations. Presently approximately 100 sites are in operation. These stations
 * measure wind, air temperature, relative humidity, snow depth, surface temperature, ground (soil) temperature, reflected short wave radiation and
 * three temperatures within the snowpack. The measurements are hourly transferred to the SLF and drive the model. SNOWPACK produces supplementary
 * information regarding the state of the snowpack at the sites of the automatic stations. The model is connected to a relational data base which
 * stores the measurements as well as the model results. Validations of SNOWPACK suggest that the calculations are reliable in terms of the mass balance
 * and the energy budget. The implemented snow metamorphism formulations yield reasonable grain types and are able to reproduce important processes
 * such as the formation of depth or surface hoar.
 *
 * @section other_uses Other uses
 * In addition to the stand-alone applications, SNOWPACK is increasingly used in a distributed way (Kuonen et al., 2010).
 * SNOWPACK has been coupled with atmospheric flow and snow drift modules as well as with spatial energy balance models. The coupled models are
 * used to investigate snow deposition and snow cover development in steep terrain (Lehning and others, 2000) and to forecast ski run conditions for racing
 * (however, the current version of SN_GUI does not include the visualization of distributed SNOWPACK calculations).
 *
 * In order to make it easier to integrate %Snowpack in other models, it has been repackaged as a library. You can therefore use the %Snowpack library in another
 * model. More details are given in \subpage libsnowpack_basics "Programming with libsnowpack".
 *
 */

/**
 * @page references References
 * In the following some important papers related to the SNOWPACK model are listed. Additional information can be found on the web:
 * www.slf.ch/lwr/prozessmodelle/aufgaben-en.html.
 * - Bartelt, P.B. and M. Lehning. 2002. A physical SNOWPACK model for Avalanche Warning Services. Part I: Numerical Model. Cold Reg. Sci. Technol, 35(3), 123-145.
 * - Colbeck, S.C., E. Akitaya, R. Armstrong, H. Gubler, J. Lafeuille, K. Lied, D. McClung and E. Morris. 1990. The international classification of seasonal snow on the ground. Wallingford, Oxon, U.K., International Commission on Snow and Ice (ICSI), International Association of Scientific Hydrology. http://www.crrel.usace.army.mil/techpub/CRREL_Reports/reports/Seasonal_Snow.pdf
 * - Lehning, M., P.B. Bartelt, R.L. Brown, C. Fierz and P. Satyawali. 2002a. A physical SNOWPACK model for the Swiss Avalanche Warning Services. Part II: Snow Microstructure. Cold Reg. Sci. Technol, 35(3), 147-167.
 * - Lehning, M., P.B. Bartelt, R.L. Brown, C. Fierz and P. Satyawali. 2002b. A physical SNOWPACK model for the Swiss Avalanche Warning Services. Part III: Meteorological Boundary Conditions, Thin Layer Formulation and Evaluation. Cold Reg. Sci. Technol, 35(3), 169-184.
 * - Lehning, M., C. Fierz and C. Lundy. 2001. An objective snow profile comparison method and its application to SNOWPACK. Cold Reg. Sci. Technol., 33, 253-261.
 * - Lehning, M., J. Doorschot, N. Raderschall and P.B. Bartelt. 2000. Combining snow drift and SNOWPACK models to estimate snow loading in avalanche slopes. In: Snow Engineering, HjorthHansen, Holand, Loset & Norem (eds.), Balkema, 113-122.
 * - Lehning, M., P. Bartelt, R.L. Brown, T. Russi, U. Stöckli and M. Zimmerli. 1999. %Snowpack model calculations for avalanche warning based upon a new network of weather and snow stations. Cold Reg. Sci. Technol., 30(1-3), 145-157.
 * - Lehning, M., C. Fierz, B. Brown, B. Jamieson,. 2003. Modelling Instability for the Snow Cover Model SNOWPACK. Ann. of Glaciol., in press.
 * - Lütschg, M., P.B. Bartelt, M. Lehning, V. Stoeckli and W. Haeberli. 2003. Numerical simulation of the interaction processes between snow cover and alpine permafrost. Proceedings of the 8th International Conference on Permafrost, 21-25 July, 2003, Zurich, Switzerland, in press.
 * - Meirold-Mautner, I. and M. Lehning. 2003. Measurement and modeling of the solar shortwave fluxes in snow on Summit/Greenland. J. of Glaciol., in press.
 * - Rasmus, S., J. Räisänen and M. Lehning. 2003. Estimating snow conditions in Finland in the late 21st century using the SNOWPACK–model with regional climate scenario data as input. Ann. of Glaciol., in press.
 * - Spreitzhofer, G, Lehning, M. and C. Fierz. 2004. SN_GUI: A graphical user interface for snowpack modelling. Environmental Modelling & Software, submitted.
 *
 */

/**
 * @page getopt_copyright BSD copyright notice
 * This copyright notice applies to files applications/snowpack/getopt.* and getopt_long.*. All other
 * files in this product are covered by the <a href="http://www.gnu.org/licenses/lgpl.txt">LGPL version 3</a> or above,
 * or <a href="http://www.gnu.org/licenses/gpl.txt">GPL version 3</a> or above unless otherwise specified.
 *
 * Copyright (c) 1987, 1993, 1994
 *	The Regents of the University of California.  All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions
 * are met:
 * 1. Redistributions of source code must retain the above copyright
 *    notice, this list of conditions and the following disclaimer.
 * 2. Redistributions in binary form must reproduce the above copyright
 *    notice, this list of conditions and the following disclaimer in the
 *    documentation and/or other materials provided with the distribution.
 * 3. All advertising materials mentioning features or use of this software
 *    must display the following acknowledgement:
 *	This product includes software developed by the University of
 *	California, Berkeley and its contributors.
 * 4. Neither the name of the University nor the names of its contributors
 *    may be used to endorse or promote products derived from this software
 *    without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE REGENTS AND CONTRIBUTORS ``AS IS'' AND
 * ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED.  IN NO EVENT SHALL THE REGENTS OR CONTRIBUTORS BE LIABLE
 * FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
 * DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS
 * OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
 * HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
 * LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY
 * OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF
 * SUCH DAMAGE.
 *
 */

/**
 * @page requirements Data requirements
 * %Snowpack performs physical modeling of the various processes taking place between the soil, snow cover and atmosphere in order to
 * simulate the evolution of the snow cover based on meteorological input data. It requires the following meteorological parameters:
 * - air temperature (TA)
 * - relative humidity (RH)
 * - wind speed (VW)
 * - incoming short wave radiation (ISWR) or reflected short wave radiation (RSWR)
 * - incoming long wave radiation (ILWR) OR surface temperature (TSS)
 * - snow height (HS) OR precipitation (HNW)
 * - ground temperature (TSG, if available)
 * - snow temperatures at various depths (TS1, TS2, etc if available and only for comparisons, see section \ref SnowSoilTemperatures)
 *
 * These parameters \b must be available at least at a hourly time step.
 *
 * @section data_recomendations Data recommendations
 * In case you have all four irradiative components under ventilated and heated conditions,  the cleanest approach in terms of energy flux calculations seems to be by using:
 * @code
 * SW_MODE = 2
 * INCOMING_LONGWAVE = 1
 * CHANGE_BC = false ;ie Neumann throughout
 * @endcode
 *
 * In case you only have reflected shortwave and snow surface temperature, using Dirichlet boundary condition would be recommended:
 * @code
 * SW_MODE = 1
 * INCOMING_LONGWAVE = 0
 * MEAS_TSS = true
 * CHANGE_BC = true
 * @endcode
 * For energy balance interpretation the change of internal energy is for that case better than the sum of fluxes.
 *
 * @section data_preparation Data preparation
 * In order to help %Snowpack handle the (sometimes broken) data sets to be used in a simulation, the <a href="https://slfsmm.indefero.net/p/meteoio">MeteoIO library</a> is used.
 * This enables %Snowpack to get data from a variety of sources (several input file formats, connection to a database, connection to a web service) and to
 * pre-process real-world data, by filtering the data on the fly and by resampling the data on the fly. Please read the MeteoIO documentation to learn about
 * the supported file formats, the available filters and resampling/re-accumulation strategies.
 *
 * @section data_checks Data checks
 * Please keep in mind that any inaccuracy on the input parameters will have an impact on
 * the quality of the simulation. Since the modelling is physically-based, manually re-constructing missing data can lead to completely wrong results if not
 * done with great care (the model performing various checks on the physical consistency of the input data, it \b will exclude data points that are not consistent
 * with the other parameters. For example, precipitation occuring simultaneously with quite dry air will be refused).
 *
 * \image html clear_sky.png "Data consistency check"
 * \image latex clear_sky.eps "Data consistency check" width=0.9\textwidth
 * For example, the figure above allows to check the following points:
 * - the (solid) precipitation are synchronized with the major snow height increase - this is consistent;
 * - the precipitation happen during hight relative humidity periods - this is consistent;
 * - during times of precipitation, the air and surface temperatures remain close - this means that the sky is cloudy, this is consistent;
 * - a few periods of low wind speed coincide with high relative humidity, which could lead to surface hoar formation - look for it in the simulation outputs;
 * - early in the season, two periods of high wind speed could have lead to drifting snow. The first one occurs during a large snow fall and therefore might be hidden in the data while the second period that experiences a strong snow height decrease could also be driven by rapidly increasing air temperatures - the precipitation might show a large undercatch while the snow might have been wind pressed;
 * - late in the data set, the snow height measurements fail for an extended period of time at a time of high wind speed - some snow drift might have gone unnoticed.
 *
 * When using spurious data or when faced with a bad behaving simulation, one should <b>first look at the consistency of the input data</b>. The vast majority of the problems can be traced back to some low quality data (either for sensor issues or spurious data manipulation at some stage).
 *
 * @section SnowSoilTemperatures Snow and/or soil temperatures
 * Up to 5 snow and/or soil temperatures, either measured or modelled or both can be monitored.
 * Measured temperatures are read in from the input file. If you use the smet format, do not forget to properly
 * label the columns as TS1, TS2, TS3, etc. If you use the snio format, refer to the documentation.
 * User defined positions (m) should be provided in the SnowpackAdvanced section of the \em "io.ini" file,
 *   for example, FIXED_POSITIONS = "0.25 0.50 -0.10":
 *   - positive values refer to heigths measured from the ground surface (snow only)
 *   - negative values refer to depths measured from either the ground surface or the snow surface in case no soil
 *       layers are present
 *   - A sensor must at least be covered by MIN_DEPTH_SUBSURF (m) snow for its temperature to be output.
 *       This parameter can be set in the SnowpackAdvanced section of the io.ini file.
 */

/**
 * @page input_formats File formats
 * Several kind of information need to be given to %Snowpack for a simulation:
 * -# the description of the place where the snow pack has to be simulated: latitutde, longitude, elevation, slope, ...
 * -# the time series of the various meteorological parameters
 * -# the intial state of the various soil and snow layers
 *
 * Very often, 1) and 2) are provided together. But this depends ultimately on the file format that is used ot provide such data (SMET, INP, etc). These two points are
 * handled by <a href="https://slfsmm.indefero.net/p/meteoio">MeteoIO</a>, so please look at its documentation.
 *
 * The intial state of the layers is now given in the SMET ascii file format. This format is described in MeteoIO's documentation of the SMET plugin.
 * Here, two files will be needed: one to contain the layers information and one that contains the temporal data relevant for hazard evaluation (this
 * one is optional and will be automatically created if not provided when starting the simulation).
 *
 * @section layers_data Layers data
 * The snow/soil layers file has the following structure:
 * @code
 * SMET 1.1 ASCII
 * [HEADER]
 * station_id       = DAV2
 * station_name     = Davos:Baerentaelli
 * latitude         = 46.701
 * longitude        = 9.82
 * altitude         = 2560
 * nodata           = -999
 * tz               = 1
 * source           = WSL-Institute for Snow and Avalanche Research SLF; CFierz, 2011-10
 * ProfileDate      = 2009-10-01T00:00
 * HS_Last          = 0.0000
 * SlopeAngle       = 38.0
 * SlopeAzi         = 0.0
 * nSoilLayerData   = 0
 * nSnowLayerData   = 0
 * SoilAlbedo       = 0.20
 * BareSoil_z0      = 0.200
 * CanopyHeight     = 0.00
 * CanopyLeafAreaIndex     = 0.00
 * CanopyDirectThroughfall = 1.00
 * WindScalingFactor       = 1.19
 * ErosionLevel     = 0
 * TimeCountDeltaHS = 0.000000
 * fields           = timestamp Layer_Thick  T  Vol_Frac_I  Vol_Frac_W  Vol_Frac_V  Vol_Frac_S Rho_S Conduc_S HeatCapac_S  rg  rb  dd  sp  mk mass_hoar ne CDot metamo
 * [DATA]
 * 2009-09-19T02:30 0.003399 273.15 0.579671 0.068490 0.351839 0.000000 0.0 0.0 0.0 1.432384 1.028390 0.000000 1.000000 22
 * @endcode
 *
 * @section hazard_data Hazard data
 * The hazards file has the following structure:
 * @code
 * SMET 1.1 ASCII
 * [HEADER]
 * station_id       = DAV2
 * station_name     = Davos:Baerentaelli
 * latitude         = 46.701
 * longitude        = 9.82
 * altitude         = 2560
 * nodata           = -999
 * tz               = 1
 * ProfileDate      = 2012-06-11T17:30
 * fields           = timestamp SurfaceHoarIndex DriftIndex ThreeHourNewSnow TwentyFourHourNewSnow
 * [DATA]
 * 2010-06-08T18:00       -999       -999   0.000000   0.000000
 * 2010-06-08T18:30       -999       -999   0.000000   0.000000
 * 2010-06-08T19:00       -999       -999   0.000000   0.000000
 * 2010-06-08T19:30       -999       -999   0.000000   0.000000
 * 2010-06-08T20:00       -999       -999   0.000000   0.000000
 * 2010-06-08T20:30       -999       -999   0.000000   0.000000
 * 2010-06-08T21:00       -999       -999   0.000000   0.000000
 * ...
 * 2010-06-11T17:30       -999       -999   0.000000   0.000000
 * @endcode
 * As can be seen in this example, the various indices as well as the snow statistics are given every half an hour in reverse chronological order until
 * the profile date.
 */

/**
 * @page output_formats File formats
 * %Snowpack creates various output files:
 * - the current state of its soil and snow layers in <i>".sno"</i> files, which structure is described in \subpage input_formats "input formats";
 * - the current state of its hazard relevant data in <i>".haz"</i> files, which structure is described in \subpage input_formats "input formats";
 * - a time serie of snow profile in <i>".pro"</i> files;
 * - a time serie of the meteorological data as used in the model in <i>".met"</i> files.
 *
 * @section Profiles_data Profiles data
 * The time resolved snow profiles are stored in <i>".pro"</i> files structured as following:
 * @code
 * [STATION_PARAMETERS]
 * StationName      = Davos:Baerentaelli
 * Latitude         = 46.701
 * Longitude        = 9.82
 * Altitude         = 2560
 * SlopeAngle= 0.00
 * SlopeAzi= 0.00
 *
 * [HEADER]
 * #2012-06-11T16:37, Snowpack DEFAULT version 20120611.193 run by "bavay" (research mode)
 * 0500,Date
 * 0501,nElems,height [> 0: top, < 0: bottom of elem.] (cm)
 * 0502,nElems,element density (kg m-3)
 * 0503,nElems,element temperature (degC)
 * 0506,nElems,liquid water content by volume (%)
 * 0508,nElems,dendricity (1)
 * 0509,nElems,sphericity (1)
 * 0510,nElems,coordination number (1)
 * 0511,nElems,bond size (mm)
 * 0512,nElems,grain size (mm)
 * 0513,nElems,grain type (Swiss Code F1F2F3)
 * 0515,nElems,ice volume fraction (%)
 * 0516,nElems,air volume fraction (%)
 * 0517,nElems,stress in (kPa)
 * 0518,nElems,viscosity (GPa s)
 * 0519,nElems,soil volume fraction (%)
 * 0520,nElems,temperature gradient (K m-1)
 * 0521,nElems,thermal conductivity (W K-1 m-1)
 * 0522,nElems,absorbed shortwave radiation (W m-2)
 * 0523,nElems,viscous deformation rate (1.e-6 s-1)
 * 0530,nElems,position (cm) and minimum stability indices:
 *           profile type, stability class, z_Sdef, Sdef, z_Sn38, Sn38, z_Sk38, Sk38
 * 0531,nElems,deformation rate stability index Sdef
 * 0532,nElems,natural stability index Sn38
 * 0533,nElems,stability index Sk38
 * 0534,nElems,hand hardness either (N) or index steps (1)
 * 0535,nElems,optical equivalent grain size (mm)
 * 0601,nElems,snow shear strength (kPa)
 * 0602,nElems,grain size difference (mm)
 * 0603,nElems,hardness difference (1)
 * 0604,nElems,ssi
 * 0605,nElems,inverse texture index ITI (Mg m-4)
 *
 * [DATA]
 * @endcode
 * The each data line starts with a code as described in the header followed by the number of elements (except for the date line) and
 * for each element, the value of the matching parameter. For example, the lines:
 * @code
 * 0500,10.12.1995 12:30
 * 0501,31,27.21,29.07,30.62,31.57,33.30,35.25,37.46,39.82,40.92,42.86,44.22,45.74,47.41,49.15,50.63,52.46,54.58
 * 0502,17,277.7,274.2,268.6,267.0,258.4,248.4,233.5,218.1,207.8,225.1,185.9,176.0,162.5,155.0,127.7,122.7,114.4
 * @endcode
 * provide the date and time (line starting with 0500), then the elements heights for each of the 17 elements (line starting with 0501) and the elements densities (line starting with 0502).
 *
 * @section Met_data Meteorological data
 * The time series of meteorological data as used by the model are stored in <i>".met"</i> files structured as following:
 * @code
 * [STATION_PARAMETERS]
 * StationName= Weissfluhjoch:StudyPlot_MST
 * Latitude= 46.83
 * Longitude= 9.81
 * Altitude= 2540
 * SlopeAngle= 0.00
 * SlopeAzi= 0.00
 * DepthTemp= 0
 *
 * [HEADER]
 * #2012-06-11T16:37, Snowpack DEFAULT version 20120611.193 run by "bavay" (research mode)
 * ,,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56,57,58,59,60,61,62,63,64,65,66,67,68,69,70,71,72,73,74,75,76,77,78,79,80,81,82,83,84,85,86,87,88,89,90,91,92,93,94,95,96,97,98,99,100
 * ID,Date,Sensible heat,Latent heat,Outgoing longwave radiation,Incoming longwave radiation,Net absorbed longwave radiation,Reflected shortwave radiation,Incoming shortwave radiation,Net absorbed shortwave radiation,Modelled surface albedo,Air temperature,Modeled surface temperature,Measured surface temperature,Temperature at bottom of snow or soil pack,Heat flux at bottom of snow or soil pack,Ground surface temperature,Heat flux at ground surface,Heat advected to the surface by liquid precipitation,Global solar radiation (horizontal),Global solar radiation on slope,Direct solar radiation on slope,Diffuse solar radiation on slope,Measured surface albedo,Relative humidity,Wind speed,Max wind speed at snow station or wind speed at ridge station,Wind direction at snow station,Precipitation rate at surface (solid only),Modelled snow depth (vertical),Enforced snow depth (vertical),Surface hoar size,24h Drift index (vertical),Height of new snow HN (24h vertical),3d sum of daily height of new snow (vertical),Total snowpack mass,Eroded mass,Rain rate,Surface runoff (without soil infiltration),Sublimation,Evaporation,Temperature 1 (modelled),Temperature 1 (measured),Temperature 2 (modelled),Temperature 2 (measured),Temperature 3 (modelled),Temperature 3 (measured),Temperature 4 (modelled),Temperature 4 (measured),Temperature 5 (modelled),Temperature 5 (measured),Measured snow depth HS or Solute load at soil surface,SWE (of snowpack),Liquid Water Content (of snowpack),Profile type,Stability class,z_Sdef,Deformation rate stability index Sdef,z_Sn38,Natural stability index Sn38,z_Sk38,Skier stability index Sk38,z_SSI,Structural Stability index SSI,z_S5,Stability index S5,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,Soil runoff,Internal energy change,Surface input (sum fluxes),Measured new snow density,Modeled new snow density,Crust thickness (S-slope),Measured sensible heat,Measured latent heat
 * ,,W m-2,W m-2,W m-2,W m-2,W m-2,W m-2,W m-2,W m-2,1,degC,degC,degC,degC,W m-2,degC,W m-2,W m-2,W m-2,W m-2,W m-2,W m-2,1,%,m s-1,m s-1,deg,kg m-2 h-1,cm,cm,mm,cm,cm,cm,kg m-2,kg m-2 h-1,kg m-2 h-1,kg m-2,kg m-2,kg m-2,degC,degC,degC,degC,degC,degC,degC,degC,degC,degC,cm or kg m-2,kg m-2,kg m-2,-,-,cm,1,cm,1,cm,1,cm,1,cm,1,,,,,,,,,,,,,,,,,,,,,,,,,,,,,kg m-2,kJ m-2,kJ m-2,kg m-3,kg m-3,cm,W m-2,W m-2
 *
 * [DATA]
 * 0203,01.11.1995 00:30,0.795426,-4.160588,308.899297,293.706000,-15.193297,0.000000,0.000000,0.000000,0.090000,0.000000,-0.100000,0.200000,-0.100000,-999.000000,-0.100000,-999.000000,0.000000,0.000000,0.000000,0.000000,0.000000,-999.000000,95.800000,0.800000,0.800000,278.200000,0.000000,0.00,0.00,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,,,,,,,,,,,0.00,0.000000,0.000000,-1,-1,0.0,6.00,0.0,6.00,0.0,6.00,0.0,6.00,0.0,0.00,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,-999.000000,-16.702613,-0.0,-151.3,0.000000,,
 * @endcode
 * Data lines start with an id, followed by the date and the other fields, as shown in the header.
 *
 */

/**
 * @page configuration Configuring a simulation
 * The configuration for a given simulation is kept in a <i>".ini"</i> file (see http://en.wikipedia.org/wiki/INI_file). This is an ascii file that contains
 * keys/values structured by sections. This can be easily edited with a simple text editor. More information about the structure of the file and how to generally deal
 * with it can be found in MeteoIO's documentation (section "How to build your io.ini configuration file"). However, it is recommended to use the inishell tool for
 * generating the configuration file for %Snowpack in order to prevent missing important keys, etc
 *
 * @section inishell_config The inishell tool
 * It is highly recommended to use the <a href="https://slfsmm.indefero.net/p/inishell">Inishell</a> tool to generate these ini files
 * in order to reduce editing errors. This tool also allows you to edit an existing file in order to change the configuration.
 * \image html inishell.png "inishell overview"
 * \image latex inishell.eps "inishell overview" width=0.9\textwidth
 * Each ini file section is shown in a separate tab, each key on its own line. Keys that appear in red are mandatory and will trigger a warning box if trying to
 * visualize/save an ini file before providing a value for these keys. A help text describing the key is shown on the right. Once all has been configured,
 * the configuration can be saved to an ini file, ready to be interpreted by %Snowpack.
 *
 * @section advanced_cfg Advanced configuration
 * The configuration files being an ascii format (<a href="https://en.wikipedia.org/wiki/INI_file">INI format</a>), it is possible to manually
 * copy/paste whole sections of such files between simulations in order to run similar simulations without the need to re-type the whole configuration.
 * It is also possible to use inishell in order to find out which keys are available, which values they can take and what they mean and manually type
 * them in an existing ini file.
 *
 * The %Snowpack_advanced section contains settings that previously required to edit the source code and recompile the model. Since these settings
 * deeply transform the operation of the model, please <b>refrain from using them</b> if you are not absolutely sure of what you are doing.
 *
 */

/**
 * @page sngui_config The sngui tool
 * The simulation outputs are saved in \a ".pro" files for the time resolved profiles and \a ".met" files for the meteorological data time series
 * (see section \subpage output_formats "File formats"). These files can be processed with some scripts, relying on GNU plot for generating graphs
 * but are usually viewed with a graphical application: <a href="http://slfsmm.indefero.net/p/sngui/">sngui</a>. This java application can be
 * downloaded after registering on the web site.
 * \image html sngui_overview_small.png "sngui overview"
 * \image latex sngui_overview.eps "sngui overview" width=0.9\textwidth
 *
 */

/**
 * @page libsnowpack_basics Programming with libsnowpack
 * %Snowpack is now distributed as a very simple application that delegates most of the work to a library, libsnowpack. This C++ library can easily be integrated in
 * other models/applications. In order to do so, the following header have to be included:
 * @code
 * #include <snowpack/libsnowpack.h>
 * #include <meteoio/MeteoIO.h>
 * @endcode
 * Usually, MeteoIO is used to get the meteorological data and the meteoio meteo data class (mio::MeteoData) is converted to %Snowpack meteo data class (CurrentMeteo).
 * The %Snowpack specific configuration options are stored in a SnowpackConfig class that is passed to the various other objects.
 *
 * The main computation is performed by the Snowpack class that needs the following data for its Snowpack::runSnowpackModel call:
 * - the meteo data for the current timestamp as a CurrentMeteo object;
 * - the information about the specific location where the snowpack is simulated contained in a SnowStation object;
 * - the boundary conditions in a BoundCond object;
 * - the surface fluxes in a SurfaceFluxes object.
 *
 * In order to initialize some of these objects from data stored in files, a helper class has been designed: SnowpackIO. Of interest are the following calls:
 * SnowpackIO::readSnowCover, SnowpackIO::writeSnowCover, SnowpackIO::writeProfile, SnowpackIO::writeTimeSeries.
 *
 * In order to compute hazard relevant data, the Hazard class has been designed. The stability data is computed by the Stability class. Some information has
 * to be exchanged between the SnowpackIO object and the Hazard and/or Stability objects. This is handled by the SN_SNOWSOIL_DATA and ZwischenData classes.
 *
 * In order to ease debugging, these classes redefine the "<<" operator to show in a compact way their relevant content. For example,
 * @code
 * SnowStation station;
 * std::cout << station;
 * @endcode
 * shows the relevant parameters of "station":
 * @code
 * <SnowStation>
 * <station>
 * <Coords>
 * Altitude        -999
 * Lat/Long        (-999°0'0.000000" , -999°0'0.000000")
 * Lat/Long        (-999 , -999)
 * X/Y_coords      (-999 , -999)
 * I/J_indices     (-999 , -999)
 * Projection      NULL NULL
 * EPSG            -1
 * </Coords>
 * ID: MST96
 * Name: MST96
 * Slope: -999 bearing: -999
 * </station>
 * 0 element(s) and 0 node(s). Soil=false canopy=false
 * Soil:   SoilNode=0 depth=0 BareSoil_z0=0
 * Albedo: mAlbedo=0 cAlbedo=0 SoilAlb=0
 * Snow:   Measured HS=0 Calculated HS=0 New snow=0 of density=0
 * Energy: ColdContent=0 dIntEnergy=0 SubSurfaceMelt=x SubSurfaceFrze=x
 * Snowdrift:      sector=0 windward=0 ErosionLevel=0 ErosionMass=0
 * Stability:      S_d(0)=0 S_n(0)=0 S_s(0)=0 S_1=0 S_2=0 S_4(0)=0 S_5(0)=0
 * Kt= 0
 * </SnowStation>
 * @endcode
 *
 */

/**
 * @page coding_style Coding style
 * @section coding_sty Recomended coding style
 * The recommended coding style for MeteoIO is the <A HREF="http://www.kernel.org/doc/Documentation/CodingStyle">Kernel coding style</A> with a few exceptions:
 * - we don't enforce strict 80 characters line width. try to remain reasonable, but don't necessarily cut everything off at 80 characters
 * - try to intelligently use spaces to visually group elements of a complex formula. If the formula can be split into meaningful elements,
 *   please do it (using some "const double element = " constructs).
 * - try to properly qualify variables: for example, if a variable will not be changed, will never be negative and always integer,
 *   then use "const unsigned int". When some specific types are used for some standard library calls, try to properly use these types (for example, "size_t")
 * - use C++ method naming convention: a method name starts with lowercase letter and each individual word in a name starts capitalized.
 *   Usually, no underscores are used in a method. For example, a method that would return the lapse rate contained in an object would be named "getLapseRate()"
 * - qualify variables and parameters with "const" when appropriate (see <A HREF="http://jriddell.org/const-in-cpp.html">const-in-cpp</A>).
 *
 * A few important points to emphasize (from the <A HREF="http://www.kernel.org/doc/Documentation/CodingStyle">Kernel coding style</A>):
 * - Functions should be short and sweet, and do just one thing.  They should fit on one or two screenfuls of text, and do one thing and do that well.
 * - If you have a complex function, and you suspect that a less-than-gifted first-year high-school student might not even understand
 *   what the function is all about, you should adhere to the maximum limits all the more closely.  Use helper functions with descriptive names.
 * - Comments are good, but there is also a danger of over-commenting.  NEVER try to explain HOW your code works in a comment:
 *   it's much better to write the code so that the _working_ is obvious, and it's a waste of time to explain badly written code.
 *
 * @section code_indentation Indentation
 * Since every user has his/her own preference for the ideal indentation width, please use <A HREF="http://www.emacswiki.org/emacs/SmartTabs">"smart tabs"</A>.
 * That practically means:
 * - indent with tabs
 * - align with spaces
 *
 * This way, each developer can set his/her indentation size as he/she wishes without forcing his/her choice to others...
 *
 * @section containers Memory management and Containers
 * Please do NOT manage memory manually but use <A HREF="https://secure.wikimedia.org/wikipedia/en/wiki/Standard_Template_Library">Standard Template Library (STL)
 * </A> <A HREF="http://www.cplusplus.com/reference/stl/">containers</A> instead.
 * This dramatically reduces memory errors (ie: <A HREF="https://secure.wikimedia.org/wikipedia/en/wiki/Segmentation_fault">segfaults</A>), often
 * offers more performance and provides you with lots of <A HREF="http://www.cplusplus.com/reference/algorithm/">associated algorithms</A>
 * (like sorting, search, filling, etc).
 *
 * When you need your own data class, please design it based on these STL containers (like grid2DObject is based on std::vector). Basically, this means
 * that you will replace mallocs and arrays by vectors (for 1d, 2d, 3d grids), maps (for multiple key/value pairs), lists (for unordered table), etc
 *
 * @section exceptions_handinling Exceptions handling
 * The recommended C++ usage should be followed: <b>"throw by value, catch by reference"</b> (as specified in <i>C++ Coding Standards: 101 Rules, Guidelines,
 * and Best Practices</i>, Herb Sutter, Andrei Alexandrescu, 2004, Addison-Wesley Professional). Moreover, we should consider catching by
 * <b>const reference</b> and not even declaring a variable if not doing anything with it: something like `catch(const IOException&)` would often be enough.
 *
 */

/**
 * @page adding_extra_models Adding extra models
 * Various processes can already be simulated using different models as configured by the user. This result is achieved by providing a specific model of
 * the process of interest, together with the proper entry in a std::map container that links a model keyword with its implementation. In order to look at the
 * required steps, we will take as an example the hand hardness implementation in the Stability class. Please keep in mind that when adding a new model to
 * a process that already has multiple available choices, only the first and the third steps are required, the other one being already done.
 *
 * @section model_implementation Model implementation
 * A method has to be implemented in the class with the same prototype as the original method. In our example, the original method (setHandHardnessMONTI)
 * has the following prototype:
 * @code
 * double setHandHardnessMONTI(const ElementData& Edata);
 * @endcode
 * so any alternative implementation must use the same prototype. If some parameters would be ignored by some implementation, simply comment out the unused variable:
 * @code
 * double my_method(const double& param1, const double /*unused_param*/);
 * @endcode
 *
 * @section function_pointer Function pointer typedef
 * All these methods sharing the same prototype, a generic function pointer type (actually, a method pointer) can be defined:
 * @code
 * typedef double (Stability::*StabMemFn)(const ElementData&);
 * @endcode
 *
 * @section model_map Model map
 * Once an alternative implementation has been written (and properly declared in the header file), it must be "registered" in the model map. In our exmaple, this map
 * is defined in the header %file:
 * @code
 * static std::map<std::string, StabMemFn> mapHandHardness;
 * @endcode
 * and statically filled in the initStaticData() method as following:
 * @code
 * const bool Stability::__init = Stability::initStaticData();
 * bool Stability::initStaticData() {
 * 	mapHandHardness["MONTI"]    = &Stability::setHandHardnessMONTI;
 * 	mapHandHardness["BELLAIRE"]  = &Stability::setHandHardnessBELLAIRE;
 * 	mapHandHardness["ASARC"]    = &Stability::setHandHardnessASARC;
 * 	return true;
 * }
 * @endcode
 * This way of fillinf the map ensures that it will be initialized only once and for all, making it consistent and efficient.
 *
 * @section model_user_choice User model configuration
 * The user selection of model must be connected with the proper model implementation. The user selects the model he wants through a key in his
 * configuration file. We need to read this key and activate the proper implementation, knowing that the proper key <-> implementation matching is done
 * through the map:
 * @code
 * cfg.getValue("HARDNESS_MODEL", "SnowpackAdvanced", hardness_model);
 * map<string, StabMemFn>::const_iterator it1 = mapHandHardness.find(hardness_model);
 * if (it1 == mapHandHardness.end()) throw InvalidArgumentException("Unknown hardness model: "+hardness_model, AT);
 * @endcode
 * This means that in the section "SnowpackAdvanced" of his ini file, the key "HARDNESS_MODEL" must contain one of the strings given in the mapHandHardness
 * above (ie. either "DEFAULT" or "ASARC" or "MONTI").
 *
 * @section calling_model Model call
 * Finally, the process model has to be called where needed. A helper macro can be defined as
 * @code
 * #define CALL_MEMBER_FN(object,ptrToMember)  ((object).*(ptrToMember))
 * @endcode
 * and in the code, each time the hand hardness has to be computed, the call becomes:
 * @code
 * hardness = CALL_MEMBER_FN(*this, mapHandHardness[hardness_model])(EMS[e]);
 * @endcode
 *
 */

#endif







