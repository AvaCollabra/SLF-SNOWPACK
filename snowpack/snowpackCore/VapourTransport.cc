/*
 *  SNOWPACK stand-alone
 *
 *  Copyright WSL Institute for Snow and Avalanche Research SLF, DAVOS, SWITZERLAND
*/
/*  This file is part of Snowpack.
	Snowpack is free software: you can redistribute it and/or modify
	it under the terms of the GNU General Public License as published by
	the Free Software Foundation, either version 3 of the License, or
	(at your option) any later version.

	Snowpack is distributed in the hope that it will be useful,
	but WITHOUT ANY WARRANTY; without even the implied warranty of
	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
	GNU General Public License for more details.

	You should have received a copy of the GNU General Public License
	along with Snowpack.  If not, see <http://www.gnu.org/licenses/>.
*/


#include <snowpack/snowpackCore/VapourTransport.h>
#include <snowpack/vanGenuchten.h>
#include <snowpack/snowpackCore/Snowpack.h>

#include <assert.h>
#include <sstream>
#include <errno.h>

//Eigen
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <Eigen/IterativeLinearSolvers>
#include <Eigen/SparseQR>
#include <Eigen/SparseCholesky>
#include <Eigen/SparseLU>
#include <Eigen/Core>

#include "omp.h"


typedef Eigen::Triplet<double> Trip;

#ifdef CLAPACK
	// Matching C data types with FORTRAN data types (taken from f2c.h):
	typedef long int integer;
	typedef double doublereal;

	// Declare the function interfaces with the LAPACK library (taken from clapack.h):
	extern "C" {
		/* Subroutine */ int dgesvd_(char *jobu, char *jobvt, integer *m, integer *n,
		doublereal *a, integer *lda, doublereal *s, doublereal *u, integer *
		ldu, doublereal *vt, integer *ldvt, doublereal *work, integer *lwork,
		integer *info);

		/* Subroutine */ int dgesdd_(char *jobz, integer *m, integer *n, doublereal *
		a, integer *lda, doublereal *s, doublereal *u, integer *ldu,
		doublereal *vt, integer *ldvt, doublereal *work, integer *lwork,
		integer *iwork, integer *info);

		/* Subroutine */ int dgtsv_(integer *n, integer *nrhs, doublereal *dl,
		doublereal *d__, doublereal *du, doublereal *b, integer *ldb, integer
		*info);
	}
#endif

using namespace mio;
using namespace std;
using namespace Eigen;

VapourTransport::VapourTransport(const SnowpackConfig& cfg)
               : WaterTransport(cfg), RichardsEquationSolver1d(cfg, false), variant(),
                 iwatertransportmodel_snow(BUCKET), iwatertransportmodel_soil(BUCKET), watertransportmodel_snow("BUCKET"), watertransportmodel_soil("BUCKET"),
                 sn_dt(IOUtils::nodata),timeStep(IOUtils::nodata),waterVaporTransport_timeStep(IOUtils::nodata),
                 hoar_thresh_rh(IOUtils::nodata), hoar_thresh_vw(IOUtils::nodata), hoar_thresh_ta(IOUtils::nodata),
                 /*hoar_density_buried(IOUtils::nodata), hoar_density_surf(IOUtils::nodata), hoar_min_size_buried(IOUtils::nodata),
                 minimum_l_element(IOUtils::nodata),*/ useSoilLayers(false), water_layer(false), enable_vapour_transport(false),enable_vapour_transport_soil(false),
                 diffusionScalingFactor_(1.0), height_of_meteo_values(0.), adjust_height_of_meteo_values(true), f(1.0),waterVaporTransport_timeStepAdjust(false)
{
	cfg.getValue("VARIANT", "SnowpackAdvanced", variant);

	// Defines whether soil layers are used
	cfg.getValue("SNP_SOIL", "Snowpack", useSoilLayers);

	//To build a thin top rain-water layer over a thin top ice layer, rocks, roads etc.
	cfg.getValue("WATER_LAYER", "SnowpackAdvanced", water_layer);

	/**
	 * @brief No surface hoar will form for rH above threshold (1)
	 * - Original calibration with the 98/99 data set: 0.9
	 * - r141: HOAR_THRESH_RH set to 0.9
	 * - r719: HOAR_THRESH_RH set to 0.97
	 */
	cfg.getValue("HOAR_THRESH_RH", "SnowpackAdvanced", hoar_thresh_rh);

	/**
	 * @brief No surface hoar will form at wind speeds above threshold (m s-1)
	 * - Original calibration with the 98/99 data set: 3.5
	 * - r141: HOAR_THRESH_VW set to 3.0
	 * - r242: HOAR_THRESH_VW set to 3.5
	 */
	cfg.getValue("HOAR_THRESH_VW", "SnowpackAdvanced", hoar_thresh_vw);

	/**
	 * @brief No surface hoar will form at air temperatures above threshold (m s-1)
	 * - Originaly, using THRESH_RAIN
	 * - r787: HOAR_THRESH_TA set to 1.2
	 */
	cfg.getValue("HOAR_THRESH_TA", "SnowpackAdvanced", hoar_thresh_ta);

	//Calculation time step in seconds as derived from CALCULATION_STEP_LENGTH
	const double calculation_step_length = cfg.get("CALCULATION_STEP_LENGTH", "Snowpack");
	sn_dt = M_TO_S(calculation_step_length);

	//Enable vapour transport
	cfg.getValue("ENABLE_VAPOUR_TRANSPORT", "SnowpackAdvanced", enable_vapour_transport);
	if(enable_vapour_transport)
	{
		cfg.getValue("ENABLE_VAPOUR_TRANSPORT_SOIL", "SnowpackAdvanced", enable_vapour_transport_soil);
		
		// the water vapor subtime step
		cfg.getValue("WATER_VAPOR_TRANSPORT_TIMESTEP", "SnowpackAdvanced", waterVaporTransport_timeStep);
		waterVaporTransport_timeStep=std::min(sn_dt,waterVaporTransport_timeStep);
		
		// the water vapor transport scheme, f=1 fully implicit, f=0.5 Crank-Nicolson
		cfg.getValue("WATER_VAPOR_TRANSPORT_IMPLICIT_FACTOR", "SnowpackAdvanced", f);

		//cfg.getValue("WATER_VAPOR_TRANSPORT_TIMESTEP_ADJUST", "SnowpackAdvanced", waterVaporTransport_timeStepAdjust);
		if(f<1.0) waterVaporTransport_timeStepAdjust=true;
	}

	//reading diffusionScalingFactor_
	//cfg.getValue("DIFFUSION_SCALING_FACTOR", "SnowpackAdvanced", diffusionScalingFactor_);

	//Water transport model snow
	cfg.getValue("WATERTRANSPORTMODEL_SNOW", "SnowpackAdvanced", watertransportmodel_snow);
	iwatertransportmodel_snow=UNDEFINED;
	if (watertransportmodel_snow=="BUCKET") {
		iwatertransportmodel_snow=BUCKET;
	} else if (watertransportmodel_snow=="NIED") {
		iwatertransportmodel_snow=NIED;
	} else if (watertransportmodel_snow=="RICHARDSEQUATION") {
		iwatertransportmodel_snow=RICHARDSEQUATION;
	}

	//Water transport model soil
	cfg.getValue("WATERTRANSPORTMODEL_SOIL", "SnowpackAdvanced", watertransportmodel_soil);
	iwatertransportmodel_soil=UNDEFINED;
	if (watertransportmodel_soil=="BUCKET") {
		iwatertransportmodel_soil=BUCKET;
	} else if (watertransportmodel_soil=="NIED") {
		iwatertransportmodel_soil=NIED;
	} else if (watertransportmodel_soil=="RICHARDSEQUATION") {
		iwatertransportmodel_soil=RICHARDSEQUATION;
	}

	cfg.getValue("HEIGHT_OF_METEO_VALUES", "Snowpack", height_of_meteo_values);
	cfg.getValue("ADJUST_HEIGHT_OF_METEO_VALUES", "SnowpackAdvanced", adjust_height_of_meteo_values);

}

void VapourTransport::compTransportMass(const CurrentMeteo& Mdata, double& ql,
                                       SnowStation& Xdata, SurfaceFluxes& Sdata, const double& surfaceVaporPressure)
{

	// First, consider no soil with no snow on the ground
	if (!useSoilLayers && Xdata.getNumberOfNodes() == Xdata.SoilNode+1) {
		return;
	}

    try {
	    LayerToLayer(Mdata, Xdata, Sdata, ql, surfaceVaporPressure);
	    WaterTransport::adjustDensity(Xdata);
    } catch(const exception&)
    {
	    prn_msg( __FILE__, __LINE__, "err", Mdata.date, "Error in transportVapourMass()");
	    throw;
    }
}

void VapourTransport::LayerToLayer(const CurrentMeteo& Mdata, SnowStation& Xdata, SurfaceFluxes& Sdata, double& ql, const double& surfaceVaporPressure)
 {

	if (enable_vapour_transport)
	{
		                                                                   
		//if(!enable_vapour_transport_soil || Xdata.SoilNode<Xdata.getNumberOfNodes()-1)
		//{
			compSurfaceSublimation(Mdata, ql, Xdata, Sdata);
			ql=0;
		//}

		if(!enable_vapour_transport_soil && Xdata.SoilNode==Xdata.getNumberOfNodes()-1) return;
		 
		const size_t nN = Xdata.getNumberOfNodes();
		size_t nE = nN-1;
		vector<NodeData>& NDS = Xdata.Ndata;
		vector<ElementData>& EMS = Xdata.Edata;
		size_t e = nE;
		std::vector<double> deltaM(nE, 0.);//Calculate the limited layer mass change
		std::vector<double> totalMassChange(nE, 0.);// store the total mass change
		std::vector<double> vaporFluxDiv(nE, 0.);// store the vapor flux divergence
		std::vector<double> oldVaporDenEl(nE, 0.);// old water vapor density for element
		std::vector<double> oldVaporDenNode(nN, 0.);// old water vapor density for node

		std::vector<double> factor_(nE, 1.);// this is for source term in vapor transport equation
		for(size_t i=0; i<Xdata.SoilNode; i++)
		{
			factor_[i]=0.;
			//EMS[i].VG.SetVGParamsSoil();
		}
		
		std::vector<double> D_(nE, 0.);
		for (size_t i=0; i<=nE-1; i++)
		{
			double theta_air = std::max(EMS[i].theta[AIR],0.0);
			double tortuosity = pow(theta_air,7./3.)/pow(1-EMS[i].theta[SOIL], 2.);    
			double D_vapSoil = tortuosity * theta_air * Constants::diffusion_coefficient_in_air; 
			//D_vapSoil = 1.5*theta_air *Constants::diffusion_coefficient_in_air;
			// const double D_vapSoil = Constants::diffusion_coefficient_in_air * pow(aaa,10./3.)/nnn/nnn; // based on jury1983

			// based on Colbeck
			//double Dsnow = Constants::diffusion_coefficient_in_snow;
			
			/*
			// based on Calonne et al. 2014
			double Dsnow = 0.;
			if(EMS[i].Rho < 450.)
			{
				Dsnow= ( (3./2.)*(1.-EMS[i].Rho/Constants::density_ice)-(1./2.) ) * Constants::diffusion_coefficient_in_air;
			}
			else
			{
				Dsnow= ( (3./2.)*(1.-450./Constants::density_ice)-(1./2.) ) * (1-(EMS[i].Rho-450.)/(Constants::density_ice-450.)) * Constants::diffusion_coefficient_in_air;
			}
			*/
			
			// based on A. C. Hansen and W. E. Foslien 2015: https://www.the-cryosphere.net/9/1857/2015/
			//double Dsnow = (EMS[i].theta[ICE]*EMS[i].theta[AIR]+1.)*Constants::diffusion_coefficient_in_air;
			
			
			// based on Foslien (1994)
			double Dsnow = EMS[i].theta[ICE]*theta_air*Constants::diffusion_coefficient_in_air+
						   theta_air*Constants::diffusion_coefficient_in_air/
						   (EMS[i].theta[ICE]*Constants::conductivity_air/Constants::conductivity_ice+
							EMS[i].theta[ICE]*Constants::lh_sublimation*Constants::diffusion_coefficient_in_air*dRhov_dT(EMS[i].Te)/Constants::conductivity_ice+
							theta_air
						   );
			/*
			double Dsnow = EMS[i].theta[ICE]*EMS[i].theta[AIR]*Constants::diffusion_coefficient_in_air+
						   EMS[i].theta[AIR]*Constants::diffusion_coefficient_in_air/
						   (EMS[i].theta[ICE]*Constants::conductivity_air/Constants::conductivity_ice+
							EMS[i].theta[ICE]*Constants::lh_sublimation*Constants::diffusion_coefficient_in_air*dRhov_dT(EMS[i].Te)/Constants::conductivity_ice+
							EMS[i].theta[AIR]
						   );
						   */
						  
			
			D_[i] = factor_[i]*Dsnow + (1.0-factor_[i])*D_vapSoil;
			//if(i==Xdata.SoilNode-1) std::cout << "i/D_[i]/factor_[i]/EMS[i].theta[AIR]/EMS[i].theta[ICE]/EMS[i].theta[SOIL]/D_vapSoil " << i << ' ' << D_[i] << ' ' << factor_[i] << ' ' << EMS[i].theta[AIR] << ' ' << EMS[i].theta[ICE] << ' ' << EMS[i].theta[SOIL] << ' ' << D_vapSoil << "\n";					
		}				
		
		std::vector<double> hm_(nN, 0.); // mass transfer coefficient m s-1
		double water_molecular_mass = 18.0153e-3; // (kg)
		double gaz_constant = 8.31451; // (J mol-1 K-1)
		for(size_t i=0; i<nN; i++)
		{
			double Re,saturationDensity;
			saturationDensity = Atmosphere::waterVaporDensity(NDS[i].T, Atmosphere::vaporSaturationPressure(NDS[i].T));	
			hm_[i] =Constants::density_ice/saturationDensity/9.7e9; // hm_experimental, Pirmin 2012, M_mm=as_all*hm_experimental*(rhov_sat-rhov)

			//Re = 0.0;
			//hm_[i] =1.e-3*(0.566*Re+0.075); //Neumann 2009
			//hm_[i] =Constants::density_ice/saturationDensity/9.7e9; // Pirmin 2012
			//if(i>=Xdata.SoilNode) std::cout << "i/hm " << i << ' ' << pow(saturationDensity/1.88*(Constants::lh_sublimation*water_molecular_mass/gaz_constant/NDS[i].T-1.0)*
			//									   Constants::lh_sublimation/Constants::conductivity_air/NDS[i].T+1.0/Constants::diffusion_coefficient_in_air/1.88,-1.0) << "\n";
			//if(i>=Xdata.SoilNode) std::cout << "i/hm " << i << ' ' << (3.8e-5/Atmosphere::waterVaporDensity(258, Atmosphere::vaporSaturationPressure(258))-1)*
			//										0.00008659*6.0/(0.001*1.3)*(1./(pow(0.001*1.3,2.)*3.14))*(0.001*1.3*3.14)*pow(1/1.88*(Constants::lh_sublimation*water_molecular_mass/gaz_constant/258.-1.0)*
			//								   Constants::lh_sublimation/Constants::conductivity_air/254+1.0/Constants::diffusion_coefficient_in_air/1.88/Atmosphere::waterVaporDensity(258, Atmosphere::vaporSaturationPressure(258)),-1.0) << "\n";
		
			// hm_theoretical,------->M_mm=as_active*hm_theoretical*(rhov_sat-rhov)
			/*
			if(i==0)
				double rwTors_u=pow((EMS[i].theta[WATER]+EMS[i].theta[ICE])/EMS[i].theta[SOIL]+1., 1./3.);				
				hm_[i]=2.0*Constants::diffusion_coefficient_in_air/(0.002*rwTors_u*EMS[i].rg);
			}
			else if(i>0 && i<Xdata.SoilNode)
				double rwTors_u=pow((EMS[i].theta[WATER]+EMS[i].theta[ICE])/EMS[i].theta[SOIL]+1., 1./3.);
				double rwTors_d=pow((EMS[i-1].theta[WATER]+EMS[i-1].theta[ICE])/EMS[i-1].theta[SOIL]+1., 1./3.);
				hm_[i]=2.0*Constants::diffusion_coefficient_in_air/(0.5*0.002*rwTors_d*EMS[i-1].rg+0.5*0.002*rwTors_u*EMS[i].rg);
			}
			else if(i==Xdata.SoilNode && Xdata.SoilNode==nN-1 )
			{
				double rwTors_d=pow((EMS[i-1].theta[WATER]+EMS[i-1].theta[ICE])/EMS[i-1].theta[SOIL]+1., 1./3.);							
				hm_[i]=2.0*Constants::diffusion_coefficient_in_air/(0.002*rwTors_d*EMS[i-1].rg);
				//std::cout << "soil rg: " << EMS[i-1].rg << ' ' << ' ' << rwTors_d << ' ' << hm_[i] << "\n";														
			}
			else if(i==Xdata.SoilNode && Xdata.SoilNode<nN-1 )
				double rwTori_u=pow(EMS[i].theta[WATER]/EMS[i].theta[ICE]+1., 1./3.);
				double rwTors_d=pow((EMS[i-1].theta[WATER]+EMS[i-1].theta[ICE])/EMS[i-1].theta[SOIL]+1., 1./3.);							
				hm_[i]=2.0*Constants::diffusion_coefficient_in_air/(0.5*0.002*rwTors_d*EMS[i-1].rg+0.5*0.001*rwTori_u*EMS[i].ogs);
				//hm_[i]=1.88*Constants::diffusion_coefficient_in_air/(0.001*rwTori_u*EMS[i].ogs);
			}
			else if(i>Xdata.SoilNode && i<nN-1)
			{
				double rwTori_u=pow(EMS[i].theta[WATER]/EMS[i].theta[ICE]+1., 1./3.);
				double rwTori_d=pow(EMS[i-1].theta[WATER]/EMS[i-1].theta[ICE]+1., 1./3.);
				hm_[i]=2.0*Constants::diffusion_coefficient_in_air/(0.5*0.001*rwTori_d*EMS[i-1].ogs+0.5*0.001*rwTori_u*EMS[i].ogs);								
			}
			else //i==nN-1
			{
				double rwTori_d=pow(EMS[i-1].theta[WATER]/EMS[i-1].theta[ICE]+1., 1./3.);
				hm_[i]=2.0*Constants::diffusion_coefficient_in_air/(0.001*rwTori_d*EMS[i-1].ogs);
			}
			*/
		}
		
		std::vector<double> as_(nN, 0.); // the specific surface area m-1				
		for(size_t i=0; i<nN; i++)
		{
			double saturationDensity = Atmosphere::waterVaporDensity(NDS[i].T, Atmosphere::vaporSaturationPressure(NDS[i].T));										

			// as_all---------->M_mm=as_all*hm_experimental*(rhov_sat-rhov)
			if(i==0)
			{
				double rwTors_u=pow((EMS[i].theta[WATER]+EMS[i].theta[ICE])/EMS[i].theta[SOIL]+1., 1./3.);				
				double apparentTheta=EMS[i].theta[SOIL]+EMS[i].theta[ICE]+EMS[i].theta[WATER];
				//if(apparentTheta==EMS[i].theta[SOIL]) apparentTheta=0.; //dry soil
				as_[i]=6.0*apparentTheta/(0.002*rwTors_u*EMS[i].rg);
			}
			else if(i>0 && i<Xdata.SoilNode)
			{
				double rwTors_u=pow((EMS[i].theta[WATER]+EMS[i].theta[ICE])/EMS[i].theta[SOIL]+1., 1./3.);
				double rwTors_d=pow((EMS[i-1].theta[WATER]+EMS[i-1].theta[ICE])/EMS[i-1].theta[SOIL]+1., 1./3.);
				double apparentTheta=0.5*(EMS[i].theta[SOIL]+EMS[i].theta[ICE]+EMS[i].theta[WATER])
									+0.5*(EMS[i-1].theta[SOIL]+EMS[i-1].theta[ICE]+EMS[i-1].theta[WATER]);
				//if(apparentTheta==0.5*(EMS[i].theta[SOIL]+EMS[i-1].theta[SOIL])) apparentTheta=0.; //dry soil
									
				as_[i]=6.0*apparentTheta/(0.5*0.002*rwTors_d*EMS[i-1].rg+0.5*0.002*rwTors_u*EMS[i].rg);
			}
			else if(i==Xdata.SoilNode && Xdata.SoilNode==nN-1 )
			{
				double rwTors_d=pow((EMS[i-1].theta[WATER]+EMS[i-1].theta[ICE])/EMS[i-1].theta[SOIL]+1., 1./3.);
				double apparentTheta=EMS[i-1].theta[SOIL]+EMS[i-1].theta[ICE]+EMS[i-1].theta[WATER];
				//if(apparentTheta==EMS[i-1].theta[SOIL]) apparentTheta=0.; //dry soil				
				as_[i]=6.0*apparentTheta/(0.002*rwTors_d*EMS[i-1].rg);											
				//std::cout << "soil rg: " << EMS[i-1].rg << ' ' << ' ' << rwTors_d << ' ' << hm_[i] << "\n";														
			}
			else if(i==Xdata.SoilNode && Xdata.SoilNode<nN-1 )
			{
				double rwTori_u=pow(EMS[i].theta[WATER]/EMS[i].theta[ICE]+1., 1./3.);
				double rwTors_d=pow((EMS[i-1].theta[WATER]+EMS[i-1].theta[ICE])/EMS[i-1].theta[SOIL]+1., 1./3.);
				double apparentTheta=0.5*(EMS[i].theta[ICE]+EMS[i].theta[WATER])
									+0.5*(EMS[i-1].theta[SOIL]+EMS[i-1].theta[ICE]+EMS[i-1].theta[WATER]);
				as_[i]=6.0*apparentTheta/(0.5*0.002*rwTors_d*EMS[i-1].rg+0.5*0.001*rwTori_u*EMS[i].ogs);											
			}
			else if(i>Xdata.SoilNode && i<nN-1)
			{
				double rwTori_u=pow(EMS[i].theta[WATER]/EMS[i].theta[ICE]+1., 1./3.);
				double rwTori_d=pow(EMS[i-1].theta[WATER]/EMS[i-1].theta[ICE]+1., 1./3.);
				double apparentTheta=0.5*(EMS[i].theta[ICE]+EMS[i].theta[WATER])
									+0.5*(EMS[i-1].theta[ICE]+EMS[i-1].theta[WATER]);
				as_[i]=6.0*apparentTheta/(0.5*0.001*rwTori_d*EMS[i-1].ogs+0.5*0.001*rwTori_u*EMS[i].ogs);														
			}
			else //i==nN-1
			{
				double rwTori_d=pow(EMS[i-1].theta[WATER]/EMS[i-1].theta[ICE]+1., 1./3.);
				double apparentTheta=EMS[i-1].theta[ICE]+EMS[i-1].theta[WATER];
				as_[i]=6.0*apparentTheta/(0.001*rwTori_d*EMS[i-1].ogs);																		
			}						

			// as_active------->M_mm=as_active*hm_theoretical*(rhov_sat-rhov)
			/*			
			if(i==0)
			{
				as_[i]=6.0*EMS[i].theta[SOIL]*Constants::density_ice/saturationDensity/9.7e9/2.0/Constants::diffusion_coefficient_in_air;
			}
			else if(i>0 && i<Xdata.SoilNode)
			{
				as_[i]=6.0*(0.5*EMS[i-1].theta[SOIL]+0.5*EMS[i].theta[SOIL])*Constants::density_ice/saturationDensity/9.7e9/2.0/Constants::diffusion_coefficient_in_air;
			}
			else if(i==Xdata.SoilNode && Xdata.SoilNode==nN-1 )
			{
				//as_[i]=6.0*(0.5*EMS[i-1].theta[SOIL]+0.5)*Constants::density_ice/saturationDensity/9.7e9/(0.5*1.88+0.5*2.)/Constants::diffusion_coefficient_in_air;
				as_[i]=6.0*(0.5*EMS[i-1].theta[SOIL])*Constants::density_ice/saturationDensity/9.7e9/2.0/Constants::diffusion_coefficient_in_air;
			}
			else if(i==Xdata.SoilNode && Xdata.SoilNode<nN-1 )
			{
				as_[i]=0.5*(6.0*Constants::density_ice/saturationDensity/9.7e9/1.88/Constants::diffusion_coefficient_in_air)
					  +0.5*(6.0*EMS[i-1].theta[SOIL]*Constants::density_ice/saturationDensity/9.7e9/2.0/Constants::diffusion_coefficient_in_air);								
			}
			else //i==nN-1
			{
				as_[i]=6.0*Constants::density_ice/saturationDensity/9.7e9/1.88/Constants::diffusion_coefficient_in_air;
			}
			*/											
		}
				
		double min_dt =1e30;
		for(size_t i=Xdata.SoilNode; i<nE; i++)
		{
			double saturationDensity=Atmosphere::waterVaporDensity(EMS[i].Te, Atmosphere::vaporSaturationPressure(EMS[i].Te));
			double diffVaporVelocity =std::abs(-D_[i]*(NDS[i+1].rhov-NDS[i].rhov)/EMS[i].L/saturationDensity);
			double dt=EMS[i].L/diffVaporVelocity;
			min_dt = std::min(min_dt,dt);			
		}
		//std::cout << "min_dt/waterVaporTransport_timeStep " << min_dt << ' ' << waterVaporTransport_timeStep << "\n";

		timeStep=(waterVaporTransport_timeStepAdjust) ? std::min(min_dt,waterVaporTransport_timeStep) : waterVaporTransport_timeStep;		
		int nTime= int(sn_dt/timeStep)+1;
		double time = 0;		
		for( size_t l=0; l<=nTime; l++)
		{
			time=time+timeStep;
			if(time>=sn_dt)
			{
				timeStep=sn_dt-(time-timeStep);
				time=sn_dt;
			}

			if (!compDensityProfile(Mdata, Xdata, true, ql, surfaceVaporPressure, hm_, as_, D_,oldVaporDenNode)) break;
			
			//if(time==sn_dt) break;			 				
		//}			

			for (size_t i=0; i<=nE-1; i++)
			{
				double saturationDensity=Atmosphere::waterVaporDensity(EMS[i].Te, Atmosphere::vaporSaturationPressure(EMS[i].Te));				
				double saturationVaporUp = Atmosphere::waterVaporDensity(NDS[i+1].T, Atmosphere::vaporSaturationPressure(NDS[i+1].T));
				double saturationVaporDown = Atmosphere::waterVaporDensity(NDS[i].T, Atmosphere::vaporSaturationPressure(NDS[i].T));				
				double diffRhov_hm_as_Up=(f*NDS[i+1].rhov+(1-f)*oldVaporDenNode[i+1]-saturationVaporUp)*hm_[i+1]*as_[i+1];
				double diffRhov_hm_as_Down=(f*NDS[i].rhov+(1-f)*oldVaporDenNode[i]-saturationVaporDown)*hm_[i]*as_[i];					
				totalMassChange[i] =(0.5*diffRhov_hm_as_Down+0.5*diffRhov_hm_as_Up)*timeStep*EMS[i].L; //total mass change, (kg m-2 )
				//totalMassChange[i]=(0.5*hm_[i+1]*as_[i+1]+0.5*hm_[i]*as_[i])*(EMS[i].rhov-saturationDensity)*timeStep*EMS[i].L;
			}
			
			e = nE;			
			// consider the mass change due to vapour transport in snow/soil
			while (e-- > 0) {
				//const double massPhaseChange = totalMassChange[e]+deltaM[e];
				const double massPhaseChange = totalMassChange[e]+deltaM[e];//JJJJJ

				double dM = 0.;	//mass change induced by vapor flux (kg m-2)

				// Now, the mass change is limited by:
				// - we cannot remove more WATER and ICE than available
				// - we cannot add more WATER and ICE than pore space available
				if ( EMS[e].theta[SOIL] < Constants::eps ) {// there is no soil in element to keep element not to merge
					dM = std::max(  -((EMS[e].theta[WATER] - EMS[e].VG.theta_r * (1. + Constants::eps)) * Constants::density_water * EMS[e].L + (EMS[e].theta[ICE] - Snowpack::min_ice_content) * Constants::density_ice * EMS[e].L)  ,
									std::min(  (EMS[e].theta[AIR] * Constants::density_ice * EMS[e].L), massPhaseChange  )
						 ); // mass change due to difference in water vapor flux (kg m-2), at most can fill the pore space.
				} else {

					dM = std::max(  -((EMS[e].theta[WATER] - EMS[e].VG.theta_r * (1. + Constants::eps)) * Constants::density_water * EMS[e].L + EMS[e].theta[ICE] * Constants::density_ice * EMS[e].L)  ,
									std::min(  (EMS[e].theta[AIR] * Constants::density_ice * EMS[e].L), massPhaseChange  )
						 ); // mass change due to difference in water vapor flux (kg m-2), at most can fill the pore space.

				}
				

				// If there is no pore space, or, in fact, only so much pore space to accomodate the larger volume occupied by ice when all water freezes,
				// we inhibit vapour flux. This is necessary to maintain saturated conditions when present, and this is in turn necessary for the stability in the Richards equation solver.
				if(EMS[e].theta[AIR] < EMS[e].theta[WATER]*(Constants::density_water/Constants::density_ice - 1.) + Constants::eps) {
					dM = 0.;
				}

				//deltaM[e] += dM;
				deltaM[e] = dM;//JJJJJ
			}
			
			if(time==sn_dt) break;			 				
		}

		for (size_t i = 0; i < nE; i++) {
			//double topSaturatedVapor=Atmosphere::waterVaporDensity(NDS[i+1].T, Atmosphere::vaporSaturationPressure(NDS[i+1].T));
			//double botSaturatedVapor=Atmosphere::waterVaporDensity(NDS[i].T, Atmosphere::vaporSaturationPressure(NDS[i].T));
			//EMS[i].vapTrans_fluxDiff =-D_[i]*(topSaturatedVapor-botSaturatedVapor)/EMS[i].L;
			EMS[i].vapTrans_fluxDiff =-D_[i]*(NDS[i+1].rhov-NDS[i].rhov)/EMS[i].L;
		}

		double dHoar = 0.;
		for (e = 0; e < nE; e++)
		{	
			EMS[e].Qmm = 0.0;
			
			if (deltaM[e] < 0.) {
				// Mass loss: apply mass change first to water, then to ice, based on energy considerations
				// We can only do this partitioning here in this "simple" way, without checking if the mass is available, because we already limited dM above, based on available ICE + WATER.
				const double dTh_water = std::max( (EMS[e].VG.theta_r * (1. + Constants::eps) - EMS[e].theta[WATER])  ,  deltaM[e] / (Constants::density_water * EMS[e].L) );
				const double dTh_ice = ( deltaM[e] - (dTh_water * Constants::density_water * EMS[e].L) ) / (Constants::density_ice * EMS[e].L);
				EMS[e].theta[WATER] += dTh_water;
				EMS[e].theta[ICE] += dTh_ice;

				Sdata.mass[SurfaceFluxes::MS_SUBLIMATION] += dTh_water * Constants::density_water * EMS[e].L;
				Sdata.mass[SurfaceFluxes::MS_SUBLIMATION] += dTh_ice * Constants::density_ice * EMS[e].L;
				EMS[e].M += dTh_water * Constants::density_water * EMS[e].L+dTh_ice * Constants::density_ice * EMS[e].L;
				assert(EMS[e].M >= (-Constants::eps2)); //mass must be positive

				EMS[e].Qmm += (dTh_water*Constants::density_water*Constants::lh_vaporization +
								dTh_ice*Constants::density_ice*Constants::lh_sublimation)/sn_dt;//[w/m^3]

				// If present at surface, surface hoar is sublimated away
				if (e == nE-1 && deltaM[e]<0) {
					dHoar = std::max(-NDS[nN-1].hoar, deltaM[e]);
				}
			} else {		// Mass gain: add water in case temperature at or above melting point, ice otherwise
				if (EMS[e].Te >= EMS[e].meltfreeze_tk) {
					EMS[e].theta[WATER] += deltaM[e] / (Constants::density_water * EMS[e].L);
					EMS[e].Qmm += (deltaM[e]*Constants::lh_vaporization)/sn_dt/EMS[e].L;//  [w/m^3]
					Sdata.mass[SurfaceFluxes::MS_SUBLIMATION] += deltaM[e]; //
				} else {
					EMS[e].theta[ICE] += deltaM[e] / (Constants::density_ice * EMS[e].L);
					EMS[e].Qmm += (deltaM[e]*Constants::lh_sublimation)/sn_dt/EMS[e].L;// [w/m^3]
					Sdata.mass[SurfaceFluxes::MS_SUBLIMATION] += deltaM[e]; //
				}
				EMS[e].M += deltaM[e];
				assert(EMS[e].M >= (-Constants::eps2)); //mass must be positive
			}

			
			EMS[e].theta[AIR] = std::max(1. - EMS[e].theta[WATER] - EMS[e].theta[WATER_PREF] - EMS[e].theta[ICE] - EMS[e].theta[SOIL],0.);
			if(std::fabs(EMS[e].theta[AIR])<1.e-15)
			{
				EMS[e].theta[AIR]=0;
			}
			EMS[e].updDensity();
			assert(EMS[e].Rho > 0 || EMS[e].Rho==IOUtils::nodata); //density must be positive
			if (!(EMS[e].Rho > Constants::eps && EMS[e].theta[AIR] >= 0. && EMS[e].theta[WATER] <= 1.+Constants::eps && EMS[e].theta[ICE] <= 1.+Constants::eps)) {
					prn_msg(__FILE__, __LINE__, "err", Date(),
						"Volume contents: e=%d nE=%d rho=%lf ice=%lf wat=%lf wat_pref=%lf soil=%lf air=%le", e, nE, EMS[e].Rho, EMS[e].theta[ICE],
							EMS[e].theta[WATER], EMS[e].theta[WATER_PREF], EMS[e].theta[SOIL], EMS[e].theta[AIR]);
					throw IOException("Cannot evaluate mass balance in vapour transport LayerToLayer routine", AT);
			}
			
			/*
			// Numerical rounding errors were found to lead to theta[AIR] < 0, so force the other components between [0,1]:
			EMS[e].theta[ICE] = std::max(0., std::min(1. - EMS[e].theta[SOIL], EMS[e].theta[ICE]));
			EMS[e].theta[WATER] = std::max(0., std::min(1. - EMS[e].theta[SOIL], EMS[e].theta[WATER]));
			EMS[e].theta[WATER_PREF] = std::max(0., std::min(1., EMS[e].theta[WATER_PREF]));
			// Update theta[AIR] and density:
			EMS[e].theta[AIR] = (1. - EMS[e].theta[WATER] - EMS[e].theta[WATER_PREF] - EMS[e].theta[ICE] - EMS[e].theta[SOIL]);		
			EMS[e].updDensity();
			assert(EMS[e].Rho > 0 || EMS[e].Rho==IOUtils::nodata); //density must be positive
			if (!(EMS[e].Rho > Constants::eps && EMS[e].theta[AIR] >= 0.)) {
				if(EMS[e].theta[AIR] > -Constants::eps2) {
					EMS[e].theta[AIR] = 0.;
				} else {
					prn_msg(__FILE__, __LINE__, "err", Date(),
						"Volume contents: e=%d nE=%d rho=%lf ice=%lf wat=%lf wat_pref=%le air=%le  soil=%le", e, nE, EMS[e].Rho, EMS[e].theta[ICE], EMS[e].theta[WATER], EMS[e].theta[WATER_PREF], EMS[e].theta[AIR], EMS[e].salinity);
					throw IOException("Cannot evaluate mass balance in vapour transport LayerToLayer routine", AT);
				}
			}			
			*/
			
			//some useful output in case of vapor transport
			double sVaporDown = Atmosphere::waterVaporDensity(NDS[e].T, Atmosphere::vaporSaturationPressure(NDS[e].T));				
			double sVaporUp = Atmosphere::waterVaporDensity(NDS[e+1].T, Atmosphere::vaporSaturationPressure(NDS[e+1].T));
			EMS[e].vapTrans_underSaturationDegree = (0.5*(NDS[e].rhov-sVaporDown)+0.5*(NDS[e+1].rhov-sVaporUp))/(0.5*sVaporDown+0.5*sVaporUp);
			EMS[e].vapTrans_cumulativeDenChange += deltaM[e]/EMS[e].L;
			EMS[e].vapTrans_snowDenChangeRate = deltaM[e]/EMS[e].L/sn_dt;			
		}

		Sdata.hoar += dHoar;
		NDS[nN-1].hoar += dHoar;
		if (NDS[nN-1].hoar < 0.) {
			NDS[nN-1].hoar = 0.;
		}
		
	}
	else
	{
		compSurfaceSublimation(Mdata, ql, Xdata, Sdata);
		ql=0;
	}
	
}

/**
 * @brief Calculate the surface sublimation / deposition (i.e., only gas-solid). \n
 * The fraction of the latent heat flux ql that has not been used so far will be used for
 * sublimation/deposition. If positive (and above a certain cutoff level) then there
 * is a possibility that surface hoar crystal have grown. Of course, if negative
 * then we are also loosing mass from the surface.\n
 * This function additionally takes care of surface hoar formation and destruction.
 * Note that surface hoar is a nodal property, altough the corresponding mass is carried
 * by the underlying element.
 * @param *Mdata
 * @param ql Latent heat flux (W m-2)
 * @param *Xdata
 * @param *Sdata
 */
void VapourTransport::compSurfaceSublimation(const CurrentMeteo& Mdata, double& ql, SnowStation& Xdata, SurfaceFluxes& Sdata)
{
	double dL = 0., dM = 0.;     // Length and mass changes
	double M = 0.;               // Initial mass and volumetric content (water or ice)
	double dHoar = 0.;           // Actual change in hoar mass
	double cH_old;               // Temporary variable to hold height of snow

	const size_t nN = Xdata.getNumberOfNodes();
	size_t nE = nN-1;
	vector<NodeData>& NDS = Xdata.Ndata;
	vector<ElementData>& EMS = Xdata.Edata;
	const double Tss = NDS[nE].T; // Surface Temperature

	/*
	 * If ql > 0:
	 * Surface hoar is formed when surface temperature is below freezing.
	 * If no surface hoar can be formed, ql is kept and is used as boundary condition
	 * when calculating vapour flux.
	 * If there are elements and ql < 0:
	 * If ql is large enough to remove full surface elements, remove them.
	 * left over ql is used as boundary condition when calculating vapour flux.
	 *
	 * In both cases: add/subtract mass to MS_SUBLIMATION
	 */
	if (ql > Constants::eps2) { // Add Mass
		const double meltfreeze_tk = (Xdata.getNumberOfElements()>0)? Xdata.Edata[Xdata.getNumberOfElements()-1].meltfreeze_tk : Constants::meltfreeze_tk;
		if (Tss < meltfreeze_tk) { // Add Ice
			dM = ql*sn_dt/Constants::lh_sublimation;
			//if rh is very close to 1, vw too high or ta too high, surface hoar is destroyed and should not be formed
			if (!((Mdata.rh > hoar_thresh_rh) || (Mdata.vw > hoar_thresh_vw) || (Mdata.ta >= IOUtils::C_TO_K(hoar_thresh_ta)))) {
				// Under these conditions, form surface hoar
				ql = 0.;
				Sdata.mass[SurfaceFluxes::MS_SUBLIMATION] += dM;
				dHoar = dM;

				// In this case adjust properties of element, keeping snow density constant
				const double L_top = EMS[nE-1].L;
				const double theta_i0 = EMS[nE-1].theta[ICE];
				dL = dM/(EMS[nE-1].Rho); // length change
				if (nE == Xdata.SoilNode) {
					dL = 0.;
					dM = std::min(dM,EMS[nE-1].theta[AIR]*(Constants::density_ice*EMS[nE-1].L));
				}
				NDS[nE].z += dL + NDS[nE].u; NDS[nE].u = 0.0;
				EMS[nE-1].L0 = EMS[nE-1].L = L_top + dL;
				EMS[nE-1].E = EMS[nE-1].Eps = EMS[nE-1].dEps = EMS[nE-1].Eps_e = EMS[nE-1].Eps_v = EMS[nE-1].S = 0.0;
				EMS[nE-1].theta[ICE] *= L_top/EMS[nE-1].L;
				EMS[nE-1].theta[ICE] += dM/(Constants::density_ice*EMS[nE-1].L);
				EMS[nE-1].theta[ICE] = std::max(0., std::min(1., EMS[nE-1].theta[ICE]));
				EMS[nE-1].theta[WATER] *= L_top/EMS[nE-1].L;
				EMS[nE-1].theta[WATER] = std::max(0., std::min(1., EMS[nE-1].theta[WATER]));
				EMS[nE-1].theta[WATER_PREF] *= L_top/EMS[nE-1].L;
				EMS[nE-1].theta[WATER_PREF] = std::max(0., std::min(1., EMS[nE-1].theta[WATER_PREF]));

				for (size_t ii = 0; ii < Xdata.number_of_solutes; ii++) {
					EMS[nE-1].conc[ICE][ii] *= L_top*theta_i0/(EMS[nE-1].theta[ICE]*EMS[nE-1].L);
				}

				EMS[nE-1].M += dM;
				assert(EMS[nE-1].M >= (-Constants::eps2)); //mass must be positive

				// Update remaining volumetric contents and density
				EMS[nE-1].theta[AIR] = std::max(0., 1.0 - EMS[nE-1].theta[WATER] - EMS[nE-1].theta[WATER_PREF] - EMS[nE-1].theta[ICE] - EMS[nE-1].theta[SOIL]);
				EMS[nE-1].updDensity();
			}
		}
	} else if ((ql < (-Constants::eps2)) && (nE > 0)) {
		// If ql < 0, SUBLIMATE mass off
		std::vector<double> M_Solutes(Xdata.number_of_solutes, 0.); // Mass of solutes from disappearing phases
		size_t e = nE;
		while ((e > 0) && (ql < (-Constants::eps2))) {  // While energy is available
			e--;
			/*
			* Determine the amount of potential sublimation and collect some variables
			* that will be continuously used: L0 and M
			*/
			const double L0 = EMS[e].L;
			const double theta_i0 = EMS[e].theta[ICE];
			M = theta_i0*Constants::density_ice*L0;
			dM = ql*sn_dt/Constants::lh_sublimation;
			if (-dM > M) {
				// Only if mass change is sufficient to remove the full element
				dM = -M;
				// Add solutes to Storage
				for (size_t ii = 0; ii < Xdata.number_of_solutes; ii++) {
					M_Solutes[ii] += EMS[e].conc[ICE][ii]*theta_i0*L0;
				}
				EMS[e].theta[ICE] = 0.;
				dL = 0.;

				EMS[e].M += dM;
				Sdata.mass[SurfaceFluxes::MS_SUBLIMATION] += dM;
				ql -= dM*Constants::lh_sublimation/sn_dt;     // Update the energy used

				// If present at surface, surface hoar is sublimated away
				if (e == nE-1) {
					dHoar = std::max(-NDS[nN-1].hoar, dM);
				}

				// Update remaining volumetric contents and density
				EMS[e].theta[AIR] = std::max(0., 1.0 - EMS[e].theta[WATER] - EMS[e].theta[WATER_PREF] - EMS[e].theta[ICE] - EMS[e].theta[SOIL]);
				EMS[e].updDensity();
				// Merge the element if it is a snow layer. This will take care of possible left over liquid water (will be put one layer down)
				// Keep layer if it is a soil layer inside the snowpack (for example with snow farming)
				if(e>=Xdata.SoilNode) {
					if(EMS[e].theta[SOIL]<Constants::eps) {
						if (e>0) SnowStation::mergeElements(EMS[e-1], EMS[e], false, true);
						// Now reduce the number of elements by one.
						nE--;
					}
					//In case e==Xdata.SoilNode, we removed the last snow element and we should break out of the loop.
					if(e==Xdata.SoilNode) break;
				}
			} else {
				// Not enough energy anymore to remove complete element, so we should break out of the loop.
				break;
			}

			//check that thetas and densities are consistent
			assert(EMS[e].theta[SOIL] >= (-Constants::eps2) && EMS[e].theta[SOIL] <= (1.+Constants::eps2));
			assert(EMS[e].theta[ICE] >= (-Constants::eps2) && EMS[e].theta[ICE]<=(1.+Constants::eps2));
			assert(EMS[e].theta[WATER] >= (-Constants::eps2) && EMS[e].theta[WATER]<=(1.+Constants::eps2));
			assert(EMS[e].theta[WATER_PREF] >= (-Constants::eps2) && EMS[e].theta[WATER_PREF]<=(1.+Constants::eps2));
			assert(EMS[e].theta[AIR] >= (-Constants::eps2) && EMS[e].theta[AIR]<=(1.+Constants::eps2));
			assert(EMS[e].Rho >= (-Constants::eps2) || EMS[e].Rho==IOUtils::nodata); //we want positive density
		}

		// Now take care of left over solute mass.
		if (nE == Xdata.SoilNode) { // Add Solute Mass to Runoff TODO HACK CHECK
			for (size_t ii = 0; ii < Xdata.number_of_solutes; ii++) {
				Sdata.load[ii] += M_Solutes[ii]/S_TO_H(sn_dt);
			}
		} else { // Add Solute Mass to Element below
			if (EMS[e].theta[WATER] > 0.) {
				for(size_t ii = 0; ii < Xdata.number_of_solutes; ii++) {
					EMS[e].conc[WATER][ii] += M_Solutes[ii]/EMS[e].theta[WATER]/EMS[e].L;
				}
			} else if (EMS[e].theta[ICE] > 0.) {
				for (size_t ii = 0; ii < Xdata.number_of_solutes; ii++) {
					EMS[e].conc[ICE][ii] += M_Solutes[ii]/EMS[e].theta[ICE]/EMS[e].L;
				}
			} else {
				for (size_t ii = 0; ii < Xdata.number_of_solutes; ii++) {
					EMS[e].conc[SOIL][ii] += M_Solutes[ii]/EMS[e].theta[SOIL]/EMS[e].L;
				}
			}
		}
		Xdata.reduceNumberOfElements(nE);
	}

	// HACK: this code is under verification. The comment reads "surface hoar *is* destroyed, but the next line says surface hoar *may be* destroyed, depending on the sign of the latent heat flux.
	// If the code is correct, we can delete this part, if the comment is correct, we should modify the code to read: hoar = -NDS[nE].hoar;
	// Check for surface hoar destruction or formation (once upon a time ml_sn_SurfaceHoar)
	/*if ((Mdata.rh > hoar_thresh_rh) || (Mdata.vw > hoar_thresh_vw) || (Mdata.ta >= IOUtils::C_TO_K(hoar_thresh_ta))) {
		//if rh is very close to 1, vw too high or ta too high, surface hoar is destroyed
		hoar = std::min(hoar, 0.);
	}*/

	Sdata.hoar += dHoar;
	NDS[nN-1].hoar += dHoar;
	if (NDS[nN-1].hoar < 0.) {
		NDS[nN-1].hoar = 0.;
	}

	// Surface hoar cannot exist when the top element is wet
	if (nE > 0) {
		const double theta_r=((iwatertransportmodel_snow==RICHARDSEQUATION && nE-1>=Xdata.SoilNode) || (iwatertransportmodel_soil==RICHARDSEQUATION && nE-1<Xdata.SoilNode)) ? (PhaseChange::RE_theta_r) : (PhaseChange::theta_r);
		if (Xdata.Edata[nE-1].theta[WATER] > theta_r) {
			NDS[nE].hoar = 0.;
		}
	}

	// At the end also update the overall height
	cH_old = Xdata.cH;
	Xdata.cH = NDS[Xdata.getNumberOfNodes()-1].z + NDS[Xdata.getNumberOfNodes()-1].u;
	if (Xdata.mH!=Constants::undefined) Xdata.mH -= std::min(Xdata.mH - Xdata.Ground, (cH_old - Xdata.cH));	// TODO/HACK: why is this correction for Xdata.mH necessary?
}

bool VapourTransport::compDensityProfile(const CurrentMeteo& Mdata, SnowStation& Xdata, const bool& ThrowAtNoConvergence, double& ql, const double& surfaceVaporPressure, std::vector<double>& hm_, std::vector<double>& as_, const std::vector<double>& D_el, std::vector<double>& oldVaporDenNode)
{
	// Determine actual height of meteo values above Xdata.SoilNode:
	double actual_height_of_meteo_values;	// Height with reference Xdata.SoilNode
	if(!adjust_height_of_meteo_values) {
		// Case of fixed height above snow surface (e.g., weather model)
		actual_height_of_meteo_values = height_of_meteo_values + Xdata.cH - Xdata.Ground + ( (Xdata.findMarkedReferenceLayer() == Constants::undefined) ? (0.) : (Xdata.findMarkedReferenceLayer())  - Xdata.Ground);
	} else {
		// Case of fixed height above ground surface (e.g., weather station)
		actual_height_of_meteo_values = height_of_meteo_values;
	}
	
	bool bottomDirichletBCtype= false;
	
	const size_t nN = Xdata.getNumberOfNodes();
	size_t nE = nN-1;
	vector<NodeData>& NDS = Xdata.Ndata;
	vector<ElementData>& EMS = Xdata.Edata;
	const size_t nX = nN; // number of unknowns	

	//Built-in direct solvers:
	//SparseLU<SparseMatrix<double> > solver;
	//SparseQR<SparseMatrix<double> > solver;
	//SimplicialLLT<SparseMatrix<double> > solver;
	//SimplicialLDLT<SparseMatrix<double> > solver;

	//Built-in iterative solvers:	
	BiCGSTAB<SparseMatrix<double> > solver;
	//ConjugateGradient<SparseMatrix<double> > solver;
	//LeastSquaresConjugateGradient <SparseMatrix<double> > solver;
    //solver.setTolerance(1.e-45);
    //solver.setMaxIterations(1e6);
	
	SparseMatrix<double,RowMajor> A(nX,nX);
	//SparseMatrix<double> A(nX,nX);
	std::vector<Trip> tripletList(nX);
	VectorXd b(nX);
	VectorXd xx(nX);	
		
	// grid
    std::vector<double> z(nN,0.);
    for(size_t i=0; i<nN; i++)
    {   
		z[i]=NDS[i].z;
	}	
	
	// initial values
	std::vector<double> D_(nN, Constants::diffusion_coefficient_in_air);
	for(size_t i=0; i<nN; i++)
	{
		if(i==0)
		{
			D_[i]=D_el[i];
		}
		else if(i==nN-1)
		{
			D_[i]=0.5*D_el[i-1]+0.5*Constants::diffusion_coefficient_in_air;
		}
		else
		{
			D_[i]=D_el[i];
		}
	}
	
	// initial values
	// if the effective water vapor diffusivity in air is defined eps_=thata_a,
	// otherwise for the effective water vapor diffusivity in snow eps_=1.
	// for now Deff,s is available so eps_=1	
	std::vector<double> eps_(nN, 1.0);
	for(size_t i=0; i<nN; i++)
	{		
		if(i==0)
		{
			eps_[i]=std::max(EMS[i].theta[AIR],1.0);
		}
		else if(i==nN-1)
		{
			eps_[i]=0.5*std::max(EMS[i-1].theta[AIR],1.0)+0.5;
		}
		else
		{
			eps_[i]=std::max(EMS[i].theta[AIR],1.0);
		}					
	}
			
	double error_max = 0;
	do
	{
		error_max = 0;
						
        // the lower B.C.
		if(bottomDirichletBCtype){
			double elementSaturationVaporDensity=Atmosphere::waterVaporDensity(NDS[0].T, Atmosphere::vaporSaturationPressure(NDS[0].T));
			NDS[0].rhov=elementSaturationVaporDensity;
		}
		
		// diffusion equation
		double v_ij=0.;
		A.setZero();
		double dz_u,dz_d,dz_uu,dz_dd;
		double saturationDensity,hm,as,Rho_ice_water;
		double eps_n;
		for(size_t k=0; k<=nN-1; k++)
		{
			saturationDensity=Atmosphere::waterVaporDensity(NDS[k].T, Atmosphere::vaporSaturationPressure(NDS[k].T));			
																						
			if(k!=0 && k!=nN-1)
			{								
				dz_u=z[k+1]-z[k];
				dz_d=z[k]-z[k-1];
				eps_n=(0.5*eps_[k]+0.5*eps_[k-1]);
				
				b[k]= eps_n*NDS[k].rhov/timeStep+hm_[k]*as_[k]*(saturationDensity-(1-f)*NDS[k].rhov);
				b[k]+= -NDS[k].rhov*((1-f)*2.0*eps_[k]*D_[k]/dz_u/(dz_u+dz_d)+(1-f)*2.0*eps_[k-1]*D_[k-1]/dz_d/(dz_u+dz_d));
				b[k]+= -NDS[k+1].rhov*(1-f)*-2.0*eps_[k]*D_[k]/dz_u/(dz_u+dz_d);
				b[k]+= -NDS[k-1].rhov*(1-f)*-2.0*eps_[k-1]*D_[k-1]/dz_d/(dz_u+dz_d);
							    										
				v_ij=f*2.0*eps_[k]*D_[k]/dz_u/(dz_u+dz_d)+f*2.0*eps_[k-1]*D_[k-1]/dz_d/(dz_u+dz_d)+eps_n*1.0/timeStep+hm_[k]*as_[k]*f;
				tripletList.push_back(Trip(k,k,v_ij));//Set up the matrix diagonal
							  
				v_ij=f*-2.0*eps_[k]*D_[k]/dz_u/(dz_u+dz_d);	  
				tripletList.push_back(Trip(k,k+1,v_ij));//Set up the matrix upper diagonals, k+1
								  								  
				v_ij=f*-2.0*eps_[k-1]*D_[k-1]/dz_d/(dz_u+dz_d);	  
				tripletList.push_back(Trip(k,k-1,v_ij));//Set up the matrix lower diagonals, k-1

				/*
				if(k>=Xdata.SoilNode)
				{
					std::cout << "k/temporal/diffusion/hmas/b " <<  k << ' ' <<	eps_n/timeStep << ' ' 
															<< eps_[k]*f*2.0*D_[k]/dz_u/(dz_u+dz_d)+eps_[k-1]*f*2.0*D_[k-1]/dz_d/(dz_u+dz_d) << ' ' 
															<< hm_[k]*as_[k]  << ' ' << b[k] << "\n";
				}
				*/
				
			}
			if(k==nN-1)
			{	
	
				//normal top B.C. assuming satuarion condition for the uppermost node of snowpack
				//if(Xdata.SoilNode<nN-1)
				//{
					b[k]=saturationDensity;
					v_ij=1.0;
					tripletList.push_back(Trip(k,k,v_ij));//Set up the matrix diagonal
				//}
				/*
				//normal top B.C. assuming zero flux of the uppermost element of snowpack
				b[k]=0;

				v_ij=1.;
				tripletList.push_back(Trip(k,k,v_ij));//Set up the matrix diagonal

				v_ij=-1.;	  
				tripletList.push_back(Trip(k,k-1,v_ij));//Set up the matrix lower diagonals, k-1
				*/
				
				
				//normal top B.C. assuming diffusive mass flux of the uppermost element of snowpack equals to mass flux from(to) atmoshpere to(from) snowpack
				/*
				if(Xdata.SoilNode==nN-1)
				{
					double beta = SnLaws::compLatentHeat(Snowpack::EVAP_RELATIVE_HUMIDITY, Mdata, Xdata, actual_height_of_meteo_values);
					double p_vapor = Mdata.rh * Atmosphere::vaporSaturationPressure(Mdata.ta);
					double rhov_atm = Atmosphere::waterVaporDensity(Mdata.ta, p_vapor);		
					dz_d=z[k]-z[k-1];
					double Xatm = 28.9647e-3/8.31451/Mdata.ta;				
					double Xsurf = 28.9647e-3/8.31451/NDS[k].T;

					b[k]=beta/Xatm/Constants::lh_sublimation*rhov_atm;

					v_ij=beta/Xsurf/Constants::lh_sublimation+D_[k]/dz_d;
					tripletList.push_back(Trip(k,k,v_ij));//Set up the matrix diagonal

					v_ij=-D_[k]/dz_d;	  
					tripletList.push_back(Trip(k,k-1,v_ij));//Set up the matrix lower diagonals, k-1
				}
				*/
				
				/*				
				//assuming the net flux divergence of the uppermost element of snowpack is zero

				eps_n=eps_[k]; // volumetric air content for node k=nN-1

				b[k]=eps_n*NDS[k].rhov/timeStep+hm_[k]*as_[k]*saturationDensity;

				v_ij=eps_n*1.0/timeStep+hm_[k]*as_[k];
				tripletList.push_back(Trip(k,k,v_ij));//Set up the matrix diagonal
				*/
								
				/*				
				//BB
				double beta = SnLaws::compLatentHeat(Snowpack::EVAP_RELATIVE_HUMIDITY, Mdata, Xdata, actual_height_of_meteo_values);
				double p_vapor = Mdata.rh * Atmosphere::vaporSaturationPressure(Mdata.ta);
				double rhov_atm = Atmosphere::waterVaporDensity(Mdata.ta, p_vapor);		
				dz_d=z[k]-z[k-1];
				dz_dd=z[k-1]-z[k-2];
				double Xatm = 28.9647e-3/8.31451/Mdata.ta;				
				double Xsurf = 28.9647e-3/8.31451/NDS[k].T;

				eps_n=eps_[k]; // volumetric air content for node k=nN-1
				double D_nd = 0.5*D_[k-1]+0.5*D_[k-2]; // diffusivity for node k=nN-2
				double eps_nd = 0.5*eps_[k-1]+0.5*eps_[k-2]; // volumetric air content for node k=nN-2

				b[k]=eps_n*NDS[k].rhov/timeStep+hm_[k]*as_[k]*saturationDensity+eps_n*beta/Xatm/Constants::lh_sublimation/dz_d*rhov_atm;

				v_ij=eps_n*1.0/timeStep+hm_[k]*as_[k]+eps_n*beta/Xsurf/Constants::lh_sublimation/dz_d;
				tripletList.push_back(Trip(k,k,v_ij));//Set up the matrix diagonal

				v_ij=eps_nd*D_nd/dz_d/dz_dd;	  
				tripletList.push_back(Trip(k,k-1,v_ij));//Set up the matrix lower diagonals, k-1
				
				v_ij=eps_nd*-D_nd/dz_d/dz_dd;	  
				tripletList.push_back(Trip(k,k-2,v_ij));//Set up the matrix lower diagonals, k-2
				*/
				
				/*
				//BC
				double beta = SnLaws::compLatentHeat(Snowpack::EVAP_RELATIVE_HUMIDITY, Mdata, Xdata, actual_height_of_meteo_values);
				double p_vapor = Mdata.rh * Atmosphere::vaporSaturationPressure(Mdata.ta);
				double rhov_atm = Atmosphere::waterVaporDensity(Mdata.ta, p_vapor);		
				dz_d=z[k]-z[k-1];
				dz_dd=z[k-1]-z[k-2];
				double Xatm = 28.9647e-3/8.31451/Mdata.ta;				
				double Xsurf = 28.9647e-3/8.31451/NDS[k].T;				

				b[k]=eps_[k]*NDS[k].rhov/timeStep+hm_[k]*as_[k]*saturationDensity+eps_[k]*beta/Xatm/Constants::lh_sublimation/dz_d*rhov_atm;

				v_ij=eps_[k-1]*D_[k-1]/dz_d/(dz_d+dz_dd)+eps_[k]*1.0/timeStep+eps_[k]*beta/Xsurf/Constants::lh_sublimation/dz_d+hm_[k]*as_[k];
				tripletList.push_back(Trip(k,k,v_ij));//Set up the matrix diagonal

				v_ij=eps_[k-1]*-D_[k-1]/dz_d/(dz_d+dz_dd);	  
				tripletList.push_back(Trip(k,k-2,v_ij));//Set up the matrix lower diagonals, k-2
				*/
				
				/*
				//BF
				double beta = SnLaws::compLatentHeat(Snowpack::EVAP_RELATIVE_HUMIDITY, Mdata, Xdata, actual_height_of_meteo_values);
				double p_vapor = Mdata.rh * Atmosphere::vaporSaturationPressure(Mdata.ta);
				double rhov_atm = Atmosphere::waterVaporDensity(Mdata.ta, p_vapor);		
				dz_d=z[k]-z[k-1];
				double Xatm = 28.9647e-3/8.31451/Mdata.ta;				
				double Xsurf = 28.9647e-3/8.31451/NDS[k].T;				

				b[k]= NDS[k].rhov/timeStep+hm_[k]*as_[k]*saturationDensity+beta/Xatm/Constants::lh_sublimation/dz_d*rhov_atm;

				v_ij=D_[k-1]/dz_d/dz_d+1.0/timeStep+beta/Xsurf/Constants::lh_sublimation/dz_d+hm_[k]*as_[k];
				tripletList.push_back(Trip(k,k,v_ij));//Set up the matrix diagonal

				v_ij=-D_[k-1]/dz_d/dz_d;	  
				tripletList.push_back(Trip(k,k-1,v_ij));//Set up the matrix lower diagonals, k-1
				*/		
			}
			if(k==0)
			{
				if(bottomDirichletBCtype)
				{
					b[k]= saturationDensity;//NDS[k].rhov;		
					v_ij= 1.0;
					tripletList.push_back(Trip(k,k,v_ij));//Set up the matrix diagonal					
				}
				else
				{				
					dz_u=z[1]-z[0];
					dz_uu=z[2]-z[1];

					b[k]= 0.0;		
					v_ij= -1.0;
					tripletList.push_back(Trip(k,k,v_ij));//Set up the matrix diagonal
					
					v_ij= 1.0;
					tripletList.push_back(Trip(k,k+1,v_ij));//Set up the matrix upper diagonals, k+1					
				}
			}
		}

		A.setFromTriplets(tripletList.begin(), tripletList.end());
		tripletList.clear();
		A.makeCompressed();
		solver.compute(A);
		if(solver.info()!=Success) 
		{
			prn_msg(__FILE__, __LINE__, "msg", Mdata.date, "solver.compute(A).......you cannot use solver Eigen....vapor transport is ignored fro this time step");
			return false;
			//throw InvalidArgumentException("solver.compute(A).......you cannot use solver Eigen....", AT);
		}			
		//solve the equation
		xx = solver.solve(b);
		//std::cout << "estimated error: " << solver.error()      << std::endl;
		if(solver.info()!=Success) 
		{
			//throw InvalidArgumentException("solver.solve(b).........you cannot use solver Eigen....", AT);
		}				

		std::vector<double> oldrhov(nN, 0);	
		for(size_t k=0; k<=nN-1; k++)
		{	
			oldVaporDenNode[k]=NDS[k].rhov;				
			NDS[k].rhov=xx(k);
			double error = std::abs(NDS[k].rhov-oldVaporDenNode[k]);
			error_max=std::max(error_max,error);
			saturationDensity=Atmosphere::waterVaporDensity(NDS[k].T, Atmosphere::vaporSaturationPressure(NDS[k].T));
			if(NDS[k].rhov<0)
			{
				std::cout << "...........3k/nN-1/oldrhov/NDS[k].rhov/saturationDensity.................. " << k << ' ' 
				          << nN-1 << ' ' << oldVaporDenNode[k] << ' ' << NDS[k].rhov << ' ' << saturationDensity << "\n";						
				//NDS[k].rhov=1.e-9;          
			}			
		}
		
		break;
	}while(error_max>1.e-6);
		
	for (size_t e = 0; e < nE; e++)
	{
		EMS[e].rhov = (NDS[e].rhov + NDS[e+1].rhov) / 2.0;
	}
			 
	return true;
}

double VapourTransport::dRhov_dT(const double Tem)
{
	double c2,c3;
	const double water_molecular_mass = 18.0153e-3; // (kg)
	const double p_water_triple_pt = 611.73; // (Pa)
	const double gaz_constant = 8.31451; // (J mol-1 K-1)
	const double t_water_triple_pt = 273.16; // (K)
	const double l_water_sublimation = 2.838e6; // (J kg-1)
		
	// For the first derivative of vapor density with respect to temperature, we have
	if ( Tem < t_water_triple_pt ) { // for a flat ice surface
		c2 = 21.88;
		c3 = 7.66;
	} else { // for a flat water surface
		c2 = 17.27;
		c3 = 35.86;
	}
	c2 = 21.88;
	c3 = 7.66;
	
	double dRhov_dT = (water_molecular_mass*p_water_triple_pt/gaz_constant/Tem)*
				                 exp(c2*(Tem-t_water_triple_pt)/(Tem-c3))*
			                     ( -1./Tem + (c2*(Tem-c3)-c2*(Tem-t_water_triple_pt))/(Tem-c3)/(Tem-c3) );			                     
	return dRhov_dT;
}
/*
 * End of VapourTransport.cc
 */
