// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// 
// Auxiliary Routines File
// 
// Generated by KPP-3.0.0 symbolic chemistry Kinetics PreProcessor
//       (https:/github.com/KineticPreProcessor/KPP
// KPP is distributed under GPL, the general public licence
//       (http://www.gnu.org/copyleft/gpl.html)
// (C) 1995-1997, V. Damian & A. Sandu, CGRER, Univ. Iowa
// (C) 1997-2022, A. Sandu, Michigan Tech, Virginia Tech
//     With important contributions from:
//        M. Damian,   Villanova University, Philadelphia, PA, USA
//        R. Sander,   Max-Planck Institute for Chemistry, Mainz, Germany
//        M. Long,     Renaissance Fiber, LLC, North Carolina, USA
//        H. Lin,      Harvard University, Cambridge, MA, USA
//        R. Yantosca, Harvard University, Cambridge, MA, USA
// 
// File                 : chem_Util.c
// Equation file        : chem.kpp
// Output root filename : chem
// 
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "chem_Parameters.h"
#include "chem_Global.h"
#include "chem_Sparse.h"



// User INLINED Utility Functions


#include "libtrac.h"

void kppchem_output2atm(
	atm_t * atm,
	ctl_t * ctl, 
	int ip){
	
		/*Output to air parcel.. */
    if (ctl->qnt_Ch2o2 > 0)
			atm->q[ctl->qnt_Ch2o2][ip] = VAR[ind_H2O2];
    if (atm->q[ctl->qnt_Cx][ip] != 0) {
			atm->q[ctl->qnt_m][ip] *= VAR[ind_SO2] / atm->q[ctl->qnt_Cx][ip];
			atm->q[ctl->qnt_Cx][ip] = VAR[ind_SO2];
		}
}

void kppchem_init_qntvar(
	ctl_t * ctl,
	atm_t * atm,
	clim_t * clim,
	int ip
){

	if (ctl->qnt_Ch2o2 > 0)	
	if (atm->q[ctl->qnt_Ch2o2][ip] == 0)
		atm->q[ctl->qnt_Ch2o2][ip] = clim_h2o2(clim, atm->time[ip], atm->lat[ip], atm->p[ip]);	//Unit: molec/cm3
if (ctl->qnt_Cho2 > 0)
	if (atm->q[ctl->qnt_Cho2][ip] == 0)
		atm->q[ctl->qnt_Cho2][ip] =
			clim_ho2(clim, atm->time[ip], atm->lat[ip], atm->p[ip]);
}



// End INLINED Utility Functions

// Utility Functions from KPP_HOME/util/util
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// 
// UTIL - Utility functions
//   Arguments :
// 
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

/*
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
*/
double min( double x, double y )
{
  return ( x <= y ) ? x : y;
}

double max( double x, double y )
{
  return ( x >= y ) ? x : y;
}

static FILE *fpDat = 0;

int InitSaveData()
{
  fpDat = fopen("chem.dat", "w");
  if( fpDat == 0 ) {
    printf("\n Can't create file : chem.dat");
    exit(1);
  }
  return 0;
}

int SaveData()
{
int i;

  fprintf( fpDat, "%6.1f ", TIME/3600.0 );
  for( i = 0; i < NLOOKAT; i++ )
    fprintf( fpDat, "%24.16e ", C[ LOOKAT[i] ]/CFACTOR );
  fprintf( fpDat, "\n");
  return 0;
}

int CloseSaveData()
{
  fclose( fpDat );
  return 0;
}

int GenerateMatlab( char * prefix )
{
int i;
FILE *fpMatlab;
  
  fpMatlab = fopen("chem.m", "w");
  if( fpMatlab == 0 ) {
    printf("\n Can't create file : chem.m");
    exit(1);
  }

  fprintf(fpMatlab, "load chem.dat;\n");
  fprintf(fpMatlab, "%sc = chem;\n", prefix);
  fprintf(fpMatlab, "clear chem;\n");
  fprintf(fpMatlab, "%st=%sc(:,1);\n", prefix, prefix);
  fprintf(fpMatlab, "%sc(:,1)=[];\n", prefix);
  
  for( i = 0; i < NLOOKAT; i++ )
    fprintf( fpMatlab, "%s%s = %sc(:,%d);\n", 
            prefix, SPC_NAMES[LOOKAT[i]], 
            prefix, i+1 );
  
  fclose( fpMatlab );
  return 0;
}

// End Utility Functions from KPP_HOME/util/util
// End of UTIL function
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// 
// Shuffle_user2kpp - function to copy concentrations from USER to KPP
//   Arguments :
//      V_USER    - Concentration of variable species in USER's order
//      V         - Concentrations of variable species (local)
// 
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

void Shuffle_user2kpp( 
  double V_USER[NVAR],                   /* Concentration of variable species in USER's order */
  double V[NVAR]                         /* Concentrations of variable species (local) */
)
{
  V[2] = V_USER[0];
  V[1] = V_USER[1];
  V[0] = V_USER[2];
}

// End of Shuffle_user2kpp function
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// 
// Shuffle_kpp2user - function to restore concentrations from KPP to USER
//   Arguments :
//      V         - Concentrations of variable species (local)
//      V_USER    - Concentration of variable species in USER's order
// 
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

void Shuffle_kpp2user( 
  double V[NVAR],                        /* Concentrations of variable species (local) */
  double V_USER[NVAR]                    /* Concentration of variable species in USER's order */
)
{
  V_USER[0] = V[2];
  V_USER[1] = V[1];
  V_USER[2] = V[0];
}

// End of Shuffle_kpp2user function
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// 
// GetMass - compute total mass of selected atoms
//   Arguments :
//      CL        - Concentration of all species (local)
//      Mass      - value of mass balance
// 
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

void GetMass( 
  double CL[NSPEC],                      /* Concentration of all species (local) */
  double Mass[1]                         /* value of mass balance */
)
{
}

// End of GetMass function
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

