

/*! TODO:
  - move these functions to separate header file with KPP code
  - copy header file to libs/build/include ...
  - doxygen comments are missing.
  - change kppchecm* to kpp_chem*
*/

void kpp_chem_bound_cond(
  ctl_t * ctl,
  atm_t * atm,
  met_t * met0,
  met_t * met1,
  int ip);

void kpp_chem_initialize(
  ctl_t * ctl,
  clim_t * clim,
  met_t * met0,
  met_t * met1,
  atm_t * atm,
  int ip);

void kpp_chem_output2atm(
  atm_t * atm,
  ctl_t * ctl,
  int ip);

void kpp_chemgrid_mass2concen(
  atm_t * atm,
  ctl_t * ctl,
  double *mass,
  double *area,
  int *ixs,
  int *iys,
  int *izs,
  double dz,
  int ip,
  int qnt_index);

void kpp_chem_init_cqnt(
  ctl_t * ctl,
  atm_t * atm,
	clim_t * clim,
  met_t * met0,
  met_t * met1,
  int ip);

double param_mixing_calc(
  ctl_t * ctl,
  clim_t * clim,
  atm_t * atm,
  int ip);

void interparc_mixing(
  ctl_t * ctl,
  atm_t * atm,
  clim_t * clim,
  int *ixs,
  int *iys,
  int *izs);

void interparc_mixing_help(
  ctl_t * ctl,
  atm_t * atm,
  clim_t * clim,
  int *ixs,
  int *iys,
  int *izs,
  int qnt_idx);