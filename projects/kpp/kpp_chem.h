
/*! Set variable species concentration. */
#define SET_VAR(ind_spec, qnt_Cspec, M) 			\
  if (ctl->qnt_Cspec >= 0)					\
    var[ind_spec] = atm->q[ctl->qnt_Cspec][ip] * M;

#define SET_RCONST(ind_react, rconst_value) 			\
  rconst[ind_react] = rconst_value;

/*! Get variable species concentration. */
#define GET_VAR(ind_spec, qnt_index)			\
  if (qnt_index >= 0)					\
    atm->q[qnt_index][ip] = var[ind_spec] / M;

/*! Initialize concentration quantity. */
#define INIT_CQNT(qnt_index, clim_zm_t)					\
  if (qnt_index >= 0)							\
    atm->q[qnt_index][ip] =						\
      clim_zm(&clim_var_t, atm->time[ip], atm->lat[ip], atm->p[ip]);

/*! Initialize KPP chemistry. */
#ifdef _OPENACC
#pragma acc routine (kpp_chem_initialize)
#endif
void kpp_chem_initialize(
  ctl_t * ctl,
  clim_t * clim,
  met_t * met0,
  met_t * met1,
  atm_t * atm,
  double * var,
  double * fix,
  double * rconst,
  int ip);

/*! Save KPP results. */
#ifdef _OPENACC
#pragma acc routine (kpp_chem_output2atm)
#endif
void kpp_chem_output2atm(
  atm_t * atm,
  ctl_t * ctl,
  met_t * met0,
	met_t * met1,
  double * var,
  int ip);

/*! KPP integration function. */
#ifdef _OPENACC
#pragma acc routine (Rosenbrock)
#endif
int Rosenbrock(double Y[], double fix[], double rconst[], double Tstart, double Tend,
     double AbsTol[], double RelTol[],
    //  void (*ode_Fun)(double, double [], double []), 
     void (*ode_Fun)(double, double [], double [], double [], double []), 
    //  void (*ode_Jac)(double, double [], double []),
    void (*ode_Jac)(double, double [], double [], double [], double []), 
     double RPAR[], int IPAR[]);
#ifdef _OPENACC
#pragma acc routine (FunTemplate)
#endif
void FunTemplate(double, double [], double [], double [], double []); 
#ifdef _OPENACC
#pragma acc routine (JacTemplate)
#endif
void JacTemplate(double, double [], double [], double [], double []);