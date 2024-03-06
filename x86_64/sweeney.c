/* Created by Language version: 7.7.0 */
/* NOT VECTORIZED */
#define NRN_VECTORIZED 0
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "mech_api.h"
#undef PI
#define nil 0
#include "md1redef.h"
#include "section.h"
#include "nrniv_mf.h"
#include "md2redef.h"
 
#if METHOD3
extern int _method3;
#endif

#if !NRNGPU
#undef exp
#define exp hoc_Exp
extern double hoc_Exp(double);
#endif
 
#define nrn_init _nrn_init__sweeney
#define _nrn_initial _nrn_initial__sweeney
#define nrn_cur _nrn_cur__sweeney
#define _nrn_current _nrn_current__sweeney
#define nrn_jacob _nrn_jacob__sweeney
#define nrn_state _nrn_state__sweeney
#define _net_receive _net_receive__sweeney 
#define rates rates__sweeney 
#define states states__sweeney 
 
#define _threadargscomma_ /**/
#define _threadargsprotocomma_ /**/
#define _threadargs_ /**/
#define _threadargsproto_ /**/
 	/*SUPPRESS 761*/
	/*SUPPRESS 762*/
	/*SUPPRESS 763*/
	/*SUPPRESS 765*/
	 extern double *getarg();
 static double *_p; static Datum *_ppvar;
 
#define t nrn_threads->_t
#define dt nrn_threads->_dt
#define gnabar _p[0]
#define gnabar_columnindex 0
#define gl _p[1]
#define gl_columnindex 1
#define el _p[2]
#define el_columnindex 2
#define il _p[3]
#define il_columnindex 3
#define minf _p[4]
#define minf_columnindex 4
#define hinf _p[5]
#define hinf_columnindex 5
#define mtau _p[6]
#define mtau_columnindex 6
#define htau _p[7]
#define htau_columnindex 7
#define mexp _p[8]
#define mexp_columnindex 8
#define hexp _p[9]
#define hexp_columnindex 9
#define m _p[10]
#define m_columnindex 10
#define h _p[11]
#define h_columnindex 11
#define ena _p[12]
#define ena_columnindex 12
#define Dm _p[13]
#define Dm_columnindex 13
#define Dh _p[14]
#define Dh_columnindex 14
#define ina _p[15]
#define ina_columnindex 15
#define _g _p[16]
#define _g_columnindex 16
#define _ion_ena	*_ppvar[0]._pval
#define _ion_ina	*_ppvar[1]._pval
#define _ion_dinadv	*_ppvar[2]._pval
 
#if MAC
#if !defined(v)
#define v _mlhv
#endif
#if !defined(h)
#define h _mlhh
#endif
#endif
 
#if defined(__cplusplus)
extern "C" {
#endif
 static int hoc_nrnpointerindex =  -1;
 /* external NEURON variables */
 extern double celsius;
 /* declaration of user functions */
 static void _hoc_Exp(void);
 static void _hoc_rates(void);
 static void _hoc_states(void);
 static int _mechtype;
extern void _nrn_cacheloop_reg(int, int);
extern void hoc_register_prop_size(int, int, int);
extern void hoc_register_limits(int, HocParmLimits*);
extern void hoc_register_units(int, HocParmUnits*);
extern void nrn_promote(Prop*, int, int);
extern Memb_func* memb_func;
 
#define NMODL_TEXT 1
#if NMODL_TEXT
static const char* nmodl_file_text;
static const char* nmodl_filename;
extern void hoc_reg_nmodl_text(int, const char*);
extern void hoc_reg_nmodl_filename(int, const char*);
#endif

 extern void _nrn_setdata_reg(int, void(*)(Prop*));
 static void _setdata(Prop* _prop) {
 _p = _prop->param; _ppvar = _prop->dparam;
 }
 static void _hoc_setdata() {
 Prop *_prop, *hoc_getdata_range(int);
 _prop = hoc_getdata_range(_mechtype);
   _setdata(_prop);
 hoc_retpushx(1.);
}
 /* connect user functions to hoc names */
 static VoidFunc hoc_intfunc[] = {
 "setdata_sweeney", _hoc_setdata,
 "Exp_sweeney", _hoc_Exp,
 "rates_sweeney", _hoc_rates,
 "states_sweeney", _hoc_states,
 0, 0
};
#define Exp Exp_sweeney
 extern double Exp( double );
 /* declare global and static user variables */
#define ahB ahB_sweeney
 double ahB = 15.6;
#define ahA ahA_sweeney
 double ahA = 56;
#define amD amD_sweeney
 double amD = 5.3;
#define amC amC_sweeney
 double amC = 0.363;
#define amB amB_sweeney
 double amB = 126;
#define amA amA_sweeney
 double amA = 49;
#define bhC bhC_sweeney
 double bhC = 5;
#define bhB bhB_sweeney
 double bhB = 74.5;
#define bhA bhA_sweeney
 double bhA = 10;
#define bmB bmB_sweeney
 double bmB = 4.17;
#define bmA bmA_sweeney
 double bmA = 56.2;
#define vtraub vtraub_sweeney
 double vtraub = 0;
 /* some parameters have upper and lower limits */
 static HocParmLimits _hoc_parm_limits[] = {
 0,0,0
};
 static HocParmUnits _hoc_parm_units[] = {
 "vtraub_sweeney", "mV",
 "gnabar_sweeney", "mho/cm2",
 "gl_sweeney", "mho/cm2",
 "el_sweeney", "mV",
 "il_sweeney", "mA/cm2",
 "mtau_sweeney", "ms",
 "htau_sweeney", "ms",
 0,0
};
 static double delta_t = 1;
 static double h0 = 0;
 static double m0 = 0;
 static double v = 0;
 /* connect global user variables to hoc */
 static DoubScal hoc_scdoub[] = {
 "vtraub_sweeney", &vtraub_sweeney,
 "amA_sweeney", &amA_sweeney,
 "amB_sweeney", &amB_sweeney,
 "amC_sweeney", &amC_sweeney,
 "amD_sweeney", &amD_sweeney,
 "bmA_sweeney", &bmA_sweeney,
 "bmB_sweeney", &bmB_sweeney,
 "ahA_sweeney", &ahA_sweeney,
 "ahB_sweeney", &ahB_sweeney,
 "bhA_sweeney", &bhA_sweeney,
 "bhB_sweeney", &bhB_sweeney,
 "bhC_sweeney", &bhC_sweeney,
 0,0
};
 static DoubVec hoc_vdoub[] = {
 0,0,0
};
 static double _sav_indep;
 static void nrn_alloc(Prop*);
static void  nrn_init(NrnThread*, _Memb_list*, int);
static void nrn_state(NrnThread*, _Memb_list*, int);
 static void nrn_cur(NrnThread*, _Memb_list*, int);
static void  nrn_jacob(NrnThread*, _Memb_list*, int);
 
static int _ode_count(int);
 /* connect range variables in _p that hoc is supposed to know about */
 static const char *_mechanism[] = {
 "7.7.0",
"sweeney",
 "gnabar_sweeney",
 "gl_sweeney",
 "el_sweeney",
 0,
 "il_sweeney",
 "minf_sweeney",
 "hinf_sweeney",
 "mtau_sweeney",
 "htau_sweeney",
 "mexp_sweeney",
 "hexp_sweeney",
 0,
 "m_sweeney",
 "h_sweeney",
 0,
 0};
 static Symbol* _na_sym;
 
extern Prop* need_memb(Symbol*);

static void nrn_alloc(Prop* _prop) {
	Prop *prop_ion;
	double *_p; Datum *_ppvar;
 	_p = nrn_prop_data_alloc(_mechtype, 17, _prop);
 	/*initialize range parameters*/
 	gnabar = 1.445;
 	gl = 0.128;
 	el = -80.01;
 	_prop->param = _p;
 	_prop->param_size = 17;
 	_ppvar = nrn_prop_datum_alloc(_mechtype, 3, _prop);
 	_prop->dparam = _ppvar;
 	/*connect ionic variables to this model*/
 prop_ion = need_memb(_na_sym);
 nrn_promote(prop_ion, 0, 1);
 	_ppvar[0]._pval = &prop_ion->param[0]; /* ena */
 	_ppvar[1]._pval = &prop_ion->param[3]; /* ina */
 	_ppvar[2]._pval = &prop_ion->param[4]; /* _ion_dinadv */
 
}
 static void _initlists();
 static void _update_ion_pointer(Datum*);
 extern Symbol* hoc_lookup(const char*);
extern void _nrn_thread_reg(int, int, void(*)(Datum*));
extern void _nrn_thread_table_reg(int, void(*)(double*, Datum*, Datum*, NrnThread*, int));
extern void hoc_register_tolerance(int, HocStateTolerance*, Symbol***);
extern void _cvode_abstol( Symbol**, double*, int);

 void _sweeney_reg() {
	int _vectorized = 0;
  _initlists();
 	ion_reg("na", -10000.);
 	_na_sym = hoc_lookup("na_ion");
 	register_mech(_mechanism, nrn_alloc,nrn_cur, nrn_jacob, nrn_state, nrn_init, hoc_nrnpointerindex, 0);
 _mechtype = nrn_get_mechtype(_mechanism[1]);
     _nrn_setdata_reg(_mechtype, _setdata);
     _nrn_thread_reg(_mechtype, 2, _update_ion_pointer);
 #if NMODL_TEXT
  hoc_reg_nmodl_text(_mechtype, nmodl_file_text);
  hoc_reg_nmodl_filename(_mechtype, nmodl_filename);
#endif
  hoc_register_prop_size(_mechtype, 17, 3);
  hoc_register_dparam_semantics(_mechtype, 0, "na_ion");
  hoc_register_dparam_semantics(_mechtype, 1, "na_ion");
  hoc_register_dparam_semantics(_mechtype, 2, "na_ion");
 	hoc_register_cvode(_mechtype, _ode_count, 0, 0, 0);
 	hoc_register_var(hoc_scdoub, hoc_vdoub, hoc_intfunc);
 	ivoc_help("help ?1 sweeney /home/merolanr/repos/BME515_HW2/sweeney.mod\n");
 hoc_register_limits(_mechtype, _hoc_parm_limits);
 hoc_register_units(_mechtype, _hoc_parm_units);
 }
static int _reset;
static char *modelname = "";

static int error;
static int _ninits = 0;
static int _match_recurse=1;
static void _modl_cleanup(){ _match_recurse=1;}
static int rates(double);
static int states();
 
static int  states (  ) {
   rates ( _threadargscomma_ v ) ;
   m = m + mexp * ( minf - m ) ;
   h = h + hexp * ( hinf - h ) ;
   
/*VERBATIM*/
	return 0;
  return 0; }
 
static void _hoc_states(void) {
  double _r;
   _r = 1.;
 states (  );
 hoc_retpushx(_r);
}
 
static int  rates (  double _lv ) {
   double _lv2 , _lalpha , _lbeta , _lsum ;
 _lv2 = _lv - vtraub ;
   _lalpha = ( amB + amC * _lv2 ) / ( 1.0 + Exp ( _threadargscomma_ - ( amA + _lv2 ) / amD ) ) ;
   _lbeta = ( 1.0 / Exp ( _threadargscomma_ ( _lv2 + bmA ) / bmB ) ) * ( amB + amC * _lv2 ) / ( 1.0 + Exp ( _threadargscomma_ - ( amA + _lv2 ) / amD ) ) ;
   _lsum = _lalpha + _lbeta ;
   mtau = 1.0 / _lsum ;
   minf = _lalpha / _lsum ;
   _lalpha = ( 1.0 / Exp ( _threadargscomma_ ( _lv2 + bhB ) / bhC ) ) * ahB / ( 1.0 + Exp ( _threadargscomma_ - ( _lv2 + ahA ) / bhA ) ) ;
   _lbeta = ahB / ( 1.0 + Exp ( _threadargscomma_ - ( _lv2 + ahA ) / bhA ) ) ;
   _lsum = _lalpha + _lbeta ;
   htau = 1.0 / _lsum ;
   hinf = _lalpha / _lsum ;
   mexp = 1.0 - Exp ( _threadargscomma_ - dt / mtau ) ;
   hexp = 1.0 - Exp ( _threadargscomma_ - dt / htau ) ;
    return 0; }
 
static void _hoc_rates(void) {
  double _r;
   _r = 1.;
 rates (  *getarg(1) );
 hoc_retpushx(_r);
}
 
double Exp (  double _lx ) {
   double _lExp;
 if ( _lx < - 100.0 ) {
     _lExp = 0.0 ;
     }
   else {
     _lExp = exp ( _lx ) ;
     }
   
return _lExp;
 }
 
static void _hoc_Exp(void) {
  double _r;
   _r =  Exp (  *getarg(1) );
 hoc_retpushx(_r);
}
 
static int _ode_count(int _type){ hoc_execerror("sweeney", "cannot be used with CVODE"); return 0;}
 extern void nrn_update_ion_pointer(Symbol*, Datum*, int, int);
 static void _update_ion_pointer(Datum* _ppvar) {
   nrn_update_ion_pointer(_na_sym, _ppvar, 0, 0);
   nrn_update_ion_pointer(_na_sym, _ppvar, 1, 3);
   nrn_update_ion_pointer(_na_sym, _ppvar, 2, 4);
 }

static void initmodel() {
  int _i; double _save;_ninits++;
 _save = t;
 t = 0.0;
{
  h = h0;
  m = m0;
 {
   rates ( _threadargscomma_ v ) ;
   m = minf ;
   h = hinf ;
   }
  _sav_indep = t; t = _save;

}
}

static void nrn_init(NrnThread* _nt, _Memb_list* _ml, int _type){
Node *_nd; double _v; int* _ni; int _iml, _cntml;
#if CACHEVEC
    _ni = _ml->_nodeindices;
#endif
_cntml = _ml->_nodecount;
for (_iml = 0; _iml < _cntml; ++_iml) {
 _p = _ml->_data[_iml]; _ppvar = _ml->_pdata[_iml];
#if CACHEVEC
  if (use_cachevec) {
    _v = VEC_V(_ni[_iml]);
  }else
#endif
  {
    _nd = _ml->_nodelist[_iml];
    _v = NODEV(_nd);
  }
 v = _v;
  ena = _ion_ena;
 initmodel();
 }}

static double _nrn_current(double _v){double _current=0.;v=_v;{ {
   ina = gnabar * m * m * h * ( v - ena ) ;
   il = gl * ( v - el ) ;
   }
 _current += ina;
 _current += il;

} return _current;
}

static void nrn_cur(NrnThread* _nt, _Memb_list* _ml, int _type){
Node *_nd; int* _ni; double _rhs, _v; int _iml, _cntml;
#if CACHEVEC
    _ni = _ml->_nodeindices;
#endif
_cntml = _ml->_nodecount;
for (_iml = 0; _iml < _cntml; ++_iml) {
 _p = _ml->_data[_iml]; _ppvar = _ml->_pdata[_iml];
#if CACHEVEC
  if (use_cachevec) {
    _v = VEC_V(_ni[_iml]);
  }else
#endif
  {
    _nd = _ml->_nodelist[_iml];
    _v = NODEV(_nd);
  }
  ena = _ion_ena;
 _g = _nrn_current(_v + .001);
 	{ double _dina;
  _dina = ina;
 _rhs = _nrn_current(_v);
  _ion_dinadv += (_dina - ina)/.001 ;
 	}
 _g = (_g - _rhs)/.001;
  _ion_ina += ina ;
#if CACHEVEC
  if (use_cachevec) {
	VEC_RHS(_ni[_iml]) -= _rhs;
  }else
#endif
  {
	NODERHS(_nd) -= _rhs;
  }
 
}}

static void nrn_jacob(NrnThread* _nt, _Memb_list* _ml, int _type){
Node *_nd; int* _ni; int _iml, _cntml;
#if CACHEVEC
    _ni = _ml->_nodeindices;
#endif
_cntml = _ml->_nodecount;
for (_iml = 0; _iml < _cntml; ++_iml) {
 _p = _ml->_data[_iml];
#if CACHEVEC
  if (use_cachevec) {
	VEC_D(_ni[_iml]) += _g;
  }else
#endif
  {
     _nd = _ml->_nodelist[_iml];
	NODED(_nd) += _g;
  }
 
}}

static void nrn_state(NrnThread* _nt, _Memb_list* _ml, int _type){
Node *_nd; double _v = 0.0; int* _ni; int _iml, _cntml;
#if CACHEVEC
    _ni = _ml->_nodeindices;
#endif
_cntml = _ml->_nodecount;
for (_iml = 0; _iml < _cntml; ++_iml) {
 _p = _ml->_data[_iml]; _ppvar = _ml->_pdata[_iml];
 _nd = _ml->_nodelist[_iml];
#if CACHEVEC
  if (use_cachevec) {
    _v = VEC_V(_ni[_iml]);
  }else
#endif
  {
    _nd = _ml->_nodelist[_iml];
    _v = NODEV(_nd);
  }
 v=_v;
{
  ena = _ion_ena;
 { error =  states();
 if(error){fprintf(stderr,"at line 63 in file sweeney.mod:\n	SOLVE states\n"); nrn_complain(_p); abort_run(error);}
 } }}

}

static void terminal(){}

static void _initlists() {
 int _i; static int _first = 1;
  if (!_first) return;
_first = 0;
}

#if NMODL_TEXT
static const char* nmodl_filename = "/home/merolanr/repos/BME515_HW2/sweeney.mod";
static const char* nmodl_file_text = 
  ": TITLE sweeney.mod Sweeney Channel\n"
  ": Sweeney channel fast sodium channel for myelinated axon\n"
  ": set for a resting potential of -80 mV\n"
  "\n"
  "NEURON {\n"
  "	SUFFIX sweeney\n"
  "	USEION na READ ena WRITE ina\n"
  "	NONSPECIFIC_CURRENT il\n"
  "	RANGE gnabar, gl, el, ena, il\n"
  "	RANGE minf, hinf\n"
  "	RANGE mtau, htau\n"
  "	RANGE mexp, hexp\n"
  "}\n"
  "\n"
  "UNITS {\n"
  "	(mA) = (milliamp)\n"
  "	(mV) = (millivolt)\n"
  "}\n"
  "\n"
  "INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}\n"
  "\n"
  "PARAMETER {\n"
  "	v (mV)\n"
  "	celsius = 37 (degC)\n"
  "	dt (ms)\n"
  "	gnabar = 1.445 (mho/cm2)\n"
  "	gl = 0.128 (mho/cm2)\n"
  "	el = -80.01 (mV)\n"
  "	ena = 35.64 (mV)\n"
  "	vtraub = 0 (mV)\n"
  "	\n"
  "	amA = 49\n"
  "	amB = 126\n"
  "	amC = 0.363\n"
  "	amD = 5.3\n"
  "	bmA = 56.2\n"
  "	bmB = 4.17\n"
  "	ahA = 56\n"
  "	ahB = 15.6\n"
  "	bhA = 10\n"
  "	bhB = 74.5\n"
  "	bhC = 5\n"
  "	\n"
  "}\n"
  "\n"
  "STATE {\n"
  "	m\n"
  "	h\n"
  "}\n"
  "\n"
  "ASSIGNED {\n"
  "	ina (mA/cm2)\n"
  "	il (mA/cm2)\n"
  "	minf\n"
  "	hinf\n"
  "	mtau (ms)\n"
  "	htau (ms)\n"
  "	mexp\n"
  "	hexp\n"
  "}\n"
  "\n"
  "BREAKPOINT {\n"
  "	SOLVE states\n"
  "	ina = gnabar * m * m * h * (v - ena)\n"
  "	il = gl * (v - el)\n"
  "}\n"
  "\n"
  ": DERIVATIVE states {\n"
  ": 	rates(v) METHOD cnexp\n"
  ": 	m' = (minf-m)/mtau\n"
  ": 	h' = (hinf-h)/htau\n"
  ":}\n"
  "\n"
  "PROCEDURE states() {    : exact when v held constant\n"
  "	rates(v)\n"
  "	m = m + mexp * (minf - m)\n"
  "	h = h + hexp * (hinf - h)\n"
  "	VERBATIM\n"
  "	return 0;\n"
  "	ENDVERBATIM\n"
  "}\n"
  "\n"
  "UNITSOFF\n"
  "\n"
  "INITIAL {\n"
  "	rates(v)\n"
  "	m = minf\n"
  "	h = hinf\n"
  "}\n"
  "\n"
  "PROCEDURE rates(v) {\n"
  "\n"
  "	LOCAL v2, alpha, beta, sum\n"
  "	v2 = v-vtraub\n"
  "	\n"
  "	: m sodium activation system\n"
  "	alpha = (amB+amC*v2)/(1+Exp(-(amA+v2)/amD))\n"
  "	beta = (1/Exp((v2+bmA)/bmB))*(amB+amC*v2)/(1+Exp(-(amA+v2)/amD))\n"
  "	sum = alpha + beta\n"
  "	mtau = 1/sum\n"
  "	minf = alpha/sum\n"
  "	\n"
  "	: h sodium inactivation system\n"
  "	alpha = (1/Exp((v2+bhB)/bhC))*ahB/(1+Exp(-(v2+ahA)/bhA))\n"
  "	beta = ahB/(1+Exp(-(v2+ahA)/bhA))\n"
  "	sum = alpha + beta\n"
  "	htau = 1/sum\n"
  "	hinf = alpha/sum\n"
  "\n"
  "	mexp = 1 - Exp(-dt/mtau)\n"
  "	hexp = 1 - Exp(-dt/htau)\n"
  "}\n"
  "\n"
  "FUNCTION Exp(x) {\n"
  "	if (x<-100) {\n"
  "		Exp = 0\n"
  "	}else{\n"
  "		Exp = exp(x)\n"
  "	}\n"
  "}\n"
  "\n"
  "UNITSON\n"
  ;
#endif
