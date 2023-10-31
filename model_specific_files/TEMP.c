/***
   NAME
    TEMP.c
   NOTES
        Definition file for single cohort and population level model, adapted with temperature-dependent functions from Ohlberger et al. (2011).
        
        Used to model: cod, sprat, herring, rockling, mullet

        Allows for constant or changing, seasonal temperature (Pethybridge et al. 2013)

        Allows for constant or changing, seasonal resource (Pethybridge et al. 2013)

        Differences from Ohlberger et al. (2011) include:
        - Using parameter values (40 and 400) in temperature dependent equations 
          unlike in Ohlberger et al. (2011) which used 20 and 200 for ingestion
            See Table 2. a (bottom of table) for more info
        - Rm calculated with weight = bones + fat + gonads (Ohlberger et al. 2011)
          - We are using standarized weight (1 + QJ)*X instead; no maintenance costs for gonads (Kooijman, 2010)
            - Define either calculation with RM_CALCULATION directive

    The problem specific definition file for a size structured model including:

      - Zooplankton resource for SPECIES ~ FEEDING_LEVEL parameter (0-1)

    Terminology:

        X:      Irreversible mass or bone
        Y:      Reversible mass or fat
        M:      Mass or standard weight, (1+QJ)*X
        W:      Weight or total weight, (X+Y)

    Last revisions: OG AvL - 5 May 2022; TEMP BD - April 20 2023

***/

// General directives
#define IMMEDIATE_HATCHING      1               // 1: Cohort is filled at first day of spawngroup
#define EGGM_EQ_SDM             0               // 0: Following Ohlberger et al. 2011, no egg mortality; 1: egg/yolk-sac mortality follows size-dep. mortality
#define END_INV                 0               // 1: ForcedRunEnd if R0 > 1.5 or R0 < 0.5
#define R0_OUTPUT               2               // 0: None; 1: output of ReproOutputs[][]; 2: In addition, output of collapsed matrix
#define RM_CALCULATION          1               // 0: Calculate maintenace with total weight (Ohlberger et al. 2011); 1: Calculate with standarized weight (Kooijman, 2010)
#define TEMPERATURE             1               // 0: Constant temperature; 1: Seasonally varying temperature using temp_sin function 
#define RESOURCE                0               // 0: Constant resource; 1: Seasonally varying resource using resource_sin function 

// Species directives      
#define SPECIES_SINGLE          1               // 1: Single cohort (i.e. individual level), no feedback on environment (i.e. constant resource)
#define SPECIES_COHORT          0               // 1: Multiple cohorts (i.e. population level), feedback on environment
#define SPECIES_INVADE          0

#if (BIFURCATION == 0)
#undef  SPECIES_INVADE
#define SPECIES_INVADE          0
#endif

#include  "escbox.h"

#if(USE_MKL == 1)
#include  "mkl_vml.h"
#else

void    vdPowx(int n, double *x, double b, double *z)
{ 
  int   i;

  for (i=0; i<n; i++)
    z[i] = pow(x[i],b);

  return;
}
    
void    vdPowx_2(int n, double *x, double *b, double *z)
{ 
  int   i;

  for (i=0; i<n; i++)
    z[i] = pow(x[i],b[i]);

  return;
}

void    vdMul(int n, double x, double *y, double *z) 
{ 
  int   i;

  for (i=0; i<n; i++)
    z[i] = x*y[i];

  return;
}
   
void    vdMul_2(int n, double *x, double *y, double *z) 
{ 
  int   i;

  for (i=0; i<n; i++)
    z[i] = x[i]*y[i];

  return;
}

void    vdExp(int n, double *x, double *y)
{ 
  int   i;

  for (i=0; i<n; i++)
    y[i] = exp(x[i]);

  return;
}

void    vdLog(int n, double *x, double *y) // removed [i] and *
{
  int   i;

  for(i = 0; i<n;  i++)
    y[i] = log(x[i]); 

  return;
}

void    vdDiff(int n, double *x, double y, double *z) // 
{
  int   i;

  for(i = 0; i<n; i++)
    z[i] = x[i]-y; 

  return;
}

void    vdDiff_2(int n, double *x, double *y, double *z) // 
{
  int   i;

  for(i = 0; i<n; i++)
    z[i] = x[i]-y[i]; 

  return;
}

void    vdDiff_3(int n, double x, double *y, double *z) // 
{
  int   i;

  for(i = 0; i<n; i++)
    z[i] = x-y[i]; 

  return;
}

void    vdAdd(int n, double *x, double y, double *z) // removed * from double y
{
  int   i;

  for(i = 0; i<n; i++)
    z[i] = x[i]+y; 

  return;
}

void    vdDiv(int n, double x, double *y, double *z) //
{
  int   i;

  for(i=0; i<n; i++)
    z[i] = x/y[i];

  return;
}

void    vdDiv_2(int n, double *x, double *y, double *z) //
{
  int   i;

  for(i=0; i<n; i++)
    z[i] = x[i]/y[i];

  return;
}
#endif  

/*
 *======================================================================
 *
 *              LABELING ENVIRONMENT AND I-STATE VARIABLES
 *
 *======================================================================
 */

#define time           (env[0])

#define plankton       (env[1]) 
#define TotalEggs      (env[2])

#define SPECIES           0

#define age             (i_state( 0)) // second value of second row of .isf. First value is # individuals (refer to EBTmanual.pdf)
#define bone            (i_state( 1))
#define fat             (i_state( 2))
#define gonads          (i_state( 3))

// Constants not set in updateIDcards()
#define IDnbegin        (i_const( 0))
#define IDlifefec       (i_const( 1))
#define IDspawning      (i_const( 2))

// Constant for temporary storage
#define IDscratch       (i_const( 3))

// Size measures set in updateIDcards()
#define IDmass          (i_const( 4))
#define IDweight        (i_const( 5))
#define IDlength        (i_const( 6))
#define IDfatratio      (i_const( 7))

// Energetic rates set in updateIDcards()
#define IDingestPs      (i_const( 8))
#define IDmaint         (i_const( 9))
#define IDnet_energy    (i_const(10))
#define IDkappa         (i_const(11))

// Mortality & reproduction rates set in updateIDcards()
#define IDfecundity     (i_const(12))
#define IDmort          (i_const(13)) // Initial mortality, change 0.0957458243 to 0.0 in if you want no initial mortality

#define MILLI           0.001
#define MICRO           1.0E-6

#define MINCONDITION    1.0E-6
#define FOODTINY        1.0E-10

#define MAXCOHORTS      10000         
#define INVADEDENSITY   1000
#define INVADEPERIOD    30            
#define MINSURVIVAL     1.0E-9

#define DIV_1_400       0.0025        // Used in temperature equation

/*
 *======================================================================
 *
 *          DEFINING AND LABELLING CONSTANTS AND PARAMETERS
 *
 *======================================================================
 */

// Species-independent parameters
#define YEAR            parameter[ 0]
#define VOLBOT_RATIO    parameter[ 1]           // Volume : Bottom ratio
#define VOLUME          parameter[ 2]           // Reference volume
#define BOTTOM          (VOLUME/VOLBOT_RATIO)   // Benthic surface area

// Species parameters: those of its resource first
#define PARONE 3

#define TEMP            parameter[PARONE+ 0]   // User defined constant environmental temperature 
#define FEEDING_LEVEL   parameter[PARONE+ 1]   // User defined feeding level (0-1); 1 is amount of food for max growth

#define RPLANKTON       parameter[PARONE+ 2]   // Growth rate zooplankon
#define KPLANKTON       parameter[PARONE+ 3]   // Carrying capacity zooplankton (g*m^3)

#define SPAWNSET        parameter[PARONE+ 4]
#define SPAWNMAX        parameter[PARONE+ 5]
#define SPAWNPERIOD     parameter[PARONE+ 6]
#define SPAWNGROUP      parameter[PARONE+ 7]

#define SPAWNSTART      (SPAWNMAX-0.5*SPAWNPERIOD)
#define SPAWNEND        (SPAWNMAX+0.5*SPAWNPERIOD)

#define MUB             parameter[PARONE+ 8]   // Background mortality
#define MUSDC           parameter[PARONE+ 9]   // Size-dependent mortality constant
#define MUSDMASS        parameter[PARONE+10]   // Size-dependent mortality characteristic size
#define MUS             parameter[PARONE+11]   // Starvation mortality constant
#define EGGMORT         parameter[PARONE+12]   // Egg mortality
#define YOLKMORT        parameter[PARONE+13]   // Mortality of yolk-sac larvae

#define FSTARTW         parameter[PARONE+14]   // Size of start fishing vulnerability; Can be removed but fine for now
#define FHALFW          parameter[PARONE+15]   // Size of 50% fishing vulnerability; Can be removed but fine for now
#define FMORT           parameter[PARONE+16]   // Fishing mortality; Can be removed but fine for now

#define LWC             parameter[PARONE+17]   // Weight-length(!) scaling constant
#define LWE             parameter[PARONE+18]   // Weight-length scaling exponent

#define EGGPERIOD       parameter[PARONE+19]   // Duration of egg period
#define AGEFEEDING      parameter[PARONE+20]   // Age at first feeding
#define BIRTHWEIGHT     parameter[PARONE+21]   // Total weigth at birth
#define MATURELEN       parameter[PARONE+22]   // Maturation length

#define QS              parameter[PARONE+23]   // Starvation condition threshold
#define QJ              parameter[PARONE+24]   // Juvenile condition target
#define QA              parameter[PARONE+25]   // Adult condition target
#define QSPAWN          parameter[PARONE+26]   // Threshold condition for spawning

#define MAINTC          parameter[PARONE+27]   // Maintenance scaling constant
#define MAINTE          parameter[PARONE+28]   // Maintenance scaling exponent

#define DIGTIMEC        parameter[PARONE+29]   // Digestion time scaling constant
#define DIGTIMEE        parameter[PARONE+30]   // Digestion time scaling exponent

#define ARATEPMAX       parameter[PARONE+31]   // Maximum attack rate
#define ARATEPWOPT      parameter[PARONE+32]   // Mass at maximum attack
#define ARATEPEXP       parameter[PARONE+33]   // Maxumim attack rate exponent

#define CONVEFF         parameter[PARONE+34]   // Conversion efficiency
#define REPROEFF        parameter[PARONE+35]   // Reproductive conversion efficiency

#define GA_M            parameter[PARONE+36]   // Allometric scalar of Qm; metabolism
#define THETA_M         parameter[PARONE+37]   // Allometric exponent of Qm
#define YMAX_M          parameter[PARONE+38]   // Allometric scalar of Tmax
#define VMAX_M          parameter[PARONE+39]   // Allometric exponent of Tmax
#define YOPT_M          parameter[PARONE+40]   // Allometric scalar of Topt
#define VOPT_M          parameter[PARONE+41]   // Allometric exponent of Topt

#define GA_A            parameter[PARONE+42]   // Allometric scalar of Qa; ingestion
#define THETA_A         parameter[PARONE+43]   // Allometric exponent of Qa
#define YMAX_A          parameter[PARONE+44]   // Allometric scalar of Tmax
#define VMAX_A          parameter[PARONE+45]   // Allometric exponent of Tmax
#define YOPT_A          parameter[PARONE+46]   // Allometric scalar of Topt
#define VOPT_A          parameter[PARONE+47]   // Allometric exponent of Topt

#define INGESTC         parameter[PARONE+48]   // Simplified ingestion scaling constant
#define INGESTE         parameter[PARONE+49]   // Simplified ingestion scaling exponent

#define T_MEAN          parameter[PARONE+50]   // Temperature sin function mean value
#define T_A             parameter[PARONE+51]   // Temperature sin function amplitude
#define T_OMEGA         parameter[PARONE+52]   // Temperature sin function phase shift

#define R_MEAN          parameter[PARONE+53]   // Resource sin function mean value
#define R_A             parameter[PARONE+54]   // Resource sin function amplitude
#define R_OMEGA         parameter[PARONE+55]   // Resource sin function phase shift

/*
 *===========================================================================
 *
 *              DEFINING GLOBAL VARIABLES AND USER-IMPLEMENTED FUNCTIONS
 *
 *===========================================================================
 */

extern cohort           abs_tols[POPULATION_NR];

static void             UpdateIDcards(double *env, population *pop);

static int              TimeInYear;

int                     ReproCohorts[POPULATION_NR];
int                     NewCohort = -1, NewCCohort = -1; // Leaving NewCCohort for now, could change name or delete variable

double                  YOYsurvival[POPULATION_NR];
double                  SizeOYO[POPULATION_NR];
double                  TotalRepro[POPULATION_NR];
double                  MeanFecundity[POPULATION_NR];

#if ((SPECIES_INVADE == 1)) // || (COD_INVADE == 1))
double                  ReproOutputs[INVADEPERIOD][INVADEPERIOD], LastInvasionStart = -1;
#if ((R0_OUTPUT == 1) || (R0_OUTPUT == 2))
static FILE             *reprofile  = NULL;
#endif
#endif

#if (LENGTHCURVES > 0)
static double           BirthTimes[LENGTHCURVES];
#endif

#if (BIFURCATION == 1)
#define VARIANCES       2                       // Output variances as CVs
//#include  "MeasureBifstats.c"
#endif

double          mass[MAXCOHORTS];
double          weight[MAXCOHORTS];
double          length[MAXCOHORTS];
double          sdmmass[MAXCOHORTS];
double          sdmort[MAXCOHORTS];
double          partial_maint[MAXCOHORTS];
double          partial_htime[MAXCOHORTS];
double          partial_ingest[MAXCOHORTS];
double          partial_aratep[MAXCOHORTS];
double          mratio[MAXCOHORTS];
double          mratiomin[MAXCOHORTS];
double          aratepE[MAXCOHORTS];
double          aratep1[MAXCOHORTS];

// Define variables to solve for temperature dependent metabolic and ingestion rates
double         rm_weight[MAXCOHORTS];
double         temporary_Tmax_m[MAXCOHORTS], Tmax_m[MAXCOHORTS], temporary_Topt_m[MAXCOHORTS], Topt_m[MAXCOHORTS], temp_Q_m[MAXCOHORTS];
double         Q_m[MAXCOHORTS], log_Q_m[MAXCOHORTS], temp_diff_m[MAXCOHORTS], temp_diff_m_add_2[MAXCOHORTS], Y_m[MAXCOHORTS], W_m[MAXCOHORTS];
double         W_m_pow_2[MAXCOHORTS], div_40_Y_m[MAXCOHORTS], add_1_div_40_Y_m[MAXCOHORTS], pow_point_5_m[MAXCOHORTS];
double         one_pow_point_5_m[MAXCOHORTS], temporary_X_m[MAXCOHORTS], temporary_X_m_2[MAXCOHORTS], X_m[MAXCOHORTS];
double         V_num_m[MAXCOHORTS], V_den_m[MAXCOHORTS], V_m[MAXCOHORTS], V_pow_X_m[MAXCOHORTS], diff_1_V_m[MAXCOHORTS], X_mult_1_V_m[MAXCOHORTS];
double         exp_X_1_V_m[MAXCOHORTS], Rm[MAXCOHORTS];

double         temporary_Tmax_a[MAXCOHORTS], Tmax_a[MAXCOHORTS], temporary_Topt_a[MAXCOHORTS], Topt_a[MAXCOHORTS], temp_Q_a[MAXCOHORTS];
double         Q_a[MAXCOHORTS], log_Q_a[MAXCOHORTS], temp_diff_a[MAXCOHORTS], temp_diff_a_add_2[MAXCOHORTS], Y_a[MAXCOHORTS], W_a[MAXCOHORTS];
double         W_a_pow_2[MAXCOHORTS], div_40_Y_a[MAXCOHORTS], add_1_div_40_Y_a[MAXCOHORTS], pow_point_5_a[MAXCOHORTS];
double         one_pow_point_5_a[MAXCOHORTS], temporary_X_a[MAXCOHORTS], temporary_X_a_2[MAXCOHORTS], X_a[MAXCOHORTS];
double         V_num_a[MAXCOHORTS], V_den_a[MAXCOHORTS], V_a[MAXCOHORTS], V_pow_X_a[MAXCOHORTS], diff_1_V_a[MAXCOHORTS], X_mult_1_V_a[MAXCOHORTS];
double         exp_X_1_V_a[MAXCOHORTS], Ra[MAXCOHORTS];

double         ENV_TEMP, ENV_RESOURCE;

/*
 *===========================================================================
 *
 * UTILITY ROUTINES
 *
 *===========================================================================
 */

double  Sigmoid(double val, double low, double half)


{
  double                result = 0.0, nx;

  // Scale the dependent range between 0 and 3
  nx = 1.5*(val - low)/(half - low);

  if      (nx <= 0.0) result  = 0.0;
  else if (nx <= 1.0) result  = nx*nx*nx/6.0;
  else if (nx <= 2.0) result  = -3.0*nx/2.0 + 3.0*nx*nx/2.0 - nx*nx*nx/3.0 + 0.5;
  else if (nx <= 3.0) result  =  9.0*nx/2.0 - 3.0*nx*nx/2.0 + nx*nx*nx/6.0 - 3.5;
  else                result  = 1.0;

  return result;
}

// Changing temperature function
double temp_sin(double T)

{
  double temp;
  double PI = 3.141592653589793238463;

  temp = T_MEAN + T_A*sin(2*PI*(T + T_OMEGA) / 250); 
  return temp;
}

// Changing resource function 
double resource_sin(double T)

{
  double resource;
  double PI = 3.141592653589793238463;

  resource = R_MEAN + R_A*sin(2*PI*(T + R_OMEGA) / 250); 
  return resource;

}

/*===========================================================================*/
#if ((SPECIES_INVADE == 1)) // || (COD_INVADE == 1))

double DominantEigenvalue(double mat[INVADEPERIOD][INVADEPERIOD], int dim)

// Constructs the next-generation matrix and calculates its dominant eigenvalue

{
  int           ii, jj, Obs, Phase;
  int           n = dim, lda = dim, ldvl = dim, ldvr = dim, info, lwork;
  double        result = -1;
  double        wkopt;
  double        *work;
  double        repout[dim][dim], NextGenMat[dim][dim], wr[dim], wi[dim], vl[dim*dim], vr[dim*dim];
  extern void   dgeev_( char* jobvl, char* jobvr, int* n, double* a,
                        int* lda, double* wr, double* wi, double* vl, int* ldvl,
                        double* vr, int* ldvr, double* work, int* lwork, int* info );

  // Collapse the ReproOutputs matrix to a (dim x dim) matrix (dim=CyclePeriod)
  memset(repout[0], 0, (dim*dim)*sizeof(double));
  for (ii = 0; ii<INVADEPERIOD; ii++)
    {
      // Column index jj refers to reproduction at age jj YEARS
      // With 4 year cycle, first add reproduction at age 4, 8, 12, etc. to reproduction at age 0.
      // Similar for age 1, 2 and 3
      for (jj = dim; jj<INVADEPERIOD; jj++)             // Process columns CyclePeriod-INVADEPERIOD
        mat[ii][jj%dim] += mat[ii][jj];

      // Average rows 0, 4, 8, 12, etc. as they are measurements of the same quantity. Same for other rows.
      // Note that 4 is the 2nd observation (=1+ii/dim) and that ii%dim corresponds to year 0, 1 , 2 or 3
      // in 4-year cycle (the phase)
      Obs   = 1 + ii/dim;
      Phase = ii%dim;
      for (jj = 0; jj<dim; jj++)                        // Process columns 0-CyclePeriod
        UpdateStats(mat[ii][jj], &(repout[Phase][jj]), NULL, NULL, Obs);
    }
  // Matrix repout is condensed ReproOutputs matrix.
  // Rows correspond to phases in cycle, column 0 corresponds to total reproduction at age 0 & 4 & 8 & 12..years,
  // column 1 to total reproduction at age 1 & 5 & 9 & 13 years, etc.
  // NextGenMat should be a phase x phase matrix. Hence, shift the column values in every row to the right
  // row steps (wrap around column dimension)
  for (ii = 0; ii<dim; ii++)
    for (jj = 0; jj<dim; jj++)
      {
        Phase = (ii+jj)%dim;
        NextGenMat[ii][Phase] = repout[ii][jj];
      }

#if (R0_OUTPUT == 2)
  if (reprofile)
    {
      (void)fprintf(reprofile, "Next generation matrix of invaders     :\n");
      for (ii=0; ii<dim; ii++)
        {
          for (jj=0; jj<dim; jj++)
            PrettyPrint(reprofile, NextGenMat[ii][jj]);
          (void)fprintf(reprofile, "\n");
        }
      (void)fprintf(reprofile, "\n");
    }
#endif

  // If period = 1, we should have a single element, corresponding to the life time reproductive effort
  if (dim == 1) return NextGenMat[0][0];

  // Query and allocate the optimal workspace
  lwork = -1;
  dgeev_( "N", "N", &n, &(NextGenMat[0][0]), &lda, wr, wi, vl, &ldvl, vr, &ldvr, &wkopt, &lwork, &info );
  lwork = (int)wkopt;
  work = (double*)malloc( lwork*sizeof(double) );
  // Solve eigenproblem
  dgeev_( "N", "N", &n, &(NextGenMat[0][0]), &lda, wr, wi, vl, &ldvl, vr, &ldvr,   work, &lwork, &info );
  // Check for convergence
  if( info > 0 )
    Warning( "The algorithm failed to compute eigenvalues" );
  else if( info < 0 )
    Warning( "The algorithm to compute eigenvalues was called with an illegal input parameter" );
  else
    {
      result = wr[0];
      for (ii=1; ii<dim; ii++) result = max(result, wr[ii]);
    }

  return result;
}

#endif

/*
 *===========================================================================
 *
 * USER INITIALIZATION ROUTINE ALLOWS OPERATIONS ON INITIAL POPULATIONS
 *
 *===========================================================================
 */

void    UserInit(int argc, char **argv, double *env, population *pop)


{
  register int          i;

  /*  switch (argc)
    {
    case 4:
      C_MUB = atof(argv[3]);
    case 3:
      KCBENTHOS = atof(argv[2]);
    default:
      break;
    } */

  ReportNote("\n    %s%s",
#if (BIFURCATION == 1)
             "Bifurcation of the size-structured consumer-resource model",
#else
             "Dynamics of the size-structured consumer-resource model",
#endif // BIFURCATION
             "\n");

#if ((SPECIES_INVADE == 1))
#if (END_INV == 1)
  ReportNote("Bifurcation is stopped when R0 > 1.5 or R0 < 0.5");
#else
  ReportNote("Bifurcation is not stopped for any threshold R0");
#endif  // END_INV
#endif  // SPECIES_INVADE

#if (SPECIES_SINGLE == 1)
  ReportNote("Only a single cohort of species accounted for");

  for (i=0; i<cohort_no[SPECIES]; i++) pop[SPECIES][i][number]=-1.0;

  pop[SPECIES][0][number] = 1.0;
  pop[SPECIES][0][age]    = 0.0;
  pop[SPECIES][0][bone]   = 1/(1+QJ)*BIRTHWEIGHT;
  pop[SPECIES][0][fat]    = (QJ)/(1+QJ)*BIRTHWEIGHT;
  pop[SPECIES][0][gonads] = 0.0;

  for (i=0; i<I_CONST_DIM; i++) popIDcard[SPECIES][0][i] = 0.0;
  // Set plankton (i.e. food) to FEEDING_LEVEL; 
  plankton = FEEDING_LEVEL;
#endif

#if (SPECIES_COHORT == 1)
  ReportNote("Multiple cohorts accounted for");
#endif
       
  for (i = 0; i < cohort_no[SPECIES]; i++)
    {
      if (!(iszero(popIDcard[SPECIES][i][IDspawning]) || isequal(popIDcard[SPECIES][i][IDspawning], 1.0)))
        popIDcard[SPECIES][i][IDspawning] = (pop[SPECIES][i][gonads] > 0.0);

      popIDcard[SPECIES][i][IDnbegin] = max(popIDcard[SPECIES][i][IDnbegin], pop[SPECIES][i][number]);
    }


/*#if (BIFURCATION == 1)
//  initMeasureBifstats(argv[1], 0);
#endif */

  // We assume cohort_limit to be 1.0 and that all timing values are in terms of integer number of days
  cohort_limit = 1.0;

  // Fill BirthTimes array with -1.0 (i.e. available entries)
#if (LENGTHCURVES > 0)
  for (i=0; i<LENGTHCURVES; i++) BirthTimes[i] = -1.0;
#endif

  ReportNote("Minimum survival extinction: %E", MINSURVIVAL);

#if (EGGM_EQ_SDM == 1)
  ReportNote("egg/yolk-sac mortality follows size-dep. mortality");
  EGGMORT  = MUSDC*exp(-BIRTHWEIGHT/MUSDMASS);
  YOLKMORT = EGGMORT;
#endif

#if (EGGM_EQ_SDM == 0)
  ReportNote("NO egg/yolk-sac mortality");
#endif

  ReportNote("CHECK_EXT: %d (0 - no action; 1 - message; 2 - message and quit)", CHECK_EXTINCTION);

#if (IMMEDIATE_HATCHING == 1)
  ReportNote("Cohort is filled at first day of spawngroup");
#endif

#if (RM_CALCULATION == 0) 
  ReportNote("Rm calculated with weight = bones + fat + gonads (Ohlberger et al. 2011)");
#else     
  ReportNote("Rm calculated with standardized weight (Kooijman, 2010)");  
#endif  // Weight used in Rm calculation

#if (TEMPERATURE == 0)
  ReportNote("Constant environmental temperature");
#else
  ReportNote("Changing environmental temperature: mean temperature = %.1f, amplitude = %.1f, omega = %.1f", T_MEAN, T_A, T_OMEGA); 
#endif

#if (RESOURCE == 0)
  ReportNote("Constant environmental resource");
#else
  ReportNote("Changing environmental resource: mean resource = %.1f, amplitude = %.1f, omega = %.1f", R_MEAN, R_A, R_OMEGA); 
#endif

  SievePop();
  UpdateIDcards(env, pop);

  return;
}

/*
 *===========================================================================
 *
 *      SPECIFICATION OF THE NUMBER AND VALUES OF BOUNDARY POINTS
 *
 *===========================================================================
 */

void    SetBpointNo(double *env, population *pop, int *bpoint_no)

{
/*#if (BIFURCATION == 1)
 measureBifstats(env, pop, 3, 3);               // Measure periodicity in species egg numbers
#endif // BIFURCATION */

  bpoint_no[SPECIES] = 0; 

  TimeInYear = floor(fmod(time, YEAR) + MILLI);

#if (SPECIES_SINGLE != 1 )
  if (TotalEggs > MICRO)
    {
      int       Hatching, CreateCohort, i;
      double    Newborns;

      // Within spawn period? NB: cohort_limit = 1.0!!!
      Hatching     = ((TimeInYear - SPAWNSTART) >= -MILLI) && ((TimeInYear - SPAWNEND) < -(1.0-MILLI));

      if (Hatching)
        {
          // On a spawngroup boundary???
          CreateCohort = ((fmod(TimeInYear - SPAWNSTART, SPAWNGROUP) < MILLI) || 
                          (fmod(TimeInYear - SPAWNSTART, SPAWNGROUP) > SPAWNGROUP - MILLI));

          // Create cohort if both conditions apply
          if (CreateCohort)
            {
              // Compute cumulative newborn produced at the end of the coming spawn group period
              Newborns  = Sigmoid(TimeInYear + SPAWNGROUP, SPAWNSTART, SPAWNMAX);

              // Substract cumulative newborn produced up to now
              Newborns -= Sigmoid(TimeInYear, SPAWNSTART, SPAWNMAX);
              Newborns *= TotalEggs;

              NewCohort = AddCohorts(pop, SPECIES, 1);

#if (IMMEDIATE_HATCHING == 0)
              pop[SPECIES][NewCohort][number] = 0.0;
#else // IMMEDIATE_HATCHING == 1
              pop[SPECIES][NewCohort][number] = Newborns;
#endif
              pop[SPECIES][NewCohort][age]    = 0.0;
              pop[SPECIES][NewCohort][bone]   = 1/(1+QJ)*BIRTHWEIGHT;
              pop[SPECIES][NewCohort][fat]    = (QJ)/(1+QJ)*BIRTHWEIGHT;
              pop[SPECIES][NewCohort][gonads] = 0.0;

              for (i=0; i<I_CONST_DIM; i++) popIDcard[SPECIES][NewCohort][i] = 0.0;

              popIDcard[SPECIES][NewCohort][IDnbegin] = Newborns;
              
#if (LENGTHCURVES > 0)
              // Only track the first and last spawn group in a single year
              if (((TimeInYear - SPAWNSTART) < 0.5*SPAWNGROUP) || ((TimeInYear - SPAWNEND) > -1.5*SPAWNGROUP))
                for (i=0; i<LENGTHCURVES; i++)
                  {
                    if (iszero(BirthTimes[i]+1.0))
                      {
                        BirthTimes[i] = time;
                        break;
                      }
                  }
#endif
            }
        }
#if (IMMEDIATE_HATCHING == 0)
      else 
#endif
        NewCohort =  -1;
    }
#endif

  return;
} 


/*===========================================================================*/

void    SetBpoints(double *env, population *pop, population *bpoints)


{
    return;
}

/*
 *=====================================================================
 *
 *           SPECIFICATION OF DERIVATIVES
 *
 *=====================================================================
 */

void    Gradient(double *env,     population *pop,     population *ofs,
                 double *envgrad, population *popgrad, population *ofsgrad,
                 population *bpoints)

{
  register int          i;
  double                kappa = 0.0, net_energy = 0.0;
  double                tot_grazingsp = 0.0;
  double                mort = 0.0;

  UpdateIDcards(env, pop);

// #if (BIFURCATION == 1)
//   if (rk_level == 1) measureMinMax(env, pop);
// #endif

  // Species derivatives
  for (i=0; i<cohort_no[SPECIES]; i++)
    {
      popgrad[SPECIES][i][number] =  0.0;
      popgrad[SPECIES][i][age]    =  1.0;
      popgrad[SPECIES][i][bone]   =  0.0;
      popgrad[SPECIES][i][fat]    = -1.0;
      popgrad[SPECIES][i][gonads] =  0.0;
      if (popIDcard[SPECIES][i][IDfatratio] <= MINCONDITION) continue;

      mort       = popIDcard[SPECIES][i][IDmort];
      kappa      = popIDcard[SPECIES][i][IDkappa];
      net_energy = popIDcard[SPECIES][i][IDnet_energy];

      popgrad[SPECIES][i][number] = -mort*pop[SPECIES][i][number];
      if (iszero(popIDcard[SPECIES][i][IDspawning]))
        {
          popgrad[SPECIES][i][bone]       = kappa*net_energy;
          popgrad[SPECIES][i][fat]        = (1-kappa)*net_energy;
        }
      else
        {
          if (net_energy > 0.0)
            {
              popgrad[SPECIES][i][bone]   = kappa*net_energy;
              popgrad[SPECIES][i][gonads] = (1-kappa)*net_energy;
              popgrad[SPECIES][i][fat]    = 0.0;
            }
          else
            {
              if ((pop[SPECIES][i][fat]    > QS*pop[SPECIES][i][bone]) ||
                  (pop[SPECIES][i][gonads] < MICRO))
                {
                  popgrad[SPECIES][i][gonads] = 0.0;
                  popgrad[SPECIES][i][fat]    = net_energy;
                }
              else
                {
                  popgrad[SPECIES][i][gonads] = net_energy;
                  popgrad[SPECIES][i][fat]    = 0.0;
                }
            }
        }
#if ((SPECIES_SINGLE != 1) && (SPECIES_INVADE !=1))
      tot_grazingsp += popIDcard[SPECIES][i][IDingestPs]*pop[SPECIES][i][number];
#endif

    }

#if(IMMEDIATE_HATCHING == 0)
  double                nx, Hatching;

  if (NewCohort >= 0)
    {
      nx = 3.0*(fmod(time, YEAR) - SPAWNSTART)/SPAWNPERIOD;

      if ((nx <= 0.0) || (nx >= 3.0))   Hatching = 0.0;
      else if (nx <= 1.0)               Hatching =                 0.5*nx*nx;
      else if (nx <= 2.0)               Hatching = -1.5 + 3.0*nx -     nx*nx;
      else                              Hatching =  4.5 - 3.0*nx + 0.5*nx*nx;

      Hatching *= TotalEggs/(SPAWNPERIOD/3.0);

      mort  = popIDcard[SPECIES][NewCohort][IDmort];
#if ( COD_FEEDBACK & 0x0100 )
      mort += popIDcard[SPECIES][NewCohort][IDpiscmort];
#endif

      popgrad[SPECIES][NewCohort][number] = Hatching - mort*pop[SPECIES][NewCohort][number];
      popgrad[SPECIES][NewCohort][age]    = 0.5;
      popgrad[SPECIES][NewCohort][bone]   = 0.0;
      popgrad[SPECIES][NewCohort][fat]    = 0.0;
      popgrad[SPECIES][NewCohort][gonads] = 0.0;
    }
#endif


  // Environment derivatives
  #if (SPECIES_SINGLE == 1) // No feedback of single cohort on resource
  envgrad[0] = 1.0;
  envgrad[1] = 0.0; // No change (constant) in PLANKTON = FEEDING_LEVEL
  envgrad[2] = 0.0;
  #elif (SPECIES_COHORT == 1) // Feedback of cohorts on resource
  envgrad[0] = 1.0;
  envgrad[1] = RPLANKTON * (KPLANKTON - plankton) - tot_grazingsp;
  envgrad[2] = 0.0;
  #endif

  return;
}

/*
 *===========================================================================
 *
 *              SPECIFICATION OF BETWEEN COHORT CYCLE DYNAMICS
 *
 *===========================================================================
 */

void    InstantDynamics(double *env, population *pop, population *ofs)

{
  register int          i, pn;

#if ((LENGTHCURVES > 0))
  int                   j;
#endif

  if (FEEDING_LEVEL < FOODTINY) FEEDING_LEVEL = FOODTINY; // Check negative food

  UpdateIDcards(env, pop);

  for (pn=0; pn<POPULATION_NR; pn++)
    for (i=0; i<cohort_no[pn]; i++)             // remove cohorts with no fat
      {
        if ((popIDcard[pn][i][IDfatratio] > MINCONDITION) &&
            (pop[pn][i][number] > abs_tols[pn][0]) &&
            (pop[pn][i][number] > MINSURVIVAL*popIDcard[pn][i][IDnbegin])) continue;

        pop[pn][i][number] = 0.0;
#if ((LENGTHCURVES))
        switch (pn)
          {
            case 0:
#if (LENGTHCURVES > 0)
              for (j=0; j<LENGTHCURVES; j++)
                {
                  if (fabs(time - BirthTimes[j] - pop[pn][i][age]) < 0.5) // cohort_limit is 1.0
                    {
                      BirthTimes[j] = -2.0;
                      break;
                    }
                }
#endif
              break;
            case 1:
              break;
            }
#endif
      }

#if(IMMEDIATE_HATCHING == 0)
  // Ensure that the newborn cohort is not being discarded
  if (NewCohort >= 0)
    pop[SPECIES][NewCohort][number] = max(pop[SPECIES][NewCohort][number], (1.0 + MILLI)*abs_tols[SPECIES][0]);
#endif

  TimeInYear = floor(fmod(time, YEAR) + MILLI);

  // Species reproduction
  if (isequal(TimeInYear, SPAWNSET) ||
      ((iszero(SPAWNSET) || isequal(SPAWNSET, YEAR)) && isequal(TimeInYear, YEAR)))
    {
      for (i=0; i<cohort_no[SPECIES]; i++)
        {
        if (pop[SPECIES][i][number] <= abs_tols[SPECIES][0]) continue;

        // The following does not lead to cohort splitting due to skipped spawning
          if ((popIDcard[SPECIES][i][IDlength]   >= MATURELEN) &&
              (popIDcard[SPECIES][i][IDfatratio] >= QSPAWN))
            {
              popIDcard[SPECIES][i][IDspawning] = 1.0;
              pop[SPECIES][i][gonads]   = pop[SPECIES][i][fat];
              pop[SPECIES][i][fat]      = QSPAWN*pop[SPECIES][i][bone];
              pop[SPECIES][i][gonads]  -= pop[SPECIES][i][fat];
            }
          else
            popIDcard[SPECIES][i][IDspawning] = 0.0;
        }
    }

  if (isequal(TimeInYear, SPAWNSTART) ||
      ((iszero(SPAWNSTART) || isequal(SPAWNSTART, YEAR)) && isequal(TimeInYear, YEAR)))
    {
      double            fecundity;

      TotalEggs        = 0.0;
      ReproCohorts[SPECIES] = 0;
      TotalRepro[SPECIES] = 0.0;

      /*
       * NB: Note that the initial values for bone and fat of newly spawned
       * eggs is set to the bone and fat values of a first-feeding
       * larvae. Hence, the conversion coefficient REPROEFF should
       * account for the energy content of the yolk sac that is used to
       * maintain the egg and yolk-sac larvae during the period of
       * non-feeding
       */
      for (i=0; i<cohort_no[SPECIES]; i++)
        {
          if (pop[SPECIES][i][number] <= abs_tols[SPECIES][0]) continue;

          if (!iszero(popIDcard[SPECIES][i][IDspawning]))
            {
              fecundity = REPROEFF*max(pop[SPECIES][i][gonads], 0.0)/BIRTHWEIGHT;

              popIDcard[SPECIES][i][IDfecundity] = fecundity;
              popIDcard[SPECIES][i][IDlifefec]  += (fecundity*pop[SPECIES][i][number]/popIDcard[SPECIES][i][IDnbegin]);

#if (SPECIES_INVADE == 1)
              int       row, col;

              // row: the year number since start of invasions
              // col: age in years of reproducing cohort
              //      (1 is added because first possible reproduction is at age less than 1 year!)
              row = (int)floor((time - (LastInvasionStart + pop[SPECIES][i][age] - 0.5*cohort_limit))/YEAR);
              col = 1 + (int)floor((pop[SPECIES][i][age]-0.5*cohort_limit)/YEAR);
              ReproOutputs[row][col]          += (fecundity*pop[SPECIES][i][number]/INVADEDENSITY);
#else
              TotalEggs                      += fecundity*pop[SPECIES][i][number];
#endif

              ReproCohorts[SPECIES]++;
              TotalRepro[SPECIES]               += pop[SPECIES][i][number];
            }
          else
            popIDcard[SPECIES][i][IDfecundity]   = 0.0;

            popIDcard[SPECIES][i][IDspawning]    = 0.0;
            pop[SPECIES][i][gonads]              = 0.0;
        }

      if (TotalRepro[SPECIES]) MeanFecundity[SPECIES] = TotalEggs/TotalRepro[SPECIES];
      else MeanFecundity[SPECIES] = 0.0;

#if ((SPECIES_INVADE == 1) && (BIFURCATION == 1))
      double            dummy = fmod(time, BifPeriod);
      if (((dummy + 0.5*cohort_limit) > (BifPeriod - 2*INVADEPERIOD*YEAR)) &&
          ((dummy + 0.5*cohort_limit) < (BifPeriod -   INVADEPERIOD*YEAR)))
        {
          TotalEggs = INVADEDENSITY;

          if (LastInvasionStart < 0)
            {
              LastInvasionStart = time;
              memset((void *)ReproOutputs, 0, (INVADEPERIOD*INVADEPERIOD)*sizeof(double));
            }
        }
#endif
    }

  return;
} 


/*
 *===========================================================================
 *
 *                      SPECIFICATION OF OUTPUT VARIABLES
 *
 *===========================================================================
 */

void    DefineOutput(double *env, population *pop, double *output)

{
  register int          i, index;
  double                num, biom;

#if (LENGTHCURVES > 0)
  int                   j;
  double                Lengths[LENGTHCURVES];

  for (j=0; j<LENGTHCURVES; j++) Lengths[j] = DBL_MAX;
#endif


  UpdateIDcards(env, pop);

#if (BIFURCATION == 1)
  TimeInYear = (int)floor(fmod(time, YEAR) + MILLI);

  LabelState(SPECIES, "par. = %.4f   T = %6.0f years, %3d days", parameter[BifParIndex],
             floor((env[0]+0.1*cohort_limit)/YEAR), TimeInYear);
 #endif // BIFURCATION

#define OUTONE 1 // OUTONE 3

#if (SPECIES_SINGLE == 1)
  // Set output column two (first column is default set to time) to FEEDING_LEVEL; 
  output[ 0] = FEEDING_LEVEL; 
  output[OUTONE+ 0] = pop[SPECIES][0][number];
  output[OUTONE+ 1] = pop[SPECIES][0][bone];
  output[OUTONE+ 2] = pop[SPECIES][0][fat];
  output[OUTONE+ 3] = pop[SPECIES][0][gonads];
  output[OUTONE+ 4] = popIDcard[SPECIES][0][IDweight];
  output[OUTONE+ 5] = popIDcard[SPECIES][0][IDlength];
  output[OUTONE+ 6] = popIDcard[SPECIES][0][IDlifefec];
  output[OUTONE+ 7] = popIDcard[SPECIES][0][IDfatratio];
  output[OUTONE+ 8] = popIDcard[SPECIES][0][IDkappa];
  
  // Adjusted by Bass
  output[OUTONE+ 9] = time/YEAR; 
  output[OUTONE+10] = popIDcard[SPECIES][0][IDmaint];
  output[OUTONE+11] = popIDcard[SPECIES][0][IDingestPs];
  output[OUTONE+12] = popIDcard[SPECIES][0][IDnet_energy];
  output[OUTONE+13] = popIDcard[SPECIES][0][IDmort];
  output[OUTONE+14] = Rm[0];
  output[OUTONE+15] = Ra[0];
  output[OUTONE+16] = ENV_TEMP;
  output[OUTONE+17] = popIDcard[SPECIES][0][IDmass];
  output[OUTONE+18] = ENV_RESOURCE;

#elif(SPECIES_INVADE !=1) // 
  output[ 0] = plankton; 
  for (i=0; i<cohort_no[SPECIES]; i++)
    {
      num  = pop[SPECIES][i][number];
      biom = (pop[SPECIES][i][bone] + pop[SPECIES][i][fat] + pop[SPECIES][i][gonads])*pop[SPECIES][i][number];

      if (pop[SPECIES][i][age]< (YEAR - MILLI))           // YOY
        {
          output[OUTONE+ 0] += num; 
          output[OUTONE+ 1] += biom; 
        }
      else
        {
          output[OUTONE+ 2] += num;                     // 1+
          output[OUTONE+ 3] += biom;
        }
      if(popIDcard[SPECIES][i][IDlength] < MATURELEN)     // Juveniles
        {
          output[OUTONE+ 4] += num; 
          output[OUTONE+ 5] += biom; 
        }
      else                                              // Adults
        {
          output[OUTONE+ 6] += num; 
          output[OUTONE+ 7] += biom; 
        }

      output[OUTONE+ 8] += num;                         // Total population
      output[OUTONE+ 9] += biom; 

#if (LENGTHCURVES > 0)
      for (j=0; j<LENGTHCURVES; j++) 
        {
          if (fabs(time - BirthTimes[j] - pop[SPECIES][i][age]) < 0.5) // cohort_limit is 1.0
            {
              Lengths[j] = popIDcard[SPECIES][i][IDlength];
              break;
            }
        }
#endif
    }

  output[OUTONE+10] = popIDcard[SPECIES][0][IDlength]; 
  output[OUTONE+11] = popIDcard[SPECIES][0][IDlifefec]; // Same value as TotalEggs += fecundity*pop[SPECIES][i][number] on line 931
  output[OUTONE+12] = time/YEAR; 

index = 14;  
#if (LENGTHCURVES > 0)
  for (j=0; j<LENGTHCURVES; j++) 
    {
      output[index++]  = Lengths[j]; // 
      BirthTimes[j]   = max(BirthTimes[j], -1.0);
    }
#endif

#endif // Moved from line 1081 so single cohort output is only 19 outputs

/*#if (BIFURCATION == 1)
//  outputMeasureBifstats(env, NULL, OUTPUT_VAR_NR, 0);

#if ((SPECIES_INVADE == 1))
 static int             first = 1, direction = 0, fftIndex = 0, outputIndex = OUTONE+12;

#if (SPECIES_INVADE == 1)
 fftIndex    = 2;
 outputIndex = OUTONE+12;
#endif

 double dummy = fmod(env[0], BifPeriod);
 if ((env[0] > 0.5*cohort_limit) && BifPeriodEnd && ((dummy < 0.5*cohort_limit) || (dummy > (BifPeriod-0.5*cohort_limit))))
   {
     int CyclePeriod = (int)(rint(periodFFT[fftIndex]/YEAR) + MILLI);
#if ((R0_OUTPUT == 1) || (R0_OUTPUT == 2))
     if (first)
       {
         strcpy(fn, currun);                            // Open additional file for R0 output
         strcat(fn, ".R0.out");
         reprofile = fopen(fn, "a");
         if (!reprofile)
           fprintf(stderr, "Failed to open file %s\n", fn);
       }
     first = 0;
     if (reprofile)
       {
         int            ii, jj;

         (void)fprintf(reprofile, "Time                                   : ");
         PrettyPrint(reprofile, env[0]);
         (void)fprintf(reprofile, "\n");
         (void)fprintf(reprofile, "Current value of bifurcation parameter : ");
         PrettyPrint(reprofile, parameter[BifParIndex]);
         (void)fprintf(reprofile, "\n");
         (void)fprintf(reprofile, "Dominant period of population cycle    : ");
         PrettyPrint(reprofile, CyclePeriod);
         (void)fprintf(reprofile, "\n");
         (void)fprintf(reprofile, "\n");
         (void)fprintf(reprofile, "Full reproduction matrix of invaders   :\n");
         for (ii=0; ii<INVADEPERIOD; ii++)
           {
             for (jj=0; jj<INVADEPERIOD; jj++)
               PrettyPrint(reprofile, ReproOutputs[ii][jj]);
             (void)fprintf(reprofile, "\n");
           }
         (void)fprintf(reprofile, "\n");
       }
#endif
     if (!CyclePeriod) CyclePeriod = 1;
     output[outputIndex]  = DominantEigenvalue(ReproOutputs, CyclePeriod);
#if ((R0_OUTPUT == 1) || (R0_OUTPUT == 2))
     (void)fprintf(reprofile, "Dominant eigenvalue                    : ");
     PrettyPrint(reprofile, output[outputIndex]);
     (void)fprintf(reprofile, "\n\n\n");
#endif
#if (END_INV == 1)
     if (!direction)
       {
          if (output[outputIndex] < 1.0) direction = 1;
          else direction = -1;
       }
     else if (direction == 1)
       {
          if (output[outputIndex] > 1.5) ForcedRunEnd = 1;
       }
     else
       {
          if (output[outputIndex] < 0.5) ForcedRunEnd = 1;
       }
#endif

     // Reset invasion start time to -1
     LastInvasionStart = -1;

#if (DEBUG > 0)
     fprintf(stderr, "R0: %.10f\n", output[outputIndex]);
     ForcedRunEnd = 1;
#endif
   }

#else // i.e. #if ((SPECIES_INVADE != 1)

#endif // #if (SPECIES_INVADE == 1)

#endif // #if (BIFURCATION == 1)

#if (DEBUG > 0)
//  fprintf(stderr, "Time: %6.0f\tCohorts[0]: %5d\tCohorts[1]: %5d\n", env[0], cohort_no[0], cohort_no[1]);
#endif
*/
  return;
}


/*==============================================================================*/

static void     UpdateIDcards(double *env, population *pop)

{
  register int          i, s, c;
  double                fatratio, l_ratio;
  double                aratep, mass_encps, ingestps, ingestp_tot; //, mass_encpc, ingestpc
  /*double                arateb, mass_encb, ingestb;
  double                aratef, mass_encf, ingestf, aoptratef; */
  double                htime, gutspace;
  double                simple_ingest;
  double                maint, net_energy, kappa, q;
  double                mort;
  double                hatchfrac; //bentime, pisctime, 
  double                eggmort, yolkmort;
  

  // Initialize arrays for use in vdpowx
  for (i=0; i<cohort_no[SPECIES]; i++)
    {
      mass[i]      = pop[SPECIES][i][bone]*(1 + QJ);
      weight[i]    = pop[SPECIES][i][bone] + max(pop[SPECIES][i][fat] + pop[SPECIES][i][gonads], 0.0);
      mratio[i]    = mass[i]/ARATEPWOPT;
      mratiomin[i] = (1 - mratio[i]);
      if (MUSDC > 0) 
        sdmmass[i] = -mass[i]/MUSDMASS;

      popIDcard[SPECIES][i][IDmass]     = mass[i];
      popIDcard[SPECIES][i][IDweight]   = weight[i];
      popIDcard[SPECIES][i][IDfatratio] = (pop[SPECIES][i][fat] + pop[SPECIES][i][gonads])/pop[SPECIES][i][bone];

#if (RESOURCE == 0) // Constant 
      ENV_RESOURCE = plankton; // or set to FEEDING_LEVEL
#else 
      ENV_RESOURCE = resource_sin(time);
#endif
  
#if (TEMPERATURE == 0) // Calculate Rm and Ra with constant or varying temperature
      ENV_TEMP = TEMP;
#else 
      ENV_TEMP = temp_sin(time);
#endif

#if (RM_CALCULATION == 0) 
      rm_weight[i] = weight[i];  // Rm calculated with weight = bones + fat + gonads (Ohlberger et al. 2011)
#else     
      rm_weight[i] = mass[i];    // Rm calculated with standardized weight
#endif  
    }

  if (cohort_no[SPECIES])
    {
      #if (RM_CALCULATION == 0)      
      // Maintenance calculated with weight = bones + fat + gonads (Ohlberger et al. 2011)
      vdPowx(cohort_no[SPECIES], weight, MAINTE, partial_maint);  // pow(weight, MAINTE); solve for part of maintenance equation
#else   
      // Maintenance calculated with standardized weight  
      vdPowx(cohort_no[SPECIES], mass, MAINTE, partial_maint);   // pow(mass, MAINTE); solve for part of maintenance equation
#endif  

      vdPowx(cohort_no[SPECIES], mass, LWE, length);          // pow(mass, LWE)
      vdPowx(cohort_no[SPECIES], mass, DIGTIMEE, partial_htime);      // pow(mass, DIGTIMEE); solve for part of handling time equation, mass ~ standardized weight
      vdPowx(cohort_no[SPECIES], mass, INGESTE, partial_ingest); // pow(mass, INGESTC); solve for part of the simplifed ingestion equation, mass ~ standarized weight
      vdExp(cohort_no[SPECIES],  mratiomin, aratepE);           // exp(1-(mass/ARATEWOPT))
      vdMul_2(cohort_no[SPECIES], mratio, aratepE, aratep1);    // (mass/Wopt)*(exp(1-mass/Wopt))
      vdPowx(cohort_no[SPECIES], aratep1, ARATEPEXP, partial_aratep); // pow((mass/ARATEPWOPT)*exp(1-(mass/ARATEWOPT)), ARATEPEXP)      

      if (MUSDC > 0) 
        vdExp(cohort_no[SPECIES], sdmmass, sdmort);             // exp(-mass[i]/MUSDMASS) 

      // Steps to calculate Rm - temperature dependent metabolic factor 
      // Solve for Tmax_m which is function x + y

      vdPowx(cohort_no[SPECIES], rm_weight, VMAX_M, temporary_Tmax_m); // (x + y)^VMAX_M... temp_Tmax_m is z[i]
      vdMul(cohort_no[SPECIES], YMAX_M, temporary_Tmax_m, Tmax_m); // YMAX_M * (x + y)^VMAX_M... Tmax_m is z[i]
      
      // Solve for Topt_m which is function x + y
      vdPowx(cohort_no[SPECIES], rm_weight, VOPT_M, temporary_Topt_m); // (x + y)^VOPT_M...temp_Topt_m is z[i]
      vdMul(cohort_no[SPECIES], YOPT_M, temporary_Topt_m, Topt_m); // YOPT_M * (x + y)^VOPT_M... Topt_m is z[i]
      
      // Solve for Q_m = GA_M * (x + y)^THETA_M
      vdPowx(cohort_no[SPECIES], rm_weight, THETA_M, temp_Q_m); // x + y)^THETA_M... temp_Q_m is z[i]
      vdMul(cohort_no[SPECIES], GA_M, temp_Q_m, Q_m); // GA_M * (x + y)^THETA_M .. Q_m is z[i]
      
      // Solve for Y_m = (Tmax_m - Topt_m + 2) * log(Q_m)
      vdLog(cohort_no[SPECIES], Q_m, log_Q_m); // log(Q_m) = log(GA_M * (x + y)^THETA_M) ... log_Q_m is y[i]
      vdDiff_2(cohort_no[SPECIES], Tmax_m, Topt_m, temp_diff_m); // Tmax_m - Topt_m ...temp_diff_m is z[i]
      vdAdd(cohort_no[SPECIES], temp_diff_m, 2, temp_diff_m_add_2); // Tmax_m - Topt_m + 2 ... temp_diff_m_add_2 is z[i]
      vdMul_2(cohort_no[SPECIES], temp_diff_m_add_2, log_Q_m, Y_m); // (Tmax_m - Topt_m + 2) * log(Q_m) ... Y_m is z[i]

      // Solve for W_m = (Tmax_m - Topt_m) * log(Q_m)
      vdMul_2(cohort_no[SPECIES], temp_diff_m, log_Q_m, W_m); // (Tmax_m - Topt_m) * log(Q_m) .. W_m[i]

      // Solve for X_m = W_m^2 * (1/400) * (1 + (1 + (40 / Y_m))^0.5)^2
      vdPowx(cohort_no[SPECIES], W_m, 2, W_m_pow_2); // W_m^2 ...W_m_pow_2 is z[i]
      vdDiv(cohort_no[SPECIES], 40, Y_m, div_40_Y_m); // 40 / Y_m ...Div_40_Y_m is z[i]
      vdAdd(cohort_no[SPECIES], div_40_Y_m, 1, add_1_div_40_Y_m); // 1 + 40/Y_m   ..add_1_div_40_Y_m is z[i]
      vdPowx(cohort_no[SPECIES], add_1_div_40_Y_m, 0.5, pow_point_5_m); // (1 + 40/Y_m)^0.5 ...pow_point_5_m is z[i]
      vdAdd(cohort_no[SPECIES], pow_point_5_m, 1, one_pow_point_5_m); // 1 + (1 + 40/Y_m)^0.5 ...one_pow_point_5_m is z[i]
      vdPowx(cohort_no[SPECIES], one_pow_point_5_m, 2, temporary_X_m); // (1 + (1 + 40/Y_m)^0.5)^2 ..temporary_X_m is z[i]
      vdMul(cohort_no[SPECIES], DIV_1_400, temporary_X_m, temporary_X_m_2); // (1/400) * (1 + (1 + 40/Y_m)^0.5)^2 ...temporary_X_m_2 is z[i]
      vdMul_2(cohort_no[SPECIES], W_m_pow_2, temporary_X_m_2, X_m); // W_m^2 * (1/400) * (1 + (1 + (40 / Y_m))^0.5)^2... X_m is z[i]

      // Solve for V_m = (Tmax_m - TEMP) / (Tmax_m - Topt_m); TEMP is our specific environmental temperature set in parameter file
      vdDiff(cohort_no[SPECIES], Tmax_m, ENV_TEMP, V_num_m); // Tmax_m - TEMP ...V_num_m is z[i]
      vdDiff_2(cohort_no[SPECIES], Tmax_m, Topt_m, V_den_m); // Tmax_m - Topt_m ...V_den_m is z[i]
      vdDiv_2(cohort_no[SPECIES], V_num_m, V_den_m, V_m); // (Tmax_m - temp) / (Tmax_m - Topt_m) ...V_m is z[i]

      // Solve for Rm - temperature dependent metabolic factor 
      // Rm = V_m^(X_m) * (exp(X_m * 1 - V_m))) 
      vdPowx_2(cohort_no[SPECIES], V_m, X_m, V_pow_X_m); // V_m^X_m   ... V_pow_X_m is z[i]
      vdDiff_3(cohort_no[SPECIES], 1, V_m, diff_1_V_m); // 1 - V_m  ... diff_1_V_m is z[i]
      vdMul_2(cohort_no[SPECIES], X_m, diff_1_V_m, X_mult_1_V_m); // X_m * (1 - V_m)  ...X_mult_1_V_m is z[i]
      vdExp(cohort_no[SPECIES], X_mult_1_V_m, exp_X_1_V_m); // exp(X_m * (1 - V_m)) ...exp_X_1_V_m is z[i]
      vdMul_2(cohort_no[SPECIES], V_pow_X_m, exp_X_1_V_m, Rm); // V^(X) * (exp(X * (1 - V))) ... Rm is z[i]


      // Steps to calculate Ra - temperature dependent intake factor 
      // Solve for Tmax_a which is function x(1+qj)
      vdPowx(cohort_no[SPECIES], mass, VMAX_A, temporary_Tmax_a); // (standardized mass)^VMAX_a... temp_Tmax_a is z[i]
      vdMul(cohort_no[SPECIES], YMAX_A, temporary_Tmax_a, Tmax_a); // YMAX_a * (standardized mass)^VMAX_a... Tmax_a is z[i]
      
      // Solve for Topt_a which is function x + y
      vdPowx(cohort_no[SPECIES], mass, VOPT_A, temporary_Topt_a); // (standardized mass)^VOPT_a...temp_Topt_a is z[i]
      vdMul(cohort_no[SPECIES], YOPT_A, temporary_Topt_a, Topt_a); // YOPT_a * (standardized mass)^VOPT_a... Topt_a is z[i]
      
      // Solve for Q_a = GA_a * (x + y)^THETA_a
      vdPowx(cohort_no[SPECIES], mass, THETA_A, temp_Q_a); // standardized mass)^THETA_a... temp_Q_a is z[i]
      vdMul(cohort_no[SPECIES], GA_A, temp_Q_a, Q_a); // GA_a * (standardized mass)^THETA_a .. Q_a is z[i]
      
      // Solve for Y_a = (Tmax_a - Topt_a + 2) * log(Q_a)
      vdLog(cohort_no[SPECIES], Q_a, log_Q_a); // log(Q_a) = log(GA_a * (standardized mass)^THETA_a) ... log_Q_a is y[i]
      vdDiff_2(cohort_no[SPECIES], Tmax_a, Topt_a, temp_diff_a); // Tmax_a - Topt_a ...temp_diff_a is z[i]
      vdAdd(cohort_no[SPECIES], temp_diff_a, 2, temp_diff_a_add_2); // Tmax_a - Topt_a + 2 ... temp_diff_a_add_2 is z[i]
      vdMul_2(cohort_no[SPECIES], temp_diff_a_add_2, log_Q_a, Y_a); // (Tmax_a - Topt_a + 2) * log(Q_a) ... Y_a is z[i]

      // Solve for W_a = (Tmax_a - Topt_a) * log(Q_a)
      vdMul_2(cohort_no[SPECIES], temp_diff_a, log_Q_a, W_a); // (Tmax_a - Topt_a) * log(Q_a) .. W_a[i]

      // Solve for X_a = W_a^2 * (1/400) * (1 + (1 + (40 / Y_a))^0.5)^2
      vdPowx(cohort_no[SPECIES], W_a, 2, W_a_pow_2); // W_a^2 ...W_a_pow_2 is z[i]
      vdDiv(cohort_no[SPECIES], 40, Y_a, div_40_Y_a); // 40 / Y_a ...Div_40_Y_a is z[i]
      vdAdd(cohort_no[SPECIES], div_40_Y_a, 1, add_1_div_40_Y_a); // 1 + 40/Y_a   ..add_1_div_40_Y_a is z[i]
      vdPowx(cohort_no[SPECIES], add_1_div_40_Y_a, 0.5, pow_point_5_a); // (1 + 40/Y_a)^0.5 ...pow_point_5_a is z[i]
      vdAdd(cohort_no[SPECIES], pow_point_5_a, 1, one_pow_point_5_a); // 1 + (1 + 40/Y_a)^0.5 ...one_pow_point_5_a is z[i]
      vdPowx(cohort_no[SPECIES], one_pow_point_5_a, 2, temporary_X_a); // (1 + (1 + 40/Y_a)^0.5)^2 ..temporary_X_a is z[i]
      vdMul(cohort_no[SPECIES], DIV_1_400, temporary_X_a, temporary_X_a_2); // (1/400) * (1 + (1 + 40/Y_a)^0.5)^2 ...temporary_X_a_2 is z[i]
      vdMul_2(cohort_no[SPECIES], W_a_pow_2, temporary_X_a_2, X_a); // W_a^2 * (1/400) * (1 + (1 + (40 / Y_a))^0.5)^2... X_a is z[i]

      // Solve for V_a = (Tmax_a - TEMP) / (Tmax_a - Topt_a); TEMP is our specific environmental temperature set in parameter file
      vdDiff(cohort_no[SPECIES], Tmax_a, ENV_TEMP, V_num_a); // Tmax_a - TEMP ...V_num_a is z[i]
      vdDiff_2(cohort_no[SPECIES], Tmax_a, Topt_a, V_den_a); // Tmax_a - Topt_a ...V_den_a is z[i]
      vdDiv_2(cohort_no[SPECIES], V_num_a, V_den_a, V_a); // (Tmax_a - temp) / (Tmax_a - Topt_a) ...V_a is z[i]

      // Solve for Ra - temperature dependent intake factor 
      // Ra = V_a^(X_a) * (exp(X_a * 1 - V_a)))
      vdPowx_2(cohort_no[SPECIES], V_a, X_a, V_pow_X_a); // V_a^X_a   ... V_pow_X_a is z[i]
      vdDiff_3(cohort_no[SPECIES], 1, V_a, diff_1_V_a); // 1 - V_a  ... diff_1_V_a is z[i]
      vdMul_2(cohort_no[SPECIES], X_a, diff_1_V_a, X_mult_1_V_a); // X_a * (1 - V_a)  ...X_ault_1_V_a is z[i]
      vdExp(cohort_no[SPECIES], X_mult_1_V_a, exp_X_1_V_a); // exp(X_a * (1 - V_a)) ...exp_X_1_V_a is z[i]
      vdMul_2(cohort_no[SPECIES], V_pow_X_a, exp_X_1_V_a, Ra); // V^(X) * (exp(X * (1 - V))) ... Ra is z[i]

    }

#if (EGGM_EQ_SDM == 1)
  eggmort  = MUSDC*exp(-BIRTHWEIGHT/MUSDMASS);
  yolkmort = eggmort;
#else
  eggmort  = EGGMORT;
  yolkmort = YOLKMORT;
#endif
  // Update IDcards for species
  for (i=0; i<cohort_no[SPECIES]; i++)
    {
      length[i] *= LWC;
      fatratio    = popIDcard[SPECIES][i][IDfatratio];
      mort        = MUB;

      htime       = DIGTIMEC*partial_htime[i]/Ra[i]; // handling time = DIGTIMEC * MASS ^ DIGTIMEE * 1/Ra; Ra temp. adjust factor
      simple_ingest = INGESTC*partial_ingest[i]*Ra[i]; // simple_ingest = INGESTC * MASS ^ INGESTE * Ra 
      maint       = Rm[i]*MAINTC*partial_maint[i]; // maintenance;  Rm is temperature adjustment factor
      net_energy  = -maint;
      
      ingestps    = 0.0;
      kappa       = 0.0;
      aratep      = 0.0;
      
      if (pop[SPECIES][i][age] < EGGPERIOD)
        {
          // Species mortality
          mort       = eggmort;
          maint      = 0.0;
          net_energy = 0.0;
        }
      else if (pop[SPECIES][i][age] < AGEFEEDING)
        {
          // SPECIES mortality
          hatchfrac  = Sigmoid(pop[SPECIES][i][age], EGGPERIOD, EGGPERIOD+0.5*SPAWNGROUP);
          mort       = hatchfrac*yolkmort + (1-hatchfrac)*eggmort;
          maint      = 0.0;
          net_energy = 0.0;
        }
      else if (fatratio > MINCONDITION) 
        {
          if (MUSDC > 0)  // Size-dependent mortality 
            mort += MUSDC*sdmort[i];

          if (fatratio < QS) // Same condition as y < qstarv * x 
            mort += MUS*(QS/fatratio - 1.0); // s * (qstarv * x / y - 1)

          // We are including background mortality rate; specified as constant parameter
           // mort += MUB*pop[SPECIES][i][number]; // Already accounted for on line 1413

          // Size-selective fishing mortality 
          // mort   += FMORT*Sigmoid(mass[i], FSTARTW, FHALFW);

#if(SPECIES_SINGLE == 1)
          ingestps = ENV_RESOURCE*simple_ingest; // Ra accounted for in line 1358
          //ingestps = FEEDING_LEVEL*simple_ingest; // Ra accounted for in line 1358
          //ingestps   = FEEDING_LEVEL/htime; // 1/H*feedinglevel = (1/d1*w^d2)*Ra*feedinglevel; Ra accounted for in line (1356) 
          //ingestps = aratep*Kplankton/(1 + aratep*htime*Kplankton);

#elif(SPECIES_COHORT == 1)
          aratep     = ARATEPMAX*partial_aratep[i]*Ra[i];
          ingestps   = aratep*plankton/(1 + aratep*htime*plankton); 
#endif

          net_energy = (CONVEFF*ingestps - maint); // net_energy = (ke*(1/d1*w^d2)*Ra[i]*FEEDING_LEVEL) - (Rm[i]*MAINTC*maint[i])

          if (length[i] < MATURELEN) q = QJ; // Specifies energy allocation for juveniles or adults
          else q = QA;

          kappa = 0.0;
          if (net_energy > 0.0)
            {
              if (fatratio > q) // fat ratio used for these equations, same as our equations
                kappa = 1/(q+1);
              else
                // Default recovery formula changed to a quadratic dependence on fatratio (below), such that recovery in
                // reversible mass is faster
                kappa = fatratio*fatratio/((1+q)*q*q);
            }

          hatchfrac   = Sigmoid(pop[SPECIES][i][age], AGEFEEDING, AGEFEEDING+0.5*SPAWNGROUP);
          ingestps   *= hatchfrac;
          net_energy *= hatchfrac;
          maint      *= hatchfrac;
          mort       *= hatchfrac;
          mort       += (1.0 - hatchfrac)*yolkmort;
        }

      popIDcard[SPECIES][i][IDlength]     = length[i];
      popIDcard[SPECIES][i][IDmaint]      = maint;
      popIDcard[SPECIES][i][IDingestPs]   = ingestps;
      popIDcard[SPECIES][i][IDnet_energy] = net_energy;
      popIDcard[SPECIES][i][IDkappa]      = kappa;
      popIDcard[SPECIES][i][IDmort]       = mort;
    }

  return;
}  



/*=============================================================================*/
