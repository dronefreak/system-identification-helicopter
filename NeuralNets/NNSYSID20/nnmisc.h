/*
 *     HEADER FILE FOR NNMISC.C
 *
 *     Magnus Norgaard, IAU, Technical University of Denmark
 *     LastEditDate: Jan. 16, 2000.
 */  

/*
 *     Default values for training algorithms
 */  
#define TRDINFOLEVEL 0

/* Termination values */
#define TRDMAXITER 500
#define TRDCRITMIN 0
#define TRDCRITTERM 1e-7
#define TRDGRADTERM 1e-4
#define TRDPARAMTERM 1e-3

/* Weight decay */
#define TRDD 0

/* Levenberg-Marquardt parameters */
#define TRDLAMBDA 1

/* For recurrent nets */
#define TRDSKIP 0


/* Declaration of the "trparms" data-type. */
typedef struct {           
   int infolevel;
   int maxiter;
   double critmin;
   double critterm;
   double gradterm;
   double paramterm;
   matrix *D;
   double lambda;
   int skip;
} trparmstruct;


/*
 *     FUNCTION PROTOTYPES
 */  
void marqc(matrix**, int*, double*, matrix*, matrix*, matrix*, matrix*, matrix*,\
 trparmstruct*);
matrix* mat2sm(const mxArray*);
matrix* matstring2sm(const mxArray*);
void sm2mat(mxArray*, matrix*);
void mcopyi(matrix*, matrix*, matrix*, matrix*, matrix*, matrix*);
void mtanh(matrix*, matrix*, matrix*, matrix*, matrix*, matrix*);
matrix *neuvector(matrix*, int, char);
double sse(matrix*);
void deltatanh(matrix*, matrix*, matrix*, matrix*);
void mmultr1(matrix*, matrix*, matrix*);
void mmultr2(matrix*, matrix*, matrix*);
double sprod3(matrix*, matrix*);
matrix* makeindexvec(int, int, int);
void m2vreshape(matrix*, int, matrix*);
void m2vreshape2(matrix*, int, matrix*);
void v2mreshape(matrix*, matrix*, int);
void choldc(matrix*, matrix*);
void cholsl(matrix*, matrix*, matrix*);
void psi1(matrix*, int, matrix*, int, matrix*);
void psi2(matrix*, int, matrix*, int, double, matrix*);
void psi3(matrix*, matrix*, int, double);
void psi4(matrix*, int, matrix*, int, matrix*, matrix*);
void psi5(matrix*, int, matrix*, int, matrix*, matrix*, matrix*);
void vcopyi(matrix*, matrix*, int, matrix*, matrix*,int);
void vtanh(matrix*, matrix*, int, matrix*, matrix*, int);
void mvmul(matrix*, matrix*, matrix*, int);
