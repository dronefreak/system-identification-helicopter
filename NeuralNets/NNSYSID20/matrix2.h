/************************************************************************
*
*  Title     : matrix.c
*
*  Function  : header file for the matrix library.
*
*  Dependen-
*  cies      : matrix.c
*
*  Made by   : Steffen Torp & Magnus Noergaard, Servolaboratoriet DTU.
*
************************************************************************/

#define ON 1
#define OFF 0
#define RUNCHK ON             /* Run-time checks switched on/off. */

/* Inline functions with similar output as the library functions listed below. */
/* ( No run-time checks is performed, when inline functions are used )         */

#define nof_rows(ptm)                      (ptm->row)                           /* See getrows */
#define nof_cols(ptm)                      (ptm->col)                           /* See getcols */
#define vec_len(ptv)                       (ptv->row+ptv->col-1)                /* See length  */
#define get_val(ptm,row_pos,col_pos)       (ptm->mat[row_pos][col_pos])         /* See mget    */
#define put_val(ptm,row_pos,col_pos,value) ((ptm->mat[row_pos][col_pos])=value) /* See mput    */
#define rvget(ptv,element)                 (ptv->mat[0][element])               /* See vget    */
#define cvget(ptv,element)                 (ptv->mat[element][0])               /* See vget    */
#define rvput(ptv,element,value)           ((ptv->mat[0][element])=value)       /* See vput    */
#define cvput(ptv,element,value)           ((ptv->mat[element][0])=value)       /* See vput    */


/* Declaration of the "abstract" data-type. */

typedef struct {           /* Matrix structure for C library  */
	int row;               /* These are meant to be "private" */
	int col;               /* and should only be accessed via */
	double **mat;          /* the "member functions" below.   */
} matrix;

/* Declaration of the "member functions".   */

matrix *mmake( int, int );
matrix *mload( char*, char* );
double  sload( char*, char* );
matrix *minput1( void );
void mfree( matrix* );
void msave( matrix*, char*, char* );
void minput2( matrix* );
void mprint( matrix* );
void merror( char* );
void mtrans( matrix*, matrix* );
void minit( matrix* );
void minitx( matrix*, double);
void mrand( matrix* );
matrix *mfind( matrix*, double);
matrix *mnofind( matrix*, double);
void mdiag( matrix*, double );
void munit( matrix* );
void mset( matrix*, matrix* );
void madd( matrix*, matrix*, matrix* );
void msub( matrix*, matrix*, matrix* );
void mmul( matrix*, matrix*, matrix* );
void smul( matrix*, matrix*, double );
void mput( matrix*, int, int, double );
void vput( matrix*, int, double );
double mget( matrix*, int, int );
double vget( matrix*, int );
double trace( matrix* );
double sprod( matrix*, matrix* );
double sprod2( matrix*, matrix*, int, int, int, int );
void addcols( matrix*, matrix*, matrix* );
void addrows( matrix*, matrix*, matrix* );
void concat( matrix*, matrix*, matrix*, int, int, int, int );
void mat2mat( matrix*, int, int, matrix* );
void shift( matrix*, double );
void submat( matrix*, matrix*, int, int, int, int );
void subvec( matrix*, matrix*, int, int);
int getrows( matrix* );
int getcols( matrix* );
int length( matrix* );
int mscmp(matrix*, char*);
