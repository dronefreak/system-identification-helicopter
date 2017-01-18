#include <stdio.h>
#include <math.h>
#include "mex.h"
#include "matrix2.h"
#include "nnmisc.h"

#undef PI
/*********************************************************************************
 *                                                                               *
 *    MARQC                                                                      *
 *    -----                                                                      *
 *                                                                               *
 *    This is a C-version of the Matlab function marq.                           *
 *    Type 'help marq' from Matlab for information on                            *
 *    how to call this function.                                                 *
 *                                                                               *
 *                                                                               *
 *    Programmed by: Magnus Norgaard                                             *
 *    LastEditDate : Jan. 21, 2000                                               *
 *                                                                               *
 *********************************************************************************/
void marqc(matrix **PI_vectorpp, int *iter, double *lam,\
       matrix *NetDef, matrix *W1, matrix *W2, matrix *PHII,\
        matrix *Y, trparmstruct *trparms)
{
/*
-----------------------------------------------------------------------------------
---------------              VARIABLE DECLARATIONS                    -------------
-----------------------------------------------------------------------------------
*/
register i, j, k;
int iteration, outputs, N, Nout, layers, hidden, inputs;
int parameters1, parameters2, parameters, reduced, index1, ii, jj;
int lhids, hhids, louts, houts, index11;
double critdif, gradmax, paramdif, lambda, SSE, SSE_new, PI, PI_new, L, tmp1; 
double *ptm1, *ptm2, sum, lambda_old;
char dw;
matrix *L_hidden, *H_hidden, *L_output, *H_output, *h1, *h2, *y1, *y2;
matrix *E, *E_new, *E_vector, *E_new_vector, *W1_new, *W2_new, *PHI, *D, *Dtmp;
matrix *PI_vector, *tmp, *Htmp;
matrix *theta, *thtmp, *theta_index, *theta_red, *theta_red_new, *PSI, *G, *H, *h;
matrix *all, *index0, *index7, *onesvec, *tmp0, *tmp2, *tmp3, *index, *index2;

/*
-----------------------------------------------------------------------------------
---------------             NETWORK INITIALIZATIONS                   -------------
-----------------------------------------------------------------------------------
 */
outputs   = getrows(Y);                  /* # of outputs                         */
N         = getcols(Y);                  /* # of data                            */
hidden    = getrows(W1);                 /* # of hidden units                    */
inputs    = getcols(W1);                 /* Number of inputs                     */
inputs    = inputs-1;                 
Nout      = N*outputs;                   /* N*outputs                            */ 
h1        = mmake(hidden,N);             /* Argument to hidden layer act. fcts   */
h2        = mmake(outputs,N);            /* Argument to hidden layer act. fcts   */
iteration = 1;                           /* Initialize iteration counter         */
dw        = 1;                           /* Flag telling that the weights are new*/
onesvec   = mmake(1,N);                  /* Vector of all ones                   */
minitx(onesvec,1.0);
y1        = mmake(hidden+1,N); minit(y1);/* Hidden layer outputs                 */
mat2mat(y1,hidden,0,onesvec);            /* Add a row of ones (bias to outputs)  */
y2        = mmake(outputs,N); minit(y2); /* Output layer output                  */
E         = mmake(outputs,N);            /* Prediction error matrix              */
E_new     = mmake(outputs,N);            /* A priori E                           */
E_vector  = mmake(N*outputs,1);          /* Prediction error vector              */
E_new_vector = mmake(N*outputs,1);       /* A priori E_vector                    */
PHI       = mmake(inputs+1,N);           /* Matrix of input vectors (incl. bias) */
addrows(PHI,PHII,onesvec);
parameters1= hidden*(inputs+1);          /* # of input-to-hidden weights         */
parameters2= outputs*(hidden+1);         /* # of hidden-to-output weights        */
parameters = parameters1+parameters2;    /* Total # of weights                   */
theta      = mmake(parameters,1);        /* Vector containing all weights        */
m2vreshape(theta,0,W2);
m2vreshape(theta,parameters2,W1);
thtmp      = mnofind(theta,0.0);         /* Find non-zero entries in theta       */
reduced    = getrows(thtmp);             /* # of non-zero elements               */
theta_index = mmake(reduced,1);          /* Indices to weights <> 0              */
submat(theta_index,thtmp,0,reduced-1,0,0);
theta_red = mmake(reduced,1);            /* Reduced parameter vector             */
for(i=0;i<reduced;i++)                   /* theta_red = theta(theta_index)       */
  cvput(theta_red,i,cvget(theta,(int)cvget(theta_index,i)));
theta_red_new = mmake(reduced,1);        /* A priori update of parameters        */
W1_new    = mmake(hidden,inputs+1);      /* A priori updated W1                  */
W2_new    = mmake(outputs,hidden+1);     /* A priori updated W2                  */
PSI       = mmake(parameters,Nout);      /* Der. of each output wrt. each weight */
minit(PSI);
G         = mmake(reduced,1);            /* Gradient vector                      */
H         = mmake(reduced,reduced);      /* Hessian matrix                       */
Htmp      = mmake(reduced,reduced);      /* Matrix used by the linear sys solver */
h         = mmake(reduced,1);            /* Update vector                        */
all       = mmake(N,1);                  /* Index vector (0:N-1)                 */
for(k=0;k<N;k++) cvput(all,k,(double)k);
index0    = mmake(1,1);                  /* Index vector (0)                     */
put_val(index0,0,0,0);
index7    = mmake(parameters,1);         /* Index vector (0:parameters-1)        */
for(k=0;k<parameters;k++) cvput(index7,k,(double)k);
L_hidden = neuvector(NetDef,1,'L');      /* Location of linear hidden units      */
H_hidden = neuvector(NetDef,1,'H');      /* Location of tanh hidden units        */
L_output = neuvector(NetDef,2,'L');      /* Location of linear output units      */
H_output = neuvector(NetDef,2,'H');      /* Location of tanh output units        */
lhids     = getrows(L_hidden);           /* # of linear hidden units             */
hhids     = getrows(H_hidden);           /* # of tanh hidden units               */
louts     = getrows(L_output);           /* # of linear output units             */ 
houts     = getrows(H_output);           /* # of tanh output units               */
if(hhids>0) tmp0 = mmake(hhids,N);       /* Used to construct PSI                */
else tmp0 = mmake(1,1);
tmp2      = mmake(1,N);                  /* Used to construct PSI                */
tmp3      = mmake(1,N);                  /* Used to construct PSI                */
index2    = mmake(N,1);                  /* Index vector (0:N-1)*outputs         */
for(k=0;k<N;k++) cvput(index2,k,(double)k*outputs);
index     = mmake(hidden,1);             /* Index vector outputs*(hidden+1)+...  */
for(k=0;k<hidden;k++) cvput(index,k,(double)(outputs*(hidden+1)+k*(inputs+1)));
lambda_old = 0.0;
lambda    = trparms->lambda;             /* Levenberg-Marquardt parameter        */
D         = mmake(reduced,1);            /* Initialize vector cont. weight decays*/
Dtmp      = mmake(parameters,1);
if(length(trparms->D)==1)                 /* Scalar weight decay parameters       */
  for(i=0;i<reduced;i++) cvput(D,i,rvget(trparms->D,0));
else if(length(trparms->D)==2)            /* Two weight decay parameters          */
{
  for(i=0;i<parameters2;i++) cvput(Dtmp,i,rvget(trparms->D,0));
  for(i=parameters2;i<parameters;i++) cvput(Dtmp,i,rvget(trparms->D,1));
  for(i=0;i<reduced;i++) cvput(D,i,cvget(Dtmp,(int)cvget(theta_index,i)));
}
else if(length(trparms->D)==reduced){     /* Individual weight decays             */
  for(i=0;i<reduced;i++) cvput(D,i,rvget(trparms->D,i));   
}
else
  mexErrMsgTxt("tparms.D has the wrong length.");
critdif  = trparms->critterm+1.0;        /* Initialize stopping variables        */
gradmax  = trparms->gradterm+1;
paramdif = trparms->paramterm+1;
PI_vector = mmake(trparms->maxiter,1);   /* Vector containing the PI's           */


/* 
-----------------------------------------------------------------------------------
---------------                    TRAIN NETWORK                      -------------
-----------------------------------------------------------------------------------
*/

/* Clear screen on HP systems.
Uncomment the following line and comment the subsequent one */
/*printf("\x1BH\x1BJNetwork training started at %.8s\n\n",asctime(c)+11);*/

printf("\nNetwork training started.\n\n");

/*
 >>>>>>>>>>>>>>       Compute network output y2(theta)          <<<<<<<<<<<<<<< 
*/
    mmul(h1,W1,PHI);
    mtanh(y1,H_hidden,all,h1,H_hidden,all);
    mcopyi(y1,L_hidden,all,h1,L_hidden,all);
    mmul(h2,W2,y1);
    mtanh(y2,H_output,all,h2,H_output,all);
    mcopyi(y2,L_output,all,h2,L_output,all);
    msub(E,Y,y2);                         /* Training error                      */
    m2vreshape2(E_vector,0,E);            /* reshape E intor a long vector       */
    SSE=sprod3(E_vector,E_vector);        /* Sum of squared errors               */
    tmp1 = 0;
    for(i=0;i<reduced;i++)
      tmp1+=cvget(theta_red,i)*cvget(theta_red,i)*cvget(D,i);
    PI = (SSE+tmp1)/(2*N);

/* Iterate until stopping criterion is satisfied */
while (iteration<=trparms->maxiter && PI>trparms->critmin && lambda<1e7 && 
       (critdif>trparms->critterm || gradmax>trparms->gradterm || 
       paramdif>trparms->paramterm))
{
  if(dw==1)
  {
/*
 >>>>>>>>>>>>>>>>>>>>>>>>>>>   COMPUTE THE PSI MATRIX   <<<<<<<<<<<<<<<<<<<<<<<<<<
 (The derivative of each network output (y2) with respect to each weight)
*/
/* Some intermidiate computations */
    for(j=0;j<hhids;j++)
    {
      jj = (int)cvget(H_hidden,j);
      for(k=0;k<N;k++)
	   put_val(tmp0,j,k,1-get_val(y1,jj,k)*get_val(y1,jj,k));
    }

/*   ==========   Elements corresponding to the linear output units   ============*/
    for(i=0; i<louts; i++)
    {
      ii = (int)cvget(L_output,i);

      /***  The part of PSI corresponding to hidden-to-output layer weights ***/
      index1 = ii * (hidden+1);
      psi1(PSI, index1, index2, ii, y1);
      /************************************************************************/

      /**** The part of PSI corresponding to input-to-hidden layer weights ****/
      for(j=0; j<lhids; j++)
      {
	     jj = (int)cvget(L_hidden,j);
        psi2(PSI, (int)cvget(index,jj), index2, ii, get_val(W2,ii,jj), PHI);
      }

      for(j=0; j<hhids;j++)
      {
        jj = (int)cvget(H_hidden,j);
	     psi3(tmp3, tmp0, j, get_val(W2,ii,jj));
	     psi4(PSI, (int)cvget(index,jj), index2, ii, tmp3, PHI);
      }
      /************************************************************************/    
    }


    /* ============  Elements corresponding to the tanh output units   =============*/
    for(i=0; i<houts; i++)
    {
      ii = (int)cvget(H_output,i);
      index1 = ii * (hidden + 1);
      for(k=0; k<N; k++)
	     put_val(tmp2,0,k,1-get_val(y2,ii,k)*get_val(y2,ii,k));

      /* -- The part of PSI corresponding to hidden-to-output layer weights --*/
      psi4(PSI, index1, index2, ii, tmp2, y1);
      /* ---------------------------------------------------------------------*/
    
      /* -- The part of PSI corresponding to input-to-hidden layer weights ---*/
      for(j=0; j<lhids; j++)
      {
        jj = (int)cvget(L_hidden,j);
	     smul(tmp3, tmp2, get_val(W2,ii,jj));
	     psi4(PSI, (int)cvget(index,jj), index2, ii, tmp3, PHI);
      }
      
      for(j=0; j<hhids; j++)
      {
      	jj = (int)cvget(H_hidden,j);
	     psi3(tmp3, tmp0, j, get_val(W2,ii,jj));
        psi5(PSI, (int)cvget(index,jj), index2, ii, tmp3, tmp2, PHI);
      }
      /* ---------------------------------------------------------------------*/
    }
      dw = 0;
      /* -- Gradient (G = PSI_red*E_vector - D*theta_red) -- */
     	for(i=0; i<reduced; i++){
      		ii = (int)cvget(theta_index,i);
    		for(sum=0.0,k=0; k<Nout; k++) sum += get_val(PSI,ii,k)*cvget(E_vector,k);
    		cvput(G,i,sum - cvget(D,i)*cvget(theta_red,i));
      }
    	
    	/* -- Mean square error part of Hessian (PSI_red*PSI_red'  --*/
    	for(i=0; i<reduced; i++){
    		ii = (int)cvget(theta_index,i);
    		for(j=i; j<reduced; j++){
      			jj = (int)cvget(theta_index,j);
      			for(sum=0.0,k=0; k<Nout; k++)
			         sum += get_val(PSI,ii,k)*get_val(PSI,jj,k);
      			put_val(H,i,j,sum);
      			put_val(H,j,i,sum);
    		}
  	  }
     for(i=0;i<reduced;i++)                            /* Add diagonal matrix     */
       put_val(H,i,i,get_val(H,i,i)+cvget(D,i));               
  }

/*
 >>>>>>>>>>>>>>>>>>>>>>>>>>>        COMPUTE h_k        <<<<<<<<<<<<<<<<<<<<<<<<<<<
 */
 
  /* -- Hessian (H = R + lambda*I + D)  --*/
  tmp1 = lambda-lambda_old;
  for(i=0;i<reduced;i++)                            /* Add diagonal matrix     */
    put_val(H,i,i,get_val(H,i,i)+tmp1);               

  /* -- Search direction -- */
  choldc(H, Htmp);
  cholsl(Htmp,h,G);

  /* -- Compute 'apriori' iterate -- */
  madd(theta_red_new,theta_red,h);                  /* Update parameter vector */
  mcopyi(theta,theta_index,index0,theta_red_new,index7,index0);

  /* -- Put the parameters back into the weight matrices -- */
  v2mreshape(W1_new,theta,parameters2);
  v2mreshape(W2_new,theta,0);


/*
 >>>>>>>>>>>>>>       Compute network output y2(theta+h)          <<<<<<<<<<<<<<< 
*/
    mmul(h1,W1_new,PHI);
    mtanh(y1,H_hidden,all,h1,H_hidden,all);
    mcopyi(y1,L_hidden,all,h1,L_hidden,all);
    mmul(h2,W2_new,y1);
    mtanh(y2,H_output,all,h2,H_output,all);
    mcopyi(y2,L_output,all,h2,L_output,all);
    msub(E_new,Y,y2);                     /* Training error                      */
    m2vreshape2(E_new_vector,0,E_new);    /* reshape E intor a long vector       */
    SSE_new = sprod3(E_new_vector,E_new_vector); /* Sum of squared errors        */
    for(tmp1=0.0,i=0; i<reduced; i++)
    {
      tmp1+=cvget(theta_red_new,i)*cvget(theta_red_new,i)*cvget(D,i);
      PI_new = (SSE_new+tmp1)/(2*N);
    }

/*
 >>>>>>>>>>>>>>>>>>>>>>>>>>>       UPDATE  lambda     <<<<<<<<<<<<<<<<<<<<<<<<<<<<
 */
    tmp1 = 0;
    lambda_old = lambda;
    for(i=0;i<reduced;i++) tmp1+=cvget(h,i)*cvget(h,i)*(cvget(D,i)+lambda);
    L = sprod3(h,G) + tmp1;

    /* Decrease lambda if SSE has fallen 'sufficiently' */
    if(2*N*(PI - PI_new) > (0.75*L)) lambda = lambda/2;
  
    /* Increase lambda if SSE has grown 'sufficiently'  */
    else if(2*N*(PI-PI_new) <= (0.25*L)) lambda = 2*lambda;  


/*
 >>>>>>>>>>>>>>>>>>>>       UPDATES FOR NEXT ITERATION        <<<<<<<<<<<<<<<<<<<<<
 */
    /* Update only if criterion has decreased */
    if(PI_new<PI)
    {
     critdif  = PI-PI_new;                               /* Criterion difference */
     for(i=0,gradmax=0.0,ptm1=G->mat[0];i<reduced;i++){  /* Maximum gradient     */
        sum = fabs(*(ptm1++));
        if(gradmax<sum)
           gradmax = sum;
     }
     gradmax/=N;
     ptm1=theta_red_new->mat[0];
     ptm2=theta_red->mat[0];
     for(i=0,paramdif=0.0;i<reduced;i++){  /* Maximum gradient     */
        sum = fabs(*(ptm1++) - *(ptm2++));
        if(paramdif<sum)
           paramdif = sum;
     }
     lambda_old = 0.0;
     tmp = W1; W1=W1_new; W1_new=tmp;
     tmp = W2; W2=W2_new; W2_new=tmp;
     tmp = theta_red; theta_red=theta_red_new; theta_red_new = tmp;
     tmp = E_vector; E_vector = E_new_vector; E_new_vector = tmp;
     dw = 1;
     PI = PI_new;
     cvput(PI_vector,iteration-1,PI);
     switch(trparms->infolevel){                            /* Print on-line inform */
       case 1:
          printf("# %i   W=%4.3e  critdif=%3.2e  maxgrad=%3.2e  paramdif=%3.2e\n",
                                                  iteration,PI,critdif,gradmax,paramdif);
       break;
       default:
          printf("iteration # %i   W = %4.3e\r",iteration,PI);
     }
     ++iteration;
   }
}


/*
 >>>>>>>>>>    RETURN POINTERS TO RETURN ARGUMENTS & FREE MEMORY    <<<<<<<<<<<<
 */
 
/* Swap pointers if they have been messed up */
if ((iteration&1) == 0) {
	mset(W1_new,W1);
	tmp = W1; W1=W1_new; W1_new=tmp;
	mset(W2_new,W2);
     	tmp = W2; W2=W2_new; W2_new=tmp;
}
iteration=iteration-1;
if(iteration==0){
	*PI_vectorpp = mmake(1,1);
	(*PI_vectorpp)->row=0;
	(*PI_vectorpp)->col=0;
}
else {
	*PI_vectorpp = mmake(iteration,1);
	subvec(*PI_vectorpp,PI_vector,0,iteration-1);
}
*iter = iteration;
*lam  = lambda;
mfree(L_hidden); mfree(H_hidden); mfree(L_output); mfree(H_output);
mfree(h1); mfree(h2); mfree(y1); mfree(y2);
mfree(E); mfree(E_new); mfree(E_vector); mfree(E_new_vector);
mfree(W1_new); mfree(W2_new); mfree(D);mfree(Dtmp); mfree(PI_vector);mfree(Htmp);
mfree(theta); mfree(thtmp); mfree(theta_index); mfree(theta_red); mfree(theta_red_new);
mfree(PSI); mfree(G); mfree(H); mfree(h);
mfree(all); mfree(index0); mfree(index7); mfree(onesvec); mfree(tmp0); mfree(tmp2);
mfree(tmp3); mfree(index); mfree(index2); mfree(PHI);
printf("\n\nNetwork training ended.\n\n\n");
}
/*
  --------------------------------------------------------------------------------
  ----------------             END OF NETWORK TRAINING              --------------
  --------------------------------------------------------------------------------
*/

/********************************************************************************
 *                                                                              *
 *                      SOME ADDITIONAL MATRIX FUNCTIONS                        *
 *                      --------------------------------                        *
 *                                                                              *
 ********************************************************************************/

/*
 * MAT2SM
 * ------
 *
 * Conversion of a (real) Matlab matrix into the SteffenMagnus format
 *
 * INPUT : MatlabMatrix - Matrix in Matlab form
 * OUTPUT: SMmatrix     - Matrix in SM form
 *
 */
matrix *mat2sm(const mxArray *MatlabMatrix)
{
  int rows, cols, i, j;
  double *M;
  matrix *SMmatrix;

  rows = mxGetM(MatlabMatrix);
  cols = mxGetN(MatlabMatrix);
  M    = mxGetPr(MatlabMatrix);

  SMmatrix = mmake(rows,cols);
  for(i=0;i<cols;i++)
  {
    for(j=0;j<rows;j++)
    {
      put_val(SMmatrix,j,i,M[j+i*rows]);
    }
  }
  return SMmatrix;
}


matrix *matstring2sm(const mxArray *MatlabMatrix)
{
  int rows, cols, i, j;
  unsigned short *M;
  matrix *SMmatrix;

  rows = mxGetM(MatlabMatrix);
  cols = mxGetN(MatlabMatrix);
  M    = (unsigned short*)mxGetPr(MatlabMatrix);

  SMmatrix = mmake(rows,cols);
  for(i=0;i<cols;i++)
  {
    for(j=0;j<rows;j++)
    {
      put_val(SMmatrix,j,i,(double)M[j+i*rows]);
    }
  }
  return SMmatrix;
}

/*
 * SM2MAT
 * ------
 *
 * Conversion of a SteffenMagnus matrix to a (real) Matlab matrix format
 *
 * INPUTS: MatlabMatrix - Pointer to Matlab matrix 
 *         SMmatrix     - Pointer to SM matrix
 *
 */
void sm2mat(mxArray *MatlabMatrix, matrix *SMmatrix)
{
  int rows, cols, i, j;
  double *M;
  rows = mxGetM(MatlabMatrix);
  cols = mxGetN(MatlabMatrix);
  M    = mxGetPr(MatlabMatrix);

  for(i=0;i<cols;i++)
  {
    for(j=0;j<rows;j++)
    {
      M[j+i*rows]=get_val(SMmatrix,j,i);
    }
  }
  mfree(SMmatrix);
}


/* 
 * MCOPYI
 * ------
 *         Copy the elements in matrix B pointed out by the index vectors
 *         bri (rows) and bci (columns) to matrix A (to the locations
 *         specified by bri and bci)
 * 
 */
void mcopyi(matrix *A, matrix *ari, matrix *aci, matrix *B, matrix *bri,\
            matrix *bci)
{
  register int i,j, rows, cols;
  rows = ari->row;
  cols = aci->row;
  for(i=0;i<rows;i++)
  {
    for(j=0;j<cols;j++)
    {
      A->mat[(int)cvget(ari,i)][(int)cvget(aci,j)] = \
           B->mat[(int)cvget(bri,i)][(int)cvget(bci,j)];
    }
  }
}


/*
 * MTANH
 * -----
 *        Computes tanh to the elements in matrix B pointed out by the
 *        index vectors bri (rows) and bci (columns), and place the result
 *        in matrix A (in the locations specified by ari and aci).
 *
 * INPUTS: A, B - Pointers to matrices
 *         ari, aci, bri, bci - Index vectors
 *
 */
void mtanh(matrix *A, matrix *ari, matrix *aci, matrix *B, matrix *bri,\
              matrix *bci)
{
  register int i,j, rows, cols;
  rows = ari->row;
  cols = aci->row;
  for(i=0;i<rows;i++)
  {
    for(j=0;j<cols;j++)
    {
      A->mat[(int)cvget(ari,i)][(int)cvget(aci,j)] = \
           tanh(B->mat[(int)cvget(bri,i)][(int)cvget(bci,j)]);
    }
  }
}


/*
 * NEUVECTOR:
 * ----------
 *              This function finds the locations of a given type
 *              of neurons in a row of the NetDef matrix.
 *
 *  INPUTS: NetDef - Network definition string matrix
 *          layer  - Layer of interest in NetDef
 *          neu    - Neuron type to search for
 *
 */
matrix *neuvector(matrix *NetDef,int layer, char neu)
{
  matrix *tmp, *retvec, *NDef;
  int i;
  NDef = mmake(1,getcols(NetDef));
  submat(NDef,NetDef,layer-1,layer-1,0,getcols(NDef)-1);
  tmp = mfind(NDef,(double)neu);
  if(getrows(tmp)==0)
  {
    retvec = mmake(1,1);
    retvec->row=0;
    retvec->col=0;
  }
  else
  {
    retvec = mmake(getrows(tmp),1);
    for(i=0;i<getrows(tmp);i++) cvput(retvec,i,mget(tmp,i,1));
  }
  return retvec;
}



/*
 * SSE
 * ---
 *      Given an error matrix, the function finds the sum of
 *      squared error.
 *
 *      INPUT:  E   - Error matrix [ny | N]
 *      OUTPUT: SSE - Sum of squared error
 */
double sse(matrix *E)
{
  double SSE;
  register int i,j;
  SSE=0.0;
  for(j=0;j<E->col;j++)
  {
    for(i=0;i<E->row;i++)
    {
      SSE+=E->mat[i][j] * E->mat[i][j];
    }
  }
  return SSE;
}


/*
 * DELTATANH
 * ---------
 *
 * Calculates the deltas for tanh units
 *
 * INPUTS: delta   - The return argument
 *         y       - The output of the units
 *         E       - The back propagated error
 *         H_vec   - Indices to the layers tanh units
 */ 
void deltatanh(matrix *delta, matrix *y, matrix *E, matrix *H_vec)
{
	register int i, j, N, veclen;
	N = getcols(y);
        veclen=getrows(H_vec);
	for(i=0;i<veclen;i++)
	{
	  for(j=0;j<N;j++)
          {
	    put_val(delta, (int)cvget(H_vec,i), j,\
              (1 - get_val(y,(int)cvget(H_vec,i),j) * get_val(y,(int)cvget(H_vec,i),j))\
				  *get_val(E,(int)cvget(H_vec,i),j));
	  }
        }
}

/*
 * MMULTR1
 * -------
 *
 * Matrix multiplication ptm1 = ptm2'*ptm3
 *
 * Input:  *ptm1 - Pointer to result matrix (Not equal to *ptm1 or *ptm2) 
 *         *ptm2 - Pointer to first argument matrix (the one to be
 *                 transposed)
 *         *ptm3 - Pointer to second argument matrix
 *
 */
void mmultr1( matrix *ptm1, matrix *ptm2, matrix *ptm3 )
{
	register int i,j,k;

#if RUNCHK

	if ((ptm1==ptm2) || (ptm1==ptm3)) \
	merror("Memory conflict error in mmultr1!");

	if ( !( (ptm2->row == ptm3->row) && \
	(ptm2->col == ptm1->row) && (ptm3->col == ptm1->col) ) ) \
	merror("Dimension mismatch error in mmultr1!");

#endif

	for ( i=0; i < ptm2->col; i++ )
	{
		for ( j=0; j < ptm3->col; j++ )
		{
			ptm1->mat[i][j] = 0.0;
			for ( k=0; k < ptm2->col; k++ )
			{
				ptm1->mat[i][j] += ptm2->mat[k][i] * ptm3->mat[k][j];
			}
		}
	}
}



/*
 * MMULTR2
 * -------
 *
 * Matrix multiplication ptm1 = ptm2*ptm3'
 *
 * Input:  *ptm1 - Pointer to result matrix (Not equal to *ptm1 or *ptm2) 
 *         *ptm2 - Pointer to first argument matrix
 *         *ptm3 - Pointer to second argument matrix (the one to
 *                 be transposed)
 *
 */
void mmultr2( matrix *ptm1, matrix *ptm2, matrix *ptm3 )
{
	register int i,j,k;

#if RUNCHK

	if ((ptm1==ptm2) || (ptm1==ptm3)) \
	merror("Memory conflict error in mmultr1!");

	if ( !( (ptm2->col == ptm3->col) && \
	(ptm2->row == ptm1->row) && (ptm3->row == ptm1->col) ) ) \
	merror("Dimension mismatch error in mmultr2!");

#endif

	for ( i=0; i < ptm2->row; i++ )
	{
		for ( j=0; j < ptm3->row; j++ )
		{
			ptm1->mat[i][j] = 0.0;
			for ( k=0; k < ptm2->col; k++ )
			{
				ptm1->mat[i][j] += ptm2->mat[i][k] * ptm3->mat[j][k];
			}
		}
	}
}



/* SPROD3 :
 * --------
 * This function calculates the scalar-product of two COLUMN vectors of equal
 * length.
 * Inputs : ptv1 - Pointer to first factor vector.
 *          ptv2 - Pointer to second factor vector.
 * Output : prod - Scalar product of the two vectors.
 */
double sprod3( matrix* ptv1, matrix* ptv2 )
{
	register double prod = 0.0;
	register int i, elements;

#if RUNCHK

	/* Check whether vectors has the same length. */
	if (ptv1->row != ptv2->row) 
	merror("Dimension mismatch error in sprod!");

#endif
	/* Calculate scalar product of the vectors.   */

	elements = ptv1->row;

	for ( i=0; i<elements; i++ ) {
		prod += (cvget(ptv1,i))*(cvget(ptv2,i));
	}

	return prod;
}




/*
 * MAKEINDEXVEC
 * ------------
 *               Create an index vector.
 *
 * INPUTS: start - First element
 *         end   - Last element
 *         step  - Step
 * OUTPUT: indexvector
 *
 */
matrix *makeindexvec(int start, int end, int step)
{
/* FUNCTION DOES NOT WORK AND IS NOT USED */
  int k, elements;
  matrix *indexvec;
  elements = (end-start)/step;
  indexvec=mmake(elements,1);
  for(k=0;k<elements;k++)
  {
    cvput(indexvec,k,(double)start);
    start+=step;
  }
  return indexvec;
}



/*
 * M2VRESHAPE
 * ----------
 *
 *         Matrix to Vector Reshape. Take the elements from the matrix M ROW wise
 *         and put them into the vector V starting from index1
 *
 * INPUTS: V      - Pointer to vector
 *         index1 - Start index in V
 *         M      - Pointer to matrix
 *
 */
void m2vreshape(matrix *V, int index1, matrix *M)
{
  int i, j;
  for(i=0;i<getrows(M);i++)
  {
    for(j=0;j<getcols(M);j++)
    {
      cvput(V,index1++,get_val(M,i,j));
    }
  }
}



/*
 * M2VRESHAPE2
 * -----------
 *
 *         Matrix to Vector Reshape. Take the elements from the matrix M COLUMN wise
 *         and put them into the vector V starting from index1
 *
 * INPUTS: V      - Pointer to vector
 *         index1 - Start index in V
 *         M      - Pointer to matrix
 *
 */
void m2vreshape2(matrix *V, int index1, matrix *M)
{
  int i, j;
  for(i=0;i<nof_cols(M);i++)
  {
    for(j=0;j<nof_rows(M);j++)
    {
      cvput(V,index1++,get_val(M,j,i));
    }
  }
}


/*
 * V2MRESHAPE
 * ----------
 *
 *         Vector to Matrix Reshape. Take the elements from the vector V, starting
 *         at index1, and insert them into the matrix M, ROW wise.
 *
 * INPUTS: M      - Pointer to matrix
 *         V      - Pointer to vector
 *         index1 - Start index in V
 *
 */
void v2mreshape(matrix *M, matrix *V, int index1)
{
  int i, j;
  for(i=0;i<nof_rows(M);i++)
  {
    for(j=0;j<nof_cols(M);j++)
    {
      put_val(M,i,j,cvget(V,index1++));
    }
  }
}



/*
 * CHOLDC
 * ------
 *         This routine constructs the Cholesky decomposition, A = L*L' of a positive
 *         definite symmetric matrix.
 *
 * INPUTS: A - Matrix to be decomposed
 *         L - Matrix of same dimensions used for storage
 *
 */
void choldc(matrix *A, matrix *L)
{
  int i, j, k, n;
  double sum;
  n = nof_rows(A);
  for(i=0; i<n; i++)
  {
    for(j=0; j<n; j++)
      {
	for(sum=get_val(A,i,j), k=i-1; k>=0; k--) sum-=get_val(L,i,k)*get_val(L,j,k);
	if(i==j)
	{
	  if( sum<= 0.0) merror("Choldc failed. Matrix not positive definite");
	  put_val(L,i,i,sqrt(sum));
	}
	else put_val(L,j,i,sum/get_val(L,i,i));
      }
  }
}

	

/*
 * CHOLSL
 * ------
 *         This routine constructs the Cholesky decomposition, A = L*L' of a positive
 *         definite symmetric matrix.
 *
 * INPUTS: A - Matrix to be decomposed
 *         L - Matrix of same dimensions used for storage
 *
 */
void cholsl(matrix *L, matrix *x, matrix *b)
{
  int i, k, n;
  double sum;
  n = nof_rows(b);
  for(i=0; i<n; i++)
  {
    for(sum=cvget(b,i), k=i-1; k>=0; k--) sum -= get_val(L,i,k)*cvget(x,k);
    cvput(x,i,sum/get_val(L,i,i));
  }
  for(i=n-1; i>=0; i--)
  {
    for(sum=cvget(x,i), k=i+1; k<n; k++) sum -= get_val(L,k,i)*cvget(x,k);
    cvput(x,i,sum/get_val(L,i,i));
  }
}


/*
 * PSI1
 * ----
 *
 */
void psi1(matrix *PSI, int index1, matrix *index2, int i, matrix *y1)
{
  int l,k, idx, hidden, N;
  hidden = nof_rows(y1);
  N      = nof_cols(y1);
  for(l=0; l<hidden; l++)
  {
    idx = index1+l;
    for(k=0; k<N; k++)
      put_val(PSI,idx,(int)cvget(index2,k)+i,get_val(y1,l,k));
  }
}


/*
 * PSI2
 * ----
 *
 */
void psi2(matrix *PSI, int index, matrix *index2, int i, double w, matrix *PHI)
{
  int l, k, inputs, N, idx;
  inputs = nof_rows(PHI);
  N      = nof_cols(PHI);
  for(l=0; l<inputs; l++)
  {
    idx = index+l;
    for(k=0; k<N; k++)
       put_val(PSI,idx,(int)cvget(index2,k)+i,w*get_val(PHI,l,k));
  }
}


/*
 * PSI3
 * ----
 *
 */
void psi3(matrix *tmp3, matrix *tmp0, int j, double w)
{
  int k, N;
  N = nof_cols(tmp0);
  for(k=0; k<N; k++) rvput(tmp3,k,get_val(tmp0,j,k)*w);
}


/*
 * PSI4
 * ----
 */
void psi4(matrix *PSI, int index, matrix *index2, int i, matrix *tmp, matrix *PHI)
{
  int l, k, inputs, N, idx;
  inputs = nof_rows(PHI);
  N      = nof_cols(PHI);
  for(l=0; l<inputs; l++)
  {
    idx = index+l;
    for(k=0; k<N; k++)
       put_val(PSI,idx,(int)cvget(index2,k)+i,rvget(tmp,k)*get_val(PHI,l,k));
  }
}


/* 
 * PSI5
 * ----
 */
void psi5(matrix *PSI, int index, matrix *index2, int i, matrix *tmp3, matrix *tmp2,\
                                                                        matrix *PHI)
{
  int l, k, inputs, N, idx;
  inputs = nof_rows(PHI);
  N      = nof_cols(PHI);
  for(l=0; l<inputs; l++)
  {
    idx = index+l;
    for(k=0; k<N; k++)
       put_val(PSI,idx,(int)cvget(index2,k)+i,rvget(tmp3,k)*rvget(tmp2,k)\
                                                              *get_val(PHI,l,k));
  }
}


/* 
 * VCOPYI
 * ------
 *         vcopyi(A,ari,ac,B,bri,bc)
 *
 *         From matrix B copy the column 'bc' and the rows specified by 'bri' to
 *         column 'ac' and the rows specified by 'ari' in matrix A.
 *
 * ARGUMENTS: A, B     - Pointers to matrices (matrix*)
 *            ari, bri - Row index vectors (matrix*)
 *            aci, bci - Column indices (int)
 */
void vcopyi(matrix *A, matrix *ari, int ac, matrix *B, matrix *bri,int bc)
{
  register int i,j, rows;
  rows = ari->row;
  for(i=0;i<rows;i++)
  {
      A->mat[(int)cvget(ari,i)][ac] = B->mat[(int)cvget(bri,i)][bc];
  }
}


/*
 * VTANH
 * -----
 *        vtanh(A, ari, ac, B, bri, bc)
 *
 *        Computes tanh to the entries in matrix B pointed out by the
 *        row index vector 'bri' and the column index 'bci', and put the
 *        result in matrix A (at the locations specified by ari and aci).
 *
 * ARGUMENTS: A, B     - Pointers to matrices (matrix*)
 *            ari, bri - Row index vectors (matrix*)
 *            aci, bci - Column indices (int)
 *
 */
void vtanh(matrix *A, matrix *ari, int ac, matrix *B, matrix *bri, int bc)
{
  register int i,j, rows;
  rows = ari->row;
  for(i=0;i<rows;i++)
  {
      A->mat[(int)cvget(ari,i)][ac] = tanh(B->mat[(int)cvget(bri,i)][bc]);
  }
}


/*
 * MVMUL
 * -----
 *         mmul(ptm1, ptm2, ptm3, col3)
 *
 * Matrix-vector multiplication: ptm1 = ptm2*ptm3(:,col3)
 *
 * Arguments:  *ptm1 - Pointer to column vector (Not equal to *ptm1 or *ptm2)
 *             *ptm2 - Pointer to left matrix
 *             *ptm3 - Pointer to matrix from which the vector is extracted.
 *             col3  - Index to column of matrix ptm multiplied with ptm2
 *
 */
void mvmul( matrix *ptm1, matrix *ptm2, matrix *ptm3, int col3 )
{
	register int i,j,k;

	for ( i=0; i < ptm2->row; i++ ){
		ptm1->mat[i][0] = 0.0;
		for ( k=0; k < ptm2->col; k++ ){
			ptm1->mat[i][0] += ptm2->mat[i][k] * ptm3->mat[k][col3];
		}
	}
}

