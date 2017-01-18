/*
 *     INCLUDE HEADERS
 */
#include <stdio.h>
#include <math.h>
#include "mex.h"
#include "matrix2.h"
#include "nnmisc.h"

void nnssif(matrix**, int*, double*, matrix*, int, matrix*, matrix*, matrix*,\
 trparmstruct*, matrix*, matrix*);



/*********************************************************************************
 *                                                                               *
 *    NNSSIF                                                                     *
 *    ------                                                                     *
 *                                                                               *
 *    This is a CMEX-version of the Matlab function nnssif.                      *
 *    Type 'help nnssif' from Matlab for information on                          *
 *    how to call this function.                                                 *
 *                                                                               *
 *                                                                               *
 *    Programmed by: Magnus Norgaard                                             *
 *    LastEditDate : Jan. 20, 2000                                               *
 *                                                                               *
 *********************************************************************************/
void nnssif(matrix **NSSEvecpp, int *iter, double *lam,\
	matrix *NetDef, int nx, matrix *W1, matrix *W2, matrix *obsidx,\
	trparmstruct *trparms, matrix *Y, matrix *U)
{
/*
----------------------------------------------------------------------------------- 
---------------              VARIABLE DECLARATIONS                    ------------- 
----------------------------------------------------------------------------------- 
*/ 
register i, j, k, t; 
int outputs, N, Nout, layers, dummy, hidden, inputs, iteration; 
int parameters1, parameters2, parameters, reduced, index1, ii, jj; 
int lhids, hhids, louts, houts, hid1, hid2, Nny, nxu, skipstart;
int Ndat, N2, nu, ny, index5, index6, dummy2, skip;
double lambda, SSE, SSE_new, NSSE, NSSE_new, L, tmp1, sum, dummy3; 
double critdif, gradmax, paramdif, *ptm1, *ptm2, lambda_old;
char dw, stateflag; 
matrix *L_hidden, *H_hidden, *L_output, *H_output, *h1, *h2, *y1, *y2; 
matrix *E, *Evec, *Evec_new, *W1_new, *W2_new, *PHI, *D, *Dtmp; 
matrix *NSSEvec, *miter, *tmp, *Htmp, *PSIx, *Yhat; 
matrix *theta, *thtmp, *theta_index, *theta_red, *theta_red_new, *PSI, *G, *H, *h; 
matrix *all, *index0, *index7, *onesvec, *tmp0, *tmp2, *tmp3, *index, *index2;
matrix *rowidx, *nrowidx, *dxdy1, *Ahat, *Khat, *dy1de, *dy1dx, *Y2, *dummy1, *C;
matrix *AKC;
 
 
/* 
----------------------------------------------------------------------------------- 
---------------             NETWORK INITIALIZATIONS                   ------------- 
----------------------------------------------------------------------------------- 
 */
Ndat      = getcols(Y);                  /* # of data                            */
ny        = getrows(Y);                  /* # of outputs                         */
nu        = getrows(U);                  /* # of controls                        */
N         = Ndat - 1;                    /* Size of training set                 */
Nny       = N*ny;
nxu       = nx+nu;
skip      = trparms->skip;
N2        = N-skip;
skipstart = ny*skip;
Y2        = mmake(ny,N);                 /* Observed outputs used for training   */
Yhat      = mmake(ny,1);                 /* Output prediction                    */
hidden    = getcols(NetDef);             /* # of hidden units                    */
inputs    = nx+nu+ny;                 	 /* Number of inputs to network          */
outputs   = nx;		                 /* Always one outputs                   */
L_hidden  = neuvector(NetDef,1,'L');     /* Location of linear hidden units      */
H_hidden  = neuvector(NetDef,1,'H');     /* Location of tanh hidden units        */ 
L_output  = neuvector(NetDef,2,'L');     /* Location of linear output units      */ 
H_output  = neuvector(NetDef,2,'H');     /* Location of tanh output units        */ 
lhids     = getrows(L_hidden);           /* # of linear hidden units             */ 
hhids     = getrows(H_hidden);           /* # of tanh hidden units               */ 
louts     = getrows(L_output);           /* # of linear output units             */  
houts     = getrows(H_output);           /* # of tanh output units               */
if (louts+houts!=nx) mexErrMsgTxt("Output layer of NetDef does not match specified nx.");
miter     = mmake(1,1);                  /* Temp element                         */ 
h1        = mmake(hidden,1);             /* Argument to hidden layer act. fcts   */ 
h2        = mmake(outputs,1);            /* Argument to hidden layer act. fcts   */ 
onesvec   = mmake(1,N);                  /* Vector of all ones                   */
minitx(onesvec,1.0);
y1        = mmake(hidden+1,N);           /* Hidden layer outputs                 */
minit(y1);
mat2mat(y1,hidden,0,onesvec);            /* Add a row of ones (bias to outputs)  */ 
y2        = mmake(outputs,N);            /* Output layer output                  */ 
minit(y2);
E         = mmake(ny,1);                 /* Prediction error                     */
Evec      = mmake(1,Nny);                /* Prediction error vector              */
Evec_new  = mmake(1,Nny);                /* A priori Evec                        */ 
index     = mmake(hidden,1);             /* Index vector outputs*(hidden+1)+...  */
for(k=0;k<hidden;k++) cvput(index,k,(double)(outputs*(hidden+1)+k*(inputs+1))); 
index2    = mmake(N,1);                  /* Index vector (0:N-1)*outputs         */
for(k=0;k<N;k++) cvput(index2,k,(double)k*outputs); 
iteration = 1;                           /* Initialize iteration counter         */
dw        = 1;                           /* Flag telling that the weights are new*/ 
parameters1= hidden*(inputs+1);          /* # of input-to-hidden weights         */
parameters2= outputs*(hidden+1);         /* # of hidden-to-output weights        */ 
parameters = parameters1+parameters2;    /* Total # of weights                   */ 
rowidx    = mmake(ny,1);                 /* Row indices                          */
vput(rowidx,0,vget(obsidx,0));
for(k=1;k<ny;k++) vput(rowidx,k,vget(obsidx,k)+vget(rowidx,k-1));
for(k=0;k<ny;k++) vput(rowidx,k,vget(rowidx,k)-1.0);
if(nx-ny!=0){
	nrowidx = mmake(nx-ny,1);        /* Not row indices                      */
	for(j=0,k=0;k<nx;k++){
		stateflag=0;
		for(i=0;i<ny;i++){
			if((int)vget(rowidx,i)==k) stateflag=1;
		}
		if(stateflag==0) vput(nrowidx,j++,(double)k);
	}
}		
else{
	nrowidx=mmake(1,1);
	nrowidx->row = 0;
	nrowidx->col = 0;
}
/* Initialize weights if necessary */
if(getrows(W2)==0){
	W2->row=outputs;
	mrand(W1); smul(W1,W1,0.025);
   	mrand(W2); smul(W2,W2,0.025);
}

/* Insert zeros to ensure observability */
hid1 = (int)floor(hidden/2+0.5);
hid2 = hidden - hid1;
for(i=hid1;i<hidden;i++){
	for(j=0;j<nx;j++) mput(W1,i,j,0.0);
	for(j=0;j<ny;j++) mput(W2,(int)vget(rowidx,j),i,0.0);
}
for(i=0;i<hid1;i++){
	for(j=0;j<(nx-ny);j++) mput(W2,(int)vget(nrowidx,j),i,0.0);
}
/* Observation matrix */
C          = mmake(ny,nx); minit(C);
mput(C,0,0,1.0);
for(i=1;i<ny;i++) mput(C,i,vget(rowidx,i-1)+1,1.0);
Ahat       = mmake(nx,nx);               /* Deriv. of states wrt. past states    */
Khat       = mmake(nx,ny);               /* Deriv. of states wrt. past residuals */
AKC        = mmake(nx,nx);               /* Stores temp. results                 */
W1_new     = mmake(hidden,inputs+1);     /* A priori updated W1                  */
W2_new     = mmake(outputs,hidden+1);    /* A priori updated W2                  */ 
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
dy1dx     = mmake(hidden,nx);            /* Der. of hid. outp. wrt. past states  */
dy1de     = mmake(hidden,ny);            /* Der.of hid. outp. wrt. past pred. err*/
dxdy1     = mmake(nx,hidden);            /* Der. of state estim. wrt. hid. outp. */
PSI       = mmake(parameters,N*ny);      /* Der. of each output wrt. each weight */
PSIx      = mmake(parameters,N*nx);      /* Der. of est. states wrt. each weight */
minit(PSIx);
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
if (hhids>0) tmp0 = mmake(hhids,N);      /* Used to construct PSI                */
else tmp0 = mmake(1,1);
tmp2      = mmake(1,N);                  /* Used to construct PSI                */ 
tmp3      = mmake(1,N);                  /* Used to construct PSI                */ 
lambda    = trparms->lambda;             /* Levenberg-Marquardt parameter        */
lambda_old = 0.0;
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
NSSEvec  = mmake(trparms->maxiter,1);   /* Vector containing the PI's           */


/*
 >>>>>>>>>>>>>>>>>>>>  CONSTRUCT THE REGRESSION MATRIX PHI   <<<<<<<<<<<<<<<<<<<<<
 */
PHI = mmake(inputs+1,N);      	         /* Matrix of input vectors (incl. bias) */
minit(PHI);
mat2mat(PHI,inputs,0,onesvec);           /* Insert biases in PHI                 */
for(i=0;i<nu;i++){                       /* Insert controls in PHI               */
	for(t=0;t<N;t++)
		mput(PHI,nx+i,t,mget(U,i,t));
}
for(i=0;i<ny;i++){
	for(t=0;t<N;t++) mput(Y2,i,t,mget(Y,i,t+1));
}


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
for(t=0;t<N;t++){
	mvmul(h1,W1,PHI,t);
	vtanh(y1,H_hidden,t,h1,H_hidden,0);
	vcopyi(y1,L_hidden,t,h1,L_hidden,0);

	mvmul(h2,W2,y1,t);
	vtanh(y2,H_output,t,h2,H_output,0);
	vcopyi(y2,L_output,t,h2,L_output,0);

	for(k=0;k<(nx-ny);k++){
		i=(int)cvget(nrowidx,k);
		y2->mat[i][t]+=get_val(PHI,i+1,t);
	}
	mvmul(Yhat,C,y2,t);                             /* Output prediction     */
	
	for(i=0;i<ny;i++){
		cvput(E,i,get_val(Y2,i,t)-cvget(Yhat,i));/* Prediction error     */
		rvput(Evec,t*ny+i,cvget(E,i));          /*Store E in vector Evec*/
	}	
	if(t<N-1){
		for(i=0;i<nx;i++)
			put_val(PHI,i,t+1,get_val(y2,i,t));
		for(i=0;i<ny;i++)
			put_val(PHI,nx+nu+i,t+1,cvget(E,i));
	}
}
for(SSE=0,t=skipstart;t<Nny;t++)
	SSE+=rvget(Evec,t)*rvget(Evec,t);               /* Sum of squared errors */
for(tmp1=0,i=0;i<reduced;i++) tmp1+=cvget(theta_red,i)*cvget(theta_red,i)*cvget(D,i); 
NSSE = (SSE+tmp1)/(2*N2);                               /* Value of cost function*/


/* Iterate until stopping criterion is satisfied */
while (iteration<=trparms->maxiter && NSSE>trparms->critmin && lambda<1e7 && 
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

/*   ==========   Elements corresponding to the linear output units   ===========*/
    for(i=0; i<louts; i++)
    {
      ii = (int)cvget(L_output,i);

      /***  The part of PSIx corresponding to hidden-to-output layer weights ***/
      index1 = ii * (hidden+1);
      psi1(PSIx, index1, index2, ii, y1);
      /************************************************************************/

      /**** The part of PSIx corresponding to input-to-hidden layer weights ****/
      for(j=0; j<lhids; j++)
      {
	jj = (int)cvget(L_hidden,j);
        psi2(PSIx, (int)cvget(index,jj), index2, ii, get_val(W2,ii,jj), PHI);
      }

      for(j=0; j<hhids;j++)
      {
        jj = (int)cvget(H_hidden,j);
	psi3(tmp3, tmp0, j, get_val(W2,ii,jj));
	psi4(PSIx, (int)cvget(index,jj), index2, ii, tmp3, PHI);
      }
      /************************************************************************/    
    }


    /* ===========  Elements corresponding to the tanh output units   ===========*/
    for(i=0; i<houts; i++)
    {
      ii = (int)cvget(H_output,i);
      index1 = ii * (hidden + 1);
      for(k=0; k<N; k++)
	put_val(tmp2,0,k,1-get_val(y2,ii,k)*get_val(y2,ii,k));

      /* -- The part of PSIx corresponding to hidden-to-output layer weights --*/
      psi4(PSIx, index1, index2, ii, tmp2, y1);
      /* ---------------------------------------------------------------------*/
    
      /* -- The part of PSIx corresponding to input-to-hidden layer weights ---*/
      for(j=0; j<lhids; j++)
      {
        jj = (int)cvget(L_hidden,j);
	smul(tmp3, tmp2, get_val(W2,ii,jj));
	psi4(PSIx, (int)cvget(index,jj), index2, ii, tmp3, PHI);
      }
      
      for(j=0; j<hhids; j++)
      {
      	jj = (int)cvget(H_hidden,j);
	psi3(tmp3, tmp0, j, get_val(W2,ii,jj));
        psi5(PSIx, (int)cvget(index,jj), index2, ii, tmp3, tmp2, PHI);
      }
      /* ---------------------------------------------------------------------*/
    }


    /* 
     >>>>>>>>>>>>>>>>>>>>        Linearize network           <<<<<<<<<<<<<<<<<<<<<  
    */
    for(t=0;t<N;t++){
    	/*-- Derivative of states wrt. hidden outputs --*/
    	for(j=0;j<louts;j++){
    		i=(int)cvget(L_output,j);
    		for(k=0;k<hidden;k++) put_val(dxdy1,i,k,get_val(W2,i,k));
    	}
    	for(j=0;j<houts;j++){
    		i=(int)cvget(H_output,j);
    		for(k=0;k<hidden;k++) put_val(dxdy1,i,k,get_val(W2,i,k)*(1-\
    	                             	get_val(y2,i,t)*get_val(y2,i,t)));
    	}

      	/*-- Partial deriv. of output from each hidden unit wrt. net inputs --*/
      	for(j=0;j<lhids;j++){
      		i=(int)cvget(L_hidden,j);
		for(k=0;k<nx;k++) put_val(dy1dx,i,k,get_val(W1,i,k));    		
      		for(k=nxu;k<inputs;k++) put_val(dy1de,i,k-nxu,get_val(W1,i,k));
      	}
      	for(j=0;j<hhids;j++){
      		i=(int)cvget(H_hidden,j);
      		for(k=0;k<nx;k++) put_val(dy1dx,i,k,\
      			get_val(W1,i,k)*(1-get_val(y1,i,t)*get_val(y1,i,t)));
      		for(k=nxu;k<inputs;k++) put_val(dy1de,i,\
      			k-nxu,get_val(W1,i,k)*(1-get_val(y1,i,t)*get_val(y1,i,t)));
      	}

     	/*--Partial derivative of states w.r.t. past states and residuals --*/
     	mmul(Ahat,dxdy1,dy1dx);
     	for(k=0;k<nx-ny;k++) put_val(Ahat,(int)cvget(nrowidx,k),\
     					(int)cvget(nrowidx,k)+1,1.0);                
     	mmul(Khat,dxdy1,dy1de);
     	mmul(AKC,Khat,C);
     	msub(AKC,Ahat,AKC);


    /* 
     >>>>>>>>>>>>>>>>>>>     Filter partial derivatives        <<<<<<<<<<<<<<<<<<<<  
    */
    	if(t>=1){
    		/* PSIx = PSIx + PSIx*AKC' */
    		index5 = t*nx;
    		index6 = (t-1)*nx;
		for(i=0;i<reduced;i++){
			ii =(int)cvget(theta_index,i);
			for(j=0;j<nx;j++){
			   for(k=0;k<nx;k++){
			     PSIx->mat[ii][index5+j]+=get_val(PSIx,ii,index6+k)*\
				 	get_val(AKC,j,k);
			   }
			}
		}
	}

	
	/*PSI=PSIx*C';*/
	index5 = t*ny;
    	index6 = t*nx;
	for(i=0;i<reduced;i++){
		ii =(int)cvget(theta_index,i);
		for(j=0;j<ny;j++){
			for(sum=0,k=0;k<nx;k++){
			     sum+=get_val(PSIx,ii,index6+k)*get_val(C,j,k);
			}
			put_val(PSI,ii,index5+j,sum);
		}
	}
	
    }
	minit(PSIx);
	dw = 0;
	/* 
     	 >>>>>>>>>>>>  Gradient (G = PSI_red*E_vector - D*theta_red)  <<<<<<<<<<<<<  
         */    
     	for(i=0; i<reduced; i++){
      		ii = (int)cvget(theta_index,i);
    		for(sum=0.0,k=skipstart; k<Nny; k++)
    			sum+=get_val(PSI,ii,k)*rvget(Evec,k);
    		cvput(G,i,sum - cvget(D,i)*cvget(theta_red,i));
      	}

    	/* 
     	 >>>>>>>>>> Mean square error part of Hessian (PSI_red*PSI_red') <<<<<<<<<<  
         */
    	for(i=0; i<reduced; i++){
    		ii = (int)cvget(theta_index,i);
    		for(j=i; j<reduced; j++){
      			jj = (int)cvget(theta_index,j);
      			for(sum=0.0,k=skipstart; k<Nny; k++)
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
  tmp1 = lambda - lambda_old;
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
   >>>>>>>>>>>>>       Compute network output y2(theta+h)          <<<<<<<<<<<<<< 
  */
  for(t=0;t<N;t++){
	mvmul(h1,W1_new,PHI,t);
	vtanh(y1,H_hidden,t,h1,H_hidden,0);
	vcopyi(y1,L_hidden,t,h1,L_hidden,0);
	
	mvmul(h2,W2_new,y1,t);
	vtanh(y2,H_output,t,h2,H_output,0);
	vcopyi(y2,L_output,t,h2,L_output,0);
	
	for(k=0;k<(nx-ny);k++){
		i=(int)cvget(nrowidx,k);
		y2->mat[i][t]+=get_val(PHI,i+1,t);
	}
	mvmul(Yhat,C,y2,t);                             /* Output prediction     */
	
	for(i=0;i<ny;i++){
		cvput(E,i,get_val(Y2,i,t)-cvget(Yhat,i));/* Prediction error     */
		rvput(Evec_new,t*ny+i,cvget(E,i));       /* Store E in Evec_new  */
	}
	if(t<N-1){
		for(i=0;i<nx;i++)
			put_val(PHI,i,t+1,get_val(y2,i,t));
		for(i=0;i<ny;i++)
			put_val(PHI,nx+nu+i,t+1,cvget(E,i));
	}
  }	
  for(SSE_new=0,t=skipstart;t<Nny;t++)
	SSE_new+=rvget(Evec_new,t)*rvget(Evec_new,t);   /* Sum of squared errors */
  for(tmp1=0,i=0;i<reduced;i++) tmp1+=cvget(theta_red_new,i)*cvget(theta_red_new,i)*cvget(D,i); 
	NSSE_new = (SSE_new+tmp1)/(2*N2);                               /* Value of cost function*/



  /*
   >>>>>>>>>>>>>>>>>>>>>>>>>>>       UPDATE  lambda     <<<<<<<<<<<<<<<<<<<<<<<<<<
   */
    lambda_old = lambda;
    for(tmp1=0,i=0;i<reduced;i++) tmp1+=cvget(h,i)*cvget(h,i)*(cvget(D,i)+lambda);
    L = sprod3(h,G) + tmp1;

    /* Decrease lambda if SSE has fallen 'sufficiently' */
    if(2*N2*(NSSE - NSSE_new) > (0.75*L)) lambda = lambda/2;
  
    /* Increase lambda if SSE has grown 'sufficiently'  */
    else if(2*N2*(NSSE-NSSE_new) <= (0.25*L)) lambda = 2*lambda;  


  /*
   >>>>>>>>>>>>>>>>>>>       UPDATES FOR NEXT ITERATION        <<<<<<<<<<<<<<<<<<<<
   */
    /* Update only if criterion has decreased */
    if(NSSE_new<NSSE)
    {
     critdif  = NSSE-NSSE_new;                           /* Criterion difference */
     for(i=0,gradmax=0.0,ptm1=G->mat[0];i<reduced;i++){  /* Maximum gradient     */
        sum = fabs(*(ptm1++));
        if(gradmax<sum)
           gradmax = sum;
     }
     gradmax/=N2;
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
     tmp = Evec; Evec = Evec_new; Evec_new = tmp;
     dw = 1;
     NSSE = NSSE_new;
     cvput(NSSEvec,iteration-1,NSSE);
     switch(trparms->infolevel){                            /* Print on-line inform */
       case 1:
          printf("# %i   W=%4.3e  critdif=%3.2e  maxgrad=%3.2e  paramdif=%3.2e\n",
                                                  iteration,NSSE,critdif,gradmax,paramdif);
          break;
       default:
          printf("iteration # %i   W = %4.3e\r",iteration,NSSE);
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
	*NSSEvecpp = mmake(1,1);
	(*NSSEvecpp)->row=0;
	(*NSSEvecpp)->col=0;
}
else{
	*NSSEvecpp = mmake(iteration,1);
	subvec(*NSSEvecpp,NSSEvec,0,iteration-1);
}
*iter = iteration;
*lam  = lambda;
mfree(L_hidden); mfree(H_hidden); mfree(L_output); mfree(H_output); mfree(Evec);
mfree(h1); mfree(h2); mfree(y1); mfree(y2); mfree(Y2); mfree(Ahat); mfree(C);
mfree(E); mfree(Evec_new); mfree(dxdy1); mfree(Khat); mfree(dy1de); mfree(dy1dx);
mfree(W1_new); mfree(W2_new); mfree(D);mfree(Dtmp); mfree(NSSEvec);mfree(Htmp);
mfree(theta); mfree(thtmp); mfree(theta_index); mfree(theta_red); mfree(theta_red_new);
mfree(PSI); mfree(PSIx); mfree(G); mfree(H); mfree(h), mfree(Yhat);
mfree(all); mfree(index0); mfree(index7);mfree(onesvec);mfree(tmp0); mfree(tmp2);
mfree(tmp3); mfree(index); mfree(index2); mfree(PHI);
mfree(rowidx); mfree(nrowidx); mfree(AKC);

printf("\n\nNetwork training ended.\n\n\n");
}
/*
  --------------------------------------------------------------------------------
  ----------------             END OF NETWORK TRAINING              --------------
  --------------------------------------------------------------------------------
*/


/*********************************************************************************
 *                                                                               *
 *                           G A T E W A Y   R O U T I N E                       *
 *                                                                               *
 *********************************************************************************/
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  /*
   >>>>>>>>>>>>>>>>>>           VARIABLE DECLARATIONS          <<<<<<<<<<<<<<<<<<<
   */
   matrix *NSSEvec;
   matrix *NetDef, *W1, *W2, *obsidx, *U, *Y;
   double *M, lambda;
   int iter, skip, ny, N, nu, nx, hidden, k, n, a, decays;
   trparmstruct *trparms;
   mxArray  *Matmatrix;
   char *infolevelstr[] = {"infolevel", "Infolevel", "INFOLEVEL", "InfoLevel"};
   char *maxiterstr[] = {"maxiter", "MAXITER", "Maxiter", "MaxIter"};
   char *critminstr[] = {"critmin", "Critmin", "CRITMIN", "CritMin"};
   char *crittermstr[] = {"critterm", "Critterm", "CRITTERM", "CritTerm"};
   char *gradtermstr[] = {"gradterm", "Gradterm", "GRADTERM", "GradTerm"};
   char *paramtermstr[] = {"paramterm", "Paramterm", "PARAMTERM", "ParamTerm"};
   char *Dstr[] = {"D", "d"};
   char *lambdastr[] = {"lambda", "Lambda", "LAMBDA"};
   char *skipstr[] = {"skip", "Skip", "SKIP"};


  /*
   >>>>>>>>>>>>>>>>      CHECK FOR PROPER NUMBER OF ARGUMENTS      <<<<<<<<<<<<<<<
   */
   if (nrhs!=8) mexErrMsgTxt("Wrong number of input arguments");
   else if (nlhs > 6) mexErrMsgTxt("Too many output arguments");
   nu = mxGetM(prhs[7]);   /* # of control signals */
   ny  = mxGetM(prhs[6]);  /* Rows of vector Y */
   if(nu<1) mexErrMsgTxt("Wrong dimension of input vector");
   N = mxGetN(prhs[6]);      /* Columns of vector Y */
   if(N!=mxGetN(prhs[7])) mexErrMsgTxt("U and Y should have the same number of columns!");
   hidden = mxGetN(prhs[0]); /* # of hidden units */
   if(mxGetM(prhs[0])!=2) mexErrMsgTxt("Error in architecture definition!");
   if(hidden<2) mexErrMsgTxt("Use at least two hidden units!");
         

  /*
   >>>>>>>>>>>>>>>>>     CONVERT INPUT ARGUMENTS TO SM FORMAT     <<<<<<<<<<<<<<<<
   */
  NetDef  = matstring2sm(prhs[0]);     /* Network architecture */
  nx = (int)(*mxGetPr(prhs[1]));
  Y    = mat2sm(prhs[6]);        /* Vector of observed outputs  */
  U = mat2sm(prhs[7]);           /* Vector of inputs            */

  /* Initialize pseudo-observability indices if obsidx passed as [] */
  if(mxGetM(prhs[4])==0 && mxGetN(prhs[4])==0){
  	obsidx=mmake(1,ny);
  	minitx(obsidx,1.0);
  	vput(obsidx,ny-1,(double)nx-ny+1);
  }
  else
  	obsidx = mat2sm(prhs[4]);


  /* Initialize weight matrices if passed as [] */
  if(mxGetM(prhs[2])==0 || mxGetN(prhs[2])==0 || mxGetM(prhs[3])==0\
                        || mxGetN(prhs[3])==0){
	W1 = mmake(hidden,nx+nu+ny+1);
   	W2 = mmake(nx,hidden+1);
      	W2->row = 0;   /* Hack telling that the weights should be initialized */ 
   }
   else{
   	if(mxGetM(prhs[2])!=hidden) mexErrMsgTxt("W1 has the wrong dimension");
   	if(mxGetN(prhs[2])!=nx+nu+ny+1) mexErrMsgTxt("W1 has the wrong dimension");
   	if(mxGetM(prhs[3])!=nx) mexErrMsgTxt("W2 has the wrong dimension");
   	if(mxGetN(prhs[3])!=hidden+1) mexErrMsgTxt("W2 has the wrong dimension");
   	W1 = mat2sm(prhs[2]);     /* Input-to-hidden layer weights */
   	W2 = mat2sm(prhs[3]);     /* Hidden-to-output layer weights */
   }
   
 trparms = (trparmstruct*)malloc(sizeof(trparmstruct)); 
 a = 5;
 if(mxGetN(prhs[a])!=0|| mxGetM(prhs[a])!=0) {
    /* INFOLEVEL */
    trparms->infolevel   = TRDINFOLEVEL;    
    for(n=0;n<4;n++){
        if ((Matmatrix=mxGetField(prhs[a], 0, infolevelstr[n]))!=NULL){
           trparms->infolevel=(int)(*mxGetPr(Matmatrix));
           break;
        }
    }

    /* MAXITER */
    trparms->maxiter   = TRDMAXITER;    
    for(n=0;n<4;n++){
        if ((Matmatrix=mxGetField(prhs[a], 0, maxiterstr[n]))!=NULL){
           trparms->maxiter=(int)(*mxGetPr(Matmatrix));
           break;
        }
    }

    /* CRITMIN */
    trparms->critmin   = TRDCRITMIN;    
    for(n=0;n<4;n++){
        if ((Matmatrix=mxGetField(prhs[a], 0, critminstr[n]))!=NULL){
           trparms->critmin=(double)(*mxGetPr(Matmatrix));
           break;
        }
    }

    
    /* CRITTERM */
    trparms->critterm   = TRDCRITTERM;    
    for(n=0;n<4;n++){
        if ((Matmatrix=mxGetField(prhs[a], 0, crittermstr[n]))!=NULL){
           trparms->critterm=(double)(*mxGetPr(Matmatrix));
           break;
        }
    }

    /* GRADTERM */
    trparms->gradterm   = TRDGRADTERM;    
    for(n=0;n<4;n++){
        if ((Matmatrix=mxGetField(prhs[a], 0, gradtermstr[n]))!=NULL){
           trparms->gradterm=(double)(*mxGetPr(Matmatrix));
           break;
        }
    }

    /* PARAMTERM */
    trparms->paramterm   = TRDPARAMTERM;    
    for(n=0;n<4;n++){
        if ((Matmatrix=mxGetField(prhs[a], 0, paramtermstr[n]))!=NULL){
           trparms->paramterm=(double)(*mxGetPr(Matmatrix));
           break;
        }
    }

    /* Lambda */
    trparms->lambda   = TRDLAMBDA;    
    for(n=0;n<3;n++){
        if ((Matmatrix=mxGetField(prhs[a], 0, lambdastr[n]))!=NULL){
           trparms->lambda=(double)(*mxGetPr(Matmatrix));
           break;
        }
    }

    /* Skip */
    trparms->skip   = TRDSKIP;    
    for(n=0;n<3;n++){
        if ((Matmatrix=mxGetField(prhs[a], 0, skipstr[n]))!=NULL){
           trparms->skip=(int)(*mxGetPr(Matmatrix));
           break;
        }
    }


    /* D */
    for(n=0;n<2;n++){
        if ((Matmatrix=mxGetField(prhs[a], 0, Dstr[n]))!=NULL){
           decays = mxGetM(Matmatrix)*mxGetN(Matmatrix);
           trparms->D         = mmake(1,decays);
           M    = mxGetPr(Matmatrix);
           for(n=0;n<decays;n++){
              rvput(trparms->D,n,M[n]);
           }
           break;
        }
    }
    if(Matmatrix==NULL){
       trparms->D         = mmake(1,1);
       put_val(trparms->D,0,0,TRDD);
    }
}
  else
  {
    trparms->infolevel = TRDINFOLEVEL;
    trparms->maxiter   = TRDMAXITER;
    trparms->critmin   = TRDCRITMIN;
    trparms->critterm  = TRDCRITTERM;
    trparms->gradterm  = TRDGRADTERM;
    trparms->paramterm = TRDPARAMTERM;
    trparms->D         = mmake(1,1);
    put_val(trparms->D,0,0,TRDD);
    trparms->lambda    = TRDLAMBDA;
    trparms->skip      = TRDSKIP;
  }


  /*
   >>>>>>>>>>>>>>>>>>>>>>         CALL THE C-ROUTINE         <<<<<<<<<<<<<<<<<<<<<
   */
  nnssif(&NSSEvec, &iter, &lambda, NetDef, nx, W1, W2, obsidx, trparms, Y, U);


  /*
   >>>>>>>>>>>>>>>>>>>         CREATE OUTPUT MATICES            <<<<<<<<<<<<<<<<<<
   */
  plhs[0] = mxCreateDoubleMatrix(getrows(W1),getcols(W1),mxREAL);
  plhs[1] = mxCreateDoubleMatrix(getrows(W2),getcols(W2),mxREAL);
  plhs[2] = mxCreateDoubleMatrix(getrows(obsidx),getcols(obsidx),mxREAL);
  plhs[3] = mxCreateDoubleMatrix(getrows(NSSEvec),getcols(NSSEvec),mxREAL);
  plhs[4] = mxCreateDoubleMatrix(1,1,mxREAL);
  plhs[5] = mxCreateDoubleMatrix(1,1,mxREAL);

  sm2mat(plhs[0],W1);
  sm2mat(plhs[1],W2);
  sm2mat(plhs[2],obsidx);
  sm2mat(plhs[3],NSSEvec);
  M = mxGetPr(plhs[4]); M[0] = (double)iter;
  M = mxGetPr(plhs[5]); M[0] = (double)lambda;

  /*
   >>>>>>>>>>>>>>>>>>>>        FREE ARGUMENT MATRICES        <<<<<<<<<<<<<<<<<<<<<
   */
  mfree(NetDef);
  mfree(U);
  mfree(Y);
  mfree(trparms->D);
  free(trparms);
}



