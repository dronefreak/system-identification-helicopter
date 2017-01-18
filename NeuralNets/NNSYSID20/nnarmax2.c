/*
 *     INCLUDE HEADERS
 */
#include <stdio.h>
#include <math.h>
#include "mex.h"
#include "matrix2.h"
#include "nnmisc.h"

void nnarmax2(matrix**, int*, double*, matrix*, matrix*, matrix*, matrix*, trparmstruct*,\
 matrix*, matrix*);



/*********************************************************************************
 *                                                                               *
 *    NNARMAX2                                                                   *
 *    --------                                                                   *
 *                                                                               *
 *    This is a CMEX-version of the Matlab function nnarmax2.                    *
 *    Type 'help nnarmax2' from Matlab for information on                        *
 *    how to call this function.                                                 *
 *                                                                               *
 *                                                                               *
 *    Programmed by: Magnus Norgaard                                             *
 *    LastEditDate : Jan. 21, 2000                                               *
 *                                                                               *
 *********************************************************************************/
void nnarmax2(matrix **NSSEvecpp, int *iter, double *lam,\
	matrix *NetDef, matrix *NN, matrix *W1, matrix *W2, trparmstruct *trparms,\
	matrix *Y, matrix *U)
{
/*
----------------------------------------------------------------------------------- 
---------------              VARIABLE DECLARATIONS                    ------------- 
----------------------------------------------------------------------------------- 
*/ 
register i, j, k, t; 
int outputs, N, Nout, layers, dummy, hidden, inputs, iteration; 
int parameters1, parameters2, parameters, reduced, index1, ii, jj; 
int lhids, hhids, louts, houts, index11;
int Ndat, N2, na, nc, nu, nab, nabc, nmax, index5, dummy2, skip;
double lambda, SSE, SSE_new, NSSE, NSSE_new, L, tmp1, sum, dummy3; 
double critdif, gradmax, paramdif, *ptm1, *ptm2, lambda_old;
char dw; 
matrix *L_hidden, *H_hidden, *L_output, *H_output, *h1, *h2, *y1, *y2; 
matrix *E, *E_new, *W1_new, *W2_new, *PHI, *D, *Dtmp; 
matrix *NSSEvec, *miter, *tmp, *Htmp, *W1tmp; 
matrix *theta, *thtmp, *theta_index, *theta_red, *theta_red_new, *PSI, *G, *H, *h; 
matrix *all, *index0, *index7, *onesvec, *tmp0, *tmp2, *tmp3, *index, *index2;
matrix *nb, *nk, *dy2dy1, *dy2de, *dy1de, *dy2de_vec, *Y2, *dummy1;
trparmstruct *trd;
 
 
/* 
----------------------------------------------------------------------------------- 
---------------             NETWORK INITIALIZATIONS                   ------------- 
----------------------------------------------------------------------------------- 
 */
Ndat      = getcols(Y);                  /* # of data                            */
na        = vget(NN,0);                  /* Past predictions used as inputs      */
nu	  = getrows(U);			 /* # of input signals                   */ 
nc        = vget(NN,nu+1);               /* Past prediction errors used as input */
if(nu!=0){
	nb = mmake(1,nu);                /* Past controls used as inputs	 */
            subvec(nb,NN,1,nu);
	nk = mmake(1,nu);                /* Time delays                          */
	    subvec(nk,NN,nu+2,2*nu+1);
}
nmax      = na;		                 /* Oldest signal used as input          */
if(nmax<nc) nmax=nc;
for(k=0;k<nu;k++){
  i=rvget(nb,k)+rvget(nk,k)-1;
  if(nmax<i) nmax=i;
}
skip      = trparms->skip;
N         = Ndat - nmax;                 /* Size of training set                 */
N2        = N-skip;
nab       = na; 			 /* na+nb                                */
for(k=0;k<nu;k++) nab=nab+rvget(nb,k);
nabc      = nab+nc;			 /* na+nb+nc                             */
Y2        = mmake(1,N);                  /* Observed outputs used for training   */
hidden    = getcols(NetDef);             /* # of hidden units                    */
inputs    = nabc;                 	 /* Number of inputs to network          */
outputs   = 1;		                 /* Always one outputs                   */
Nout      = N*outputs;                   /* N*outputs                            */
L_hidden  = neuvector(NetDef,1,'L');     /* Location of linear hidden units      */
H_hidden  = neuvector(NetDef,1,'H');     /* Location of tanh hidden units        */ 
L_output  = neuvector(NetDef,2,'L');     /* Location of linear output units      */ 
H_output  = neuvector(NetDef,2,'H');     /* Location of tanh output units        */ 
lhids     = getrows(L_hidden);           /* # of linear hidden units             */ 
hhids     = getrows(H_hidden);           /* # of tanh hidden units               */ 
louts     = getrows(L_output);           /* # of linear output units             */  
houts     = getrows(H_output);           /* # of tanh output units               */
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
E         = mmake(outputs,N);            /* Prediction error matrix              */
E_new     = mmake(outputs,N);            /* A priori E                           */ 
index     = mmake(hidden,1);             /* Index vector outputs*(hidden+1)+...  */
for(k=0;k<hidden;k++) cvput(index,k,(double)(outputs*(hidden+1)+k*(inputs+1))); 
index2    = mmake(N,1);                  /* Index vector (0:N-1)*outputs         */
for(k=0;k<N;k++) cvput(index2,k,(double)k*outputs); 
iteration = 1;                           /* Initialize iteration counter         */
dw        = 1;                           /* Flag telling that the weights are new*/ 
parameters1= hidden*(inputs+1);          /* # of input-to-hidden weights         */
parameters2= outputs*(hidden+1);         /* # of hidden-to-output weights        */ 
parameters = parameters1+parameters2;    /* Total # of weights                   */ 

/*
 >>>>>>>>>>>>>>>>>>>>  CONSTRUCT THE REGRESSION MATRIX PHI   <<<<<<<<<<<<<<<<<<<<<
 */
PHI = mmake(nabc+1,N);      	         /* Matrix of input vectors (incl. bias) */
minit(PHI);
mat2mat(PHI,nabc,0,onesvec);
for(k=0;k<na;k++){
	for(i=0;i<Ndat-nmax;i++) mput(PHI,k,i,vget(Y,i+nmax-k-1));
}
index5 = na;                             /* Insert controls in PHI               */
for(i=0;i<nu;i++){
	for(k=0;k<vget(nb,i);k++){
		for(j=0;j<Ndat-nmax;j++){
			mput(PHI,index5+k,j,mget(U,i,nmax+j-k-vget(nk,i)));
		}
	}
	index5=index5+vget(nb,i);
}
for(t=0;t<N;t++) rvput(Y2,t,rvget(Y,t+nmax));

/*
 >>>>>>>>>>>>>>>>>  INITIALIZE WEIGHTS WITH NNARX IF NECESSARY   <<<<<<<<<<<<<<<<<<
 */
if(getrows(W2)==0){
	W2->row=1;
	mrand(W1); smul(W1,W1,0.025);
   	mrand(W2); smul(W2,W2,0.5);
   	W1tmp = mmake(hidden,nab+1);
   	mrand(W1tmp); smul(W1tmp,W1tmp,0.5);
	PHI->row=nab;
   trd = (trparmstruct*)malloc(sizeof(trparmstruct)); 
   trd->infolevel = TRDINFOLEVEL;
   trd->maxiter   = 100;
   trd->critmin   = TRDCRITMIN;
   trd->critterm  = TRDCRITTERM;
   trd->gradterm  = TRDGRADTERM;
   trd->paramterm = TRDPARAMTERM;
   trd->D         = mmake(1,1);
   put_val(trd->D,0,0,TRDD);
   trd->lambda    = TRDLAMBDA;
   trd->skip      = TRDSKIP;
	marqc(&dummy1, &dummy2, &dummy3, NetDef, W1tmp, W2, PHI, Y2, trparms);
	PHI->row=inputs+1;
   	mat2mat(W1,0,0,W1tmp);
   mfree(trd->D);
   free(trd);
	mfree(dummy1); mfree(W1tmp);
}

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
dy2de     = mmake(nc,N);                 /* Der. of output wrt. past pred. errors*/
dy1de     = mmake(hidden,nc);            /* Der.of hid. outp. wrt. past pred. err*/
dy2dy1    = mmake(1,hidden);             /* Der. of outp. wrt. hidden outp.      */
dy2de_vec = mmake(1,nc);		 /* For temp. results   		 */
PSI       = mmake(parameters,Nout);      /* Der. of each output wrt. each weight */
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
	
	rvput(E,t,rvget(Y2,t)-rvget(y2,t));          /* Prediction error        */	
	j=nc;
	if(N-t-1<nc) j=N-t-1;
	for(i=1;i<=j;i++){
		put_val(PHI,nab+i-1,t+i,rvget(E,t));
	}
}
for(SSE=0,t=skip;t<N;t++) SSE+=rvget(E,t)*rvget(E,t);/* Sum of squared errors   */
for(tmp1=0,i=0;i<reduced;i++) tmp1+=cvget(theta_red,i)*cvget(theta_red,i)*cvget(D,i); 
NSSE = (SSE+tmp1)/(2*N2);                            /* Value of cost function  */

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


    /* ===========  Elements corresponding to the tanh output units   ===========*/
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
      

    /* 
     >>>>>>>>>>>>>>>>>>>>        Linearize network           <<<<<<<<<<<<<<<<<<<<<  
    */
    for(t=0;t<N;t++){
    	/*-- Derivative of output wrt. hidden outputs --*/
    	if(louts==1) for(k=0;k<hidden;k++) rvput(dy2dy1,k,rvget(W2,k));
    	else for(k=0;k<hidden;k++) rvput(dy2dy1,k,rvget(W2,k)*(1-\
    	                             		rvget(y2,t)*rvget(y2,t)));

      	/*-- Partial deriv. of output from each hidden unit wrt. past net outp. --*/
      	for(j=0;j<lhids;j++){
      		i=(int)cvget(L_hidden,j);
      		for(k=nab;k<nabc;k++) put_val(dy1de,i,k-nab,\
      		        get_val(W1,i,k));
      	}
      	for(j=0;j<hhids;j++){
      		i=(int)cvget(H_hidden,j);
      		for(k=nab;k<nabc;k++) put_val(dy1de,i,k-nab,\
      			get_val(W1,i,k)*(1-get_val(y1,i,t)*get_val(y1,i,t)));
      	}

     	/*--Partial derivative of net output w.r.t. past net outputs --*/
     	mmul(dy2de_vec,dy2dy1,dy1de);
     	for(i=0;i<nc;i++) put_val(dy2de,i,t,rvget(dy2de_vec,i));
    }


    /* 
     >>>>>>>>>>>>>>>>>>>     Filter partial derivatives        <<<<<<<<<<<<<<<<<<<<  
    */
	for(t=0;t<N;t++){
		j=nc;
		if(t<nc) j=t;
		for(k=1;k<=j;k++){
			for(i=0;i<reduced;i++) {
			  ii =(int)cvget(theta_index,i);
			  PSI->mat[ii][t] -= get_val(dy2de,k-1,t)*get_val(PSI,ii,t-k);
			  }
		}
	}

	dw = 0;
	/* 
     	 >>>>>>>>>>>>  Gradient (G = PSI_red*E_vector - D*theta_red)  <<<<<<<<<<<<<  
         */    
     	for(i=0; i<reduced; i++){
      		ii = (int)cvget(theta_index,i);
    		for(sum=0.0,k=skip; k<N; k++) sum+=get_val(PSI,ii,k)*rvget(E,k);
    		cvput(G,i,sum - cvget(D,i)*cvget(theta_red,i));
      	}

    	/* 
     	 >>>>>>>>>> Mean square error part of Hessian (PSI_red*PSI_red') <<<<<<<<<<  
         */
    	for(i=0; i<reduced; i++){
    		ii = (int)cvget(theta_index,i);
    		for(j=i; j<reduced; j++){
      			jj = (int)cvget(theta_index,j);
      			for(sum=0.0,k=skip; k<N; k++)
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
   >>>>>>>>>>>>>       Compute network output y2(theta+h)          <<<<<<<<<<<<<< 
  */
  for(t=0;t<N;t++){
	mvmul(h1,W1_new,PHI,t);
	vtanh(y1,H_hidden,t,h1,H_hidden,0);
	vcopyi(y1,L_hidden,t,h1,L_hidden,0);
	
	mvmul(h2,W2_new,y1,t);
	vtanh(y2,H_output,t,h2,H_output,0);
	vcopyi(y2,L_output,t,h2,L_output,0);
	
	rvput(E_new,t,rvget(Y2,t)-rvget(y2,t));          /* Prediction error        */	
	j=nc;
	if(N-t-1<nc) j=N-t-1;
	for(i=1;i<=j;i++){
		put_val(PHI,nab+i-1,t+i,rvget(E_new,t));
	}
  }
  for(SSE_new=0,t=skip;t<N;t++)                        /* Sum of squared errors */
              SSE_new+=rvget(E_new,t)*rvget(E_new,t);
  for(tmp1=0.0,i=0;i<reduced;i++) tmp1+=cvget(theta_red_new,i)*cvget(theta_red_new,i)*cvget(D,i); 
  NSSE_new = (SSE_new+tmp1)/(2*N2);                    /* Value of cost function*/


  /*
   >>>>>>>>>>>>>>>>>>>>>>>>>>>       UPDATE  lambda     <<<<<<<<<<<<<<<<<<<<<<<<<<
   */
    lambda_old = lambda;
    for(tmp1=0.0,i=0;i<reduced;i++) tmp1+=cvget(h,i)*cvget(h,i)*(cvget(D,i)+lambda);
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
     tmp = E; E = E_new; E_new = tmp;
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
else {
	*NSSEvecpp = mmake(iteration,1);
	subvec(*NSSEvecpp,NSSEvec,0,iteration-1);
}
*iter = iteration;
*lam  = lambda;
mfree(L_hidden); mfree(H_hidden); mfree(L_output); mfree(H_output);
mfree(h1); mfree(h2); mfree(y1); mfree(y2); mfree(Y2);
mfree(E); mfree(E_new); mfree(dy2dy1); mfree(dy2de); mfree(dy1de); mfree(dy2de_vec);
mfree(W1_new); mfree(W2_new); mfree(D);mfree(Dtmp); mfree(NSSEvec);mfree(Htmp);
mfree(theta); mfree(thtmp); mfree(theta_index); mfree(theta_red); mfree(theta_red_new);
mfree(PSI); mfree(G); mfree(H); mfree(h);
mfree(all); mfree(index0); mfree(index7);mfree(onesvec);mfree(tmp0); mfree(tmp2);
mfree(tmp3); mfree(index); mfree(index2); mfree(PHI);
if(nu!=0){
	mfree(nb); mfree(nk);
}
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
   matrix *NSSEvec, *NN, *nb;
   matrix *NetDef, *W1, *W2, *U, *Y;
   double *M, lambda;
   int iter, skip, m, N, na, nc, nabc, nu, hidden, k, n, a, decays;
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
   if (nrhs<6 || nrhs>7) mexErrMsgTxt("Wrong number of input arguments");
   else if (nlhs > 5) mexErrMsgTxt("Too many output arguments");
   if(nrhs==7)
     nu = mxGetM(prhs[6]);   /* # of control signals */
   else
     nu=0;
   m  = mxGetM(prhs[5]);     /* Rows of vector Y */
   if(m!=1) mexErrMsgTxt("Wrong dimension of vector of observed outputs");
   N = mxGetN(prhs[5]);      /* Columns of vector Y */
   if(nu!=0) {
   if(N!=mxGetN(prhs[6])) mexErrMsgTxt("U and Y should have the same number of columns!");
   }
   hidden = mxGetN(prhs[0]); /* # of hidden units */
   if(mxGetM(prhs[0])!=2) mexErrMsgTxt("Error in architecture definition!");
   if(mxGetM(prhs[1])!=1) mexErrMsgTxt("NN should only have 1 row!");
   if(mxGetN(prhs[1])!=2*nu+2) mexErrMsgTxt("Mismatch between U and NN");
         

  /*
   >>>>>>>>>>>>>>>>>     CONVERT INPUT ARGUMENTS TO SM FORMAT     <<<<<<<<<<<<<<<<
   */
  NetDef  = matstring2sm(prhs[0]);/* Network architecture */
  NN      = mat2sm(prhs[1]);     /* Regressor structure  */
  Y    = mat2sm(prhs[5]);        /* Vector of observed outputs  */
  if(nu!=0) U = mat2sm(prhs[6]); /* Vector of inputs            */
  else{
  	U=mmake(1,1);
  	U->row=0;
  	U->col=0;
  }
  

  /*
   >>>>>>>>>>>>>>>>      CHECK FOR PROPER NUMBER OF ARGUMENTS      <<<<<<<<<<<<<<<
   */
   na        = vget(NN,0);   /* Past outputs used as input           */
   nc        = vget(NN,nu+1);/* Past prediction errors used as input */
   if(nu!=0){
   	nb = mmake(1,nu);    /* Past controls used as input          */
   	subvec(nb,NN,1,nu);
   }
   nabc      = na+nc;        /* na+nb+nc                              */
   for(k=0;k<nu;k++) nabc=nabc+rvget(nb,k);
   
   /* Initialize weight matrices if passed as [] */
   if(mxGetM(prhs[2])==0 || mxGetN(prhs[2])==0 || mxGetM(prhs[3])==0\
                         || mxGetN(prhs[3])==0){
   	W1 = mmake(hidden,nabc+1);
   	W2 = mmake(1,hidden+1);
      	W2->row = 0;   /* Hack telling that the weights should be initialized */ 

   }
   else{
   	if(mxGetM(prhs[2])!=hidden) mexErrMsgTxt("W1 has the wrong dimension");
   	if(mxGetN(prhs[2])!=nabc+1) mexErrMsgTxt("W1 has the wrong dimension");
   	if(mxGetM(prhs[3])!=1) mexErrMsgTxt("W2 has the wrong dimension");
   	if(mxGetN(prhs[3])!=hidden+1) mexErrMsgTxt("W2 has the wrong dimension");
   	W1 = mat2sm(prhs[2]);     /* Input-to-hidden layer weights */
   	W2 = mat2sm(prhs[3]);     /* Hidden-to-output layer weights */
   }

 trparms = (trparmstruct*)malloc(sizeof(trparmstruct)); 
 a = 4;
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
  nnarmax2(&NSSEvec, &iter, &lambda, NetDef, NN, W1, W2, trparms, Y, U);


  /*
   >>>>>>>>>>>>>>>>>>>         CREATE OUTPUT MATICES            <<<<<<<<<<<<<<<<<<
   */
  plhs[0] = mxCreateDoubleMatrix(getrows(W1),getcols(W1),mxREAL);
  plhs[1] = mxCreateDoubleMatrix(getrows(W2),getcols(W2),mxREAL);
  plhs[2] = mxCreateDoubleMatrix(getrows(NSSEvec),getcols(NSSEvec),mxREAL);
  plhs[3] = mxCreateDoubleMatrix(1,1,mxREAL);
  plhs[4] = mxCreateDoubleMatrix(1,1,mxREAL);

  sm2mat(plhs[0],W1);
  sm2mat(plhs[1],W2);
  sm2mat(plhs[2],NSSEvec);
  M = mxGetPr(plhs[3]); M[0] = (double)iter;
  M = mxGetPr(plhs[4]); M[0] = (double)lambda;

  /*
   >>>>>>>>>>>>>>>>>>>>        FREE ARGUMENT MATRICES        <<<<<<<<<<<<<<<<<<<<<
   */
  mfree(NetDef);
  mfree(NN);
  mfree(U);
  if(nu!=0) mfree(nb);
  mfree(Y);
  mfree(trparms->D);
  free(trparms);
}



