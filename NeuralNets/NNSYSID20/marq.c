/*
 *     INCLUDE HEADERS
 */
#include <stdio.h>
#include <math.h>
#include <time.h>
#include "mex.h"
#include "matrix2.h"
#include "nnmisc.h"


/*********************************************************************************
 *                                                                               *
 *    MARQ                                                                       *
 *    ----                                                                       *
 *                                                                               *
 *    This is a C-version of the Matlab function marq.                           *
 *    Type 'help marq' from Matlab for information on                            *
 *    how to call this function.                                                 *
 *                                                                               *
 *                                                                               *
 *    Programmed by: Magnus Norgaard                                             *
 *    LastEditDate : Jan. 18, 2000                                               *
 *                                                                               *
 *********************************************************************************/


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
   matrix *PI_vector;
   matrix *NetDef, *W1, *W2, *PHI, *Y, *L_hidden, *H_hidden;
   double *M, lambda;
   int iter, hidden, inputs, outputs, a, n, decays;
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
   if (nrhs<5 || nrhs>6)
   { 
       mexErrMsgTxt("Wrong number of input arguments");
   }
   else if (nlhs > 5)
   {
       mexErrMsgTxt("Too many output arguments");
   }


  /*
   >>>>>>>>>>>>>>>>>     CONVERT INPUT ARGUMENTS TO SM FORMAT     <<<<<<<<<<<<<<<<
   */
  NetDef  = matstring2sm(prhs[0]);
  PHI     = mat2sm(prhs[3]);
  Y       = mat2sm(prhs[4]);

  inputs  = mxGetM(prhs[3]);
  outputs = mxGetM(prhs[4]);
  L_hidden = neuvector(NetDef,1,'L');      /* Location of linear hidden units      */
  H_hidden = neuvector(NetDef,1,'H');      /* Location of tanh hidden units        */  
  hidden   = L_hidden->row*L_hidden->col + H_hidden->row*H_hidden->col;
  if(mxGetM(prhs[1])==0 || mxGetN(prhs[1])==0 || mxGetM(prhs[2])==0\
                         || mxGetN(prhs[2])==0){
   	W1 = mmake(hidden,inputs+1);
   	W2 = mmake(outputs,hidden+1);
   	mrand(W1); smul(W1,W1,0.5);
   	mrand(W2); smul(W2,W2,0.5);
  }
  else{
   	if(mxGetM(prhs[1])!=hidden) mexErrMsgTxt("W1 has the wrong dimension");
   	if(mxGetN(prhs[1])!=inputs+1) mexErrMsgTxt("W1 has the wrong dimension");
   	if(mxGetM(prhs[2])!=outputs) mexErrMsgTxt("W2 has the wrong dimension");
   	if(mxGetN(prhs[2])!=hidden+1) mexErrMsgTxt("W2 has the wrong dimension");
   	W1 = mat2sm(prhs[1]);     /* Input-to-hidden layer weights */
   	W2 = mat2sm(prhs[2]);     /* Hidden-to-output layer weights */
  }
 trparms = (trparmstruct*)malloc(sizeof(trparmstruct)); 
 a = 5;
 if (nrhs==6){
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
  marqc(&PI_vector, &iter, &lambda, NetDef, W1, W2, PHI, Y, trparms);


  /*
   >>>>>>>>>>>>>>>>>>>         CREATE OUTPUT MATICES            <<<<<<<<<<<<<<<<<<
   */
  plhs[0] = mxCreateDoubleMatrix(getrows(W1),getcols(W1),mxREAL);
  plhs[1] = mxCreateDoubleMatrix(getrows(W2),getcols(W2),mxREAL);
  plhs[2] = mxCreateDoubleMatrix(getrows(PI_vector),getcols(PI_vector),mxREAL);
  plhs[3] = mxCreateDoubleMatrix(1,1,mxREAL);
  plhs[4] = mxCreateDoubleMatrix(1,1,mxREAL);

  sm2mat(plhs[0],W1);
  sm2mat(plhs[1],W2);
  sm2mat(plhs[2],PI_vector);
  M = mxGetPr(plhs[3]); M[0] = (double)iter;
  M = mxGetPr(plhs[4]); M[0] = (double)lambda;

  /*
   >>>>>>>>>>>>>>>>>>>>        FREE ARGUMENT MATRICES        <<<<<<<<<<<<<<<<<<<<<
   */
  mfree(NetDef); mfree(PHI); mfree(Y); mfree(L_hidden); mfree(H_hidden);
  mfree(trparms->D);
  free(trparms); 
}

