#include "stdafx.h"

#ifdef MKL
#include "MKLCG.hpp"
#include "mkl_rci.h"
#include "mkl_spblas.h"
#include "mkl_blas.h"
#include "mkl_service.h"

MKLCGSymmMatrix::MKLCGSymmMatrix(const int numEqns) :
		MKLPardisoSymmMatrix(numEqns)
{
	WhoClass("MKLCGSymmMatrix");
	WhoMethod("MKLCGSymmMatrix(const int)");
	Assert(row = new SortedList<int>[numberOfEquations]);
	coefs = 0; 	// Allocation done in method allocate
	rowindex = 0; 	// Allocation done in method allocate
	colindex = 0; 	// Allocation done in method allocate
	upperFormat = 'U';

	int i;
	for(i=0;i<numberOfEquations;i++)
		specifyNonZeroLocation(i,i); // Make sure there is a diagonal component

	DoInitializeRCICGParameters=true;
}

void MKLCGSymmMatrix::InitializeRCICGParameters()
{	
	mkl_set_num_threads(num_threads_for_solver);
	/*---------------------------------------------------------------------------*/
	/* Set the desired parameters: {default value}                               */
	/*---------------------------------------------------------------------------*/
	//Problem size.  Set by dcg_init.
	//ipar[0]=numberOfEquations;

	//Type of error output.  ({6} - output to screen, != 6 - output to 
	//dcg_errors.txt and dcg_check_warnings.txt)
	//ipar[1]=6;

	//Current stage of the RCI CG computations.  Set by dcg_init
	//ipar[2]=1;

	//Current iteration number.  Initial value is 0, set by dcg_init.
	//ipar[3]=0;

	//Maximum number of iterations.  Default is min{150, numberOfEquations}
	ipar[4]=MAX_NUMBER_ITERATIONS;

	//Output error messages (0 - No, {1} - Yes)
	//ipar[5]=1;

	//Output warning messages (0 - No, {1} - Yes)
	//ipar[6]=1;

	//Perform max iteration stopping test (0 - No, {1} - Yes)
	ipar[7]=1;

	//Perform residual stopping test ({0} - No, 1 - Yes)
	//ipar[8]=0;

	//Request user-defined stopping Test (0 - No, {1} - Yes)
	ipar[9]=1;

	//Perform preconditioned conjugate gradient method ({0} - No, 1 - Yes)
	ipar[10]=1;

	//Set the relative tolerance {1.0E-6}
	dpar[0]=SOLUTION_EPSILON;

	//Set the absolute tolerance {0.0e0}
	//dpar[1]=1.E-12;

	//Square norm of the initial residual, initially 0.0
	//dpar[2]= 0.0;

	//Service variable dpar[0]*dpar[2]+dpar[1], initially 0.0
	//dpar[3]=0.0;

	//Square Norm of the current residual, initially 0.0
	//dpar[4]=0.0;

	//Square norm of residual from the previous iteration step, initially 0.0
	//dpar[5]=0.0;

	//alpha parameter of the CG method, initially 0.0
	//dpar[6]=0.0;

	//beta parameter of the CG method, dpar[4]/dpar[5], intiailly 0.0
	//dpar[7]=0.0;
}

void MKLCGSymmMatrix::solve(double *b, double *initialGuess)
{
	WhoMethod("solve(double *)");

	// Check for positive definite
	CheckAndFixZeroDiagonal();
	
	MKL_INT rci_request;
	tmp = new double[4*numberOfEquations];
	double* x;
	double ConvCheck;
	int incVar = 1;

	Assert(x = new double [numberOfEquations]);
	Assert(KDiag = new double [numberOfEquations]);
	for(int i=0; i<numberOfEquations;++i){
		KDiag[i]=getaij(i,i);
	}

	/*---------------------------------------------------------------------------*/
	/* Initialize the solver                                                     */
	/*---------------------------------------------------------------------------*/
	dcg_init(&numberOfEquations,x,b,&rci_request,ipar,dpar,tmp);
	if (rci_request!=0){
		printf("dcg_init FAILED - returned the ERROR code %d", rci_request);
		exit(1);
	}

	InitializeRCICGParameters();

	/*---------------------------------------------------------------------------*/
	/* Check the correctness and consistency of the newly set parameters         */
	/*---------------------------------------------------------------------------*/
	dcg_check(&numberOfEquations,x,b,&rci_request,ipar,dpar,tmp);
	if (rci_request!=0){
		printf("dcg_check FAILED - returned the ERROR code %d", rci_request);
		exit(1);
	}

	// Set initial guess
	if(initialGuess==0){
		for(int i=0;i<numberOfEquations;i++) x[i]=0.0;
	}
	else {
		if(verbose >= Basic)	cerr << "I have the initial Guess." << endl;
		for(int i=0;i<numberOfEquations;i++) x[i] = initialGuess[i];
	}

	/*---------------------------------------------------------------------------*/
	/* Compute the Solution by RCI (P)CG solver                                  */
	/* Reverse Communications starts here                                        */
	/*---------------------------------------------------------------------------*/
	do{
		dcg(&numberOfEquations,x,b,&rci_request,ipar,dpar,tmp);
		switch(rci_request){
			/*---------------------------------------------------------------------------*/
			/* If rci_request=0, then the solution was found with the required precision */
			/*---------------------------------------------------------------------------*/
			case 0:
				dcg_get(&numberOfEquations,x,b,&rci_request,ipar,dpar,tmp,&itercount);
				delete [] tmp;
				for(int i=0;i<numberOfEquations;++i) b[i]=x[i];
				printf("\nNumber of iterations: %d\n",itercount);
				break;

			/*---------------------------------------------------------------------------*/
			/* If rci_request=1, then compute the vector A*tmp[0]					     */
			/* and put the result in vector tmp[numberOfEquations]                       */
			/*---------------------------------------------------------------------------*/
			case 1:
				mkl_dcsrsymv(&upperFormat, &numberOfEquations, coefs, rowindex, colindex, tmp, &tmp[numberOfEquations]);
				//printf("Iteration %d square norm = %g\n", ipar[3] ,dpar[4]);
				break;

			/*---------------------------------------------------------------------------*/
			/* If rci_request=2, then apply user-defined stopping test                   */
			/*---------------------------------------------------------------------------*/
			case 2:
					dcg_get(&numberOfEquations,x,b,&rci_request,ipar,dpar,tmp,&itercount);
					
					//If ||r||_2 / ||x||_2 <= SOLUTION_EPSILON then solution has converged within given tolerance
					ConvCheck = dnrm2(&numberOfEquations,&tmp[numberOfEquations*2],&incVar) / dnrm2(&numberOfEquations, x, &incVar);
					if( ConvCheck < dpar[0] ){
						for(int i=0;i<numberOfEquations;++i) b[i]=x[i];
						delete [] tmp;
						delete [] x;
						delete [] KDiag;
						printf("Final Iteration %d ||r||_2 / ||x||_2 = %g\n", itercount ,ConvCheck);
						rci_request = 0;
					}else{
						if(ipar[3]%100 == 0) printf("Iteration %d ||r||_2 / ||x||_2 = %g\n", ipar[3] ,ConvCheck);
					}
				break;

			/*---------------------------------------------------------------------------*/
			/* If rci_request=3, then apply the preconditioner matrix on                 */
			/* tmp[numberOfEquations*2] and store to tmp[numberOfEquations*3]            */
			/*---------------------------------------------------------------------------*/
			case 3:
				//printf("Preconditioner Not Currently Implemented!  Exiting!\n");
				//exit(1);
				//Perform Simple Jacobi Preconditioning (divide residual by diagonal of 
				//the coefficient matrix)
				//vdDiv( numberOfEquations, &(tmp[numberOfEquations*2]), KDiag, &(tmp[numberOfEquations*3]) );
				for(int i=0; i<numberOfEquations; ++i){
					tmp[numberOfEquations*3+i] = tmp[numberOfEquations*2+i] / KDiag[i];
				}
				break;
			/*---------------------------------------------------------------------------*/
			/* If rci_request=anything else, then dcg subroutine failed                  */
			/* to compute the Solution vector: Solution[numberOfEquations]               */
			/*---------------------------------------------------------------------------*/
			default:
				printf("dcg FAILED - returned the ERROR code %d", rci_request);
				exit(1);
		}
	}while(rci_request != 0);
}

void MKLCGSymmMatrix::setMAX_NUMBER_ITERATIONS(int max_iter)
{
	MAX_NUMBER_ITERATIONS=max_iter;
}

void MKLCGSymmMatrix::setSOLUTION_EPSILON(double sol_epsilon)
{
	SOLUTION_EPSILON=sol_epsilon;
}


#endif