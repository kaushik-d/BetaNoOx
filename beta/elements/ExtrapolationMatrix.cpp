#include "stdafx.h"
#include "ExtrapolationMatrix.hpp"
#include "utility/utility.h"
#include "math/gaussQuadrature/gaussPointList.hpp"
#include "mesh/node.hpp"

//Prototypes
void CompletePolynomialMonomialVector(const int &numDims, 
                                      const double * coords, 
                                      const int &extrapOrder,
                                      std::vector<double> &GVec);

int CompletePolynomialGetNumMonomials(const int numDims,
                                      const int extrapOrder);

void TensorProductMonomialVector(const int &numDims, 
                                 const double * coords, 
                                 const int &extrapOrder,
                                 std::vector<double> &GVec);

int TensorProductGetNumMonomials(const int numDims,
                                 const int extrapOrder);

//Constructor for quad data extrapolation to an element's nodes, performed in the normalized space
ExtrapolationMatrix::ExtrapolationMatrix(const int &nNodes, 
                                         const int &nDims, 
                                         const GaussPointList &gpl, 
                                         const double *masterNodeCoords,  
                                         const int &extrapOrder,
                                         const int &ExtrapMonomialsSource):
Matrix(nNodes,gpl.totalNumIPs),
numUnknownCoords(nNodes), 
numKnownCoords(gpl.totalNumIPs), 
dimensionality(nDims),
order(extrapOrder),
desiredOrder(extrapOrder)
{
	/*
	Polynomials are used to approximate the field variable across the element.
	The interpolation order will either be that requested in the function
	argument or the highest order which is permitted by the number of quadrature
	points in the element.
	*/

    switch(ExtrapMonomialsSource){
    case CompletePolynomial:
        MonomialVector = &CompletePolynomialMonomialVector;
        GetNumMonomials = &CompletePolynomialGetNumMonomials;
        break;
    case TensorProduct:
        MonomialVector = &TensorProductMonomialVector;
        GetNumMonomials = &TensorProductGetNumMonomials;
        break;
    }

	int numMonomials = CheckOrder();
    vector <double> G;

	// Calculate H:
	Matrix H(numKnownCoords,numMonomials);
	double coord[3]; // maximum size of a coordinate
	switch(dimensionality){
		case 1:
			for(int i=0;i<numKnownCoords;++i){
				coord[0]=gpl.xiList[i];
				MonomialVector(dimensionality, coord, order, G);
				for(int j=0;j<numMonomials;++j){
					H[i][j]=G[j];
				}
			}
			break;
		case 2:
			for(int i=0;i<numKnownCoords;++i){
				coord[0]=gpl.xiList[i];
				coord[1]=gpl.etaList[i];
				MonomialVector(dimensionality, coord, order, G);
				for(int j=0;j<numMonomials;++j){
					H[i][j]=G[j];
				}
			}
			break;
		case 3:
			for(int i=0;i<numKnownCoords;++i){
				coord[0]=gpl.xiList[i];
				coord[1]=gpl.etaList[i];
				coord[2]=gpl.zetaList[i];
				MonomialVector(dimensionality, coord, order, G);
				for(int j=0;j<numMonomials;++j){
					H[i][j]=G[j];
				}
			}
			break;
		default:
			std::cout << "ExtrapolationMatrix.cpp: You can only have 1, 2, or 3 dimension problems" << std::endl;
			exit(1);
	}
	//H.print("H Matrix");

	//Make Q Matrix
	Matrix Q(numUnknownCoords,numMonomials);
	for(int i=0;i<numUnknownCoords;++i){
		MonomialVector(dimensionality,&masterNodeCoords[i*dimensionality],order,G);
		for(int j=0;j<numMonomials;++j){
			Q[i][j]=G[j];
		}
	}

	//Q.print("Q Matrix");
	
    PopulateMatrixFromHandQ(H,Q);

	//temp.print("Temp Extrapolation Matrix");
	//print("Actual Extrapolation Matrix");
}

//============================================================================
// A constructor for a more general mapping from a vector of nodes where data is known
// to a GaussPointList of QuadPoints where data is not known.
 ExtrapolationMatrix::ExtrapolationMatrix(const int &numDims,
                                          const std::vector<Node> &KnownPoints,
                                          const GaussPointList &UnknownGaussPoints,
                                          const int &extrapOrder,
                                          const int &ExtrapMonomialsSource):
Matrix(UnknownGaussPoints.totalNumIPs,KnownPoints.size()),
numUnknownCoords(UnknownGaussPoints.totalNumIPs), 
numKnownCoords(KnownPoints.size()), 
dimensionality(numDims),
order(extrapOrder),
desiredOrder(extrapOrder)
{
	/*
	Polynomials are used to approximate the field variable across the element.
	The interpolation order will either be that requested in the function
	argument or the highest order which is permitted by the number of quadrature
	points in the element.
	*/

    switch(ExtrapMonomialsSource){
    case CompletePolynomial:
        MonomialVector = &CompletePolynomialMonomialVector;
        GetNumMonomials = &CompletePolynomialGetNumMonomials;
        break;
    case TensorProduct:
        MonomialVector = &TensorProductMonomialVector;
        GetNumMonomials = &TensorProductGetNumMonomials;
        break;
    }

    int numMonomials = CheckOrder();
    vector <double> G;

	// Calculate H:
	Matrix H(numKnownCoords,numMonomials);
	double coord[3]; // maximum size of a coordinate
	switch(dimensionality){
		case 1:
			for(int i=0;i<numKnownCoords;++i){
				coord[0]=KnownPoints[i].x;
				MonomialVector(dimensionality, coord, order, G);
				for(int j=0;j<numMonomials;++j){
					H[i][j]=G[j];
				}
			}
			break;
		case 2:
			for(int i=0;i<numKnownCoords;++i){
				coord[0]=KnownPoints[i].x;
				coord[1]=KnownPoints[i].y;
				MonomialVector(dimensionality, coord, order, G);
				for(int j=0;j<numMonomials;++j){
					H[i][j]=G[j];
				}
			}
			break;
		case 3:
			for(int i=0;i<numKnownCoords;++i){
				coord[0]=KnownPoints[i].x;
				coord[1]=KnownPoints[i].y;
				coord[2]=KnownPoints[i].z;
				MonomialVector(dimensionality, coord, order, G);
				for(int j=0;j<numMonomials;++j){
					H[i][j]=G[j];
				}
			}
			break;
		default:
			std::cout << "ExtrapolationMatrix.cpp: You can only have 1, 2, or 3 dimension problems" << std::endl;
			exit(1);
	}
	//H.print("H Matrix");

	//Make Q Matrix
	Matrix Q(numUnknownCoords,numMonomials);
    	switch(dimensionality){
		case 1:
			for(int i=0;i<numUnknownCoords;++i){
				coord[0]=UnknownGaussPoints.xiList[i];
				MonomialVector(dimensionality, coord, order, G);
				for(int j=0;j<numMonomials;++j){
					Q[i][j]=G[j];
				}
			}
			break;
		case 2:
			for(int i=0;i<numUnknownCoords;++i){
				coord[0]=UnknownGaussPoints.xiList[i];
				coord[1]=UnknownGaussPoints.etaList[i];
				MonomialVector(dimensionality, coord, order, G);
				for(int j=0;j<numMonomials;++j){
					Q[i][j]=G[j];
				}
			}
			break;
		case 3:
			for(int i=0;i<numUnknownCoords;++i){
				coord[0]=UnknownGaussPoints.xiList[i];
				coord[1]=UnknownGaussPoints.etaList[i];
				coord[2]=UnknownGaussPoints.zetaList[i];
				MonomialVector(dimensionality, coord, order, G);
				for(int j=0;j<numMonomials;++j){
					Q[i][j]=G[j];
				}
			}
			break;
		default:
			std::cout << "ExtrapolationMatrix.cpp: You can only have 1, 2, or 3 dimension problems" << std::endl;
			exit(1);
	}
	//Q.print("Q Matrix");
	
	PopulateMatrixFromHandQ(H,Q);

	//temp.print("Temp Extrapolation Matrix");
	//print("Actual Extrapolation Matrix");
}
//============================================================================
int ExtrapolationMatrix::CheckOrder(){
    //Ensure that there are fewer monomials than coordinates where the value is known
    //Also return the number of monomials
    
	//Get the interpolation order low enough for the number of quadrature points,
	//regardless of the requested extrapolation order (if there aren't enough 
	//quad points, the H matrix won't invert)
	int numMonomials = GetNumMonomials(dimensionality,order);
    // Reduce the order of the extrapolation until there are fewer monomials than known values
	while(numMonomials>numKnownCoords){
		--order;
		numMonomials = GetNumMonomials(dimensionality,order);
        //std::cout<<"WARNING - ExtrapolationMatrix::CheckOrder:"<<std::endl;
        //std::cout<<"Specified extrapolation order too high for number of points where value is known." << std::endl;
        //std::cout<<"Automatically reducing order..."<<std::endl;
	}
    return numMonomials;
}
//============================================================================
void ExtrapolationMatrix::PopulateMatrixFromHandQ(const Matrix &H, 
                                                  const Matrix &Q){
    //There is likely room for efficiency improvements here...
    Matrix HTransposed(H);
	HTransposed.transpose();
	Matrix HTHInv(numKnownCoords,numKnownCoords);
	HTHInv = HTransposed * H;
	HTHInv.invert();
	Matrix temp(numUnknownCoords, numKnownCoords);
	temp = Q * HTHInv * HTransposed;
	for(int i=0; i<numUnknownCoords; ++i){
		for(int j=0; j<numKnownCoords; ++j){
			(*this)(i,j) = temp(i,j);
        }
    }
}

//============================================================================
// NON-CLASS FUNCTIONS
//============================================================================
int CompletePolynomialGetNumMonomials(const int numDims,
                                      const int extrapOrder)
{
    //Returns the number of monomials for the given dimensionality and order.
    int NumMonomials;
    int temp = extrapOrder+1;
    switch(numDims){
        case 1:
            NumMonomials = temp;
            break;
        case 2:
            NumMonomials = temp*(temp+1)/2;
            break;
        case 3:
            NumMonomials = temp*(temp+1)*(temp+2)/6;
            break;
        default:
            std::cout << "ERROR! - GetNumMonomials() only works for 1,2,or 3 dimensions." << std::endl;
            std::cout << numDims << " dimensions were requested.  Exiting..." << std::endl;
            exit(1);
    }
    return NumMonomials;
}
//============================================================================
void  CompletePolynomialMonomialVector(const int &numDims, 
                                       const double * coords, 
                                       const int &extrapOrder,
                                       std::vector<double> &GVec)
{
	//returns a vector of monomials in the complete polynomial up to the input extrapOrder evaluated at the input coordinate
    int NumMonomials = CompletePolynomialGetNumMonomials(numDims,extrapOrder);
    GVec.resize(NumMonomials);

    //1D is very straightforward...
	if(numDims == 1){
		GVec[0]=1.0;
        //std::cout << GVec.at(0) << "\t";
		for(int i=1;i<extrapOrder+1;++i){
			GVec[i] = coords[0]*GVec[i-1];
            //std::cout << GVec.at(i) << "\t";
		}
	}

	if(numDims == 2){
        double** PascalsTriangle = new double* [extrapOrder+1];
        int LastIndex=0;
        for(int i=0;i<extrapOrder+1;++i){
            PascalsTriangle[i] = &GVec[LastIndex];
            LastIndex+=i+1;
        }
        PascalsTriangle[0][0]=1.0;
        for(int i=1;i<extrapOrder+1;++i){
            PascalsTriangle[i][0]=coords[0]*PascalsTriangle[i-1][0];
            for(int j=1;j<i;++j){
                PascalsTriangle[i][j]=PascalsTriangle[i-1][j-1]*PascalsTriangle[i-1][j];
            }
            PascalsTriangle[i][i] = coords[1]*PascalsTriangle[i-1][i-1];
        }

        //Check
        //for(int i=0;i<extrapOrder+1;++i){
        //    for(int j=0;j<i+1;++j){
        //        std::cout << PascalsTriangle[i][j] << "  ";
        //    }
        //    std::cout << std::endl;
        //}

        delete [] PascalsTriangle;
    }

    if(numDims == 3){
        // For a 3D Pascal's Pyramid, generate it like this...
        // (This example goes to 3rd order monomials)
        // First Index:  1                     2               3         4  
        //                                                                  
        //               1                     z              z^2       z^3 
        //            x     y               zx   zy       z^2x   z^2y       
        //         x^2  xy   y^2        zx^2  zxy  zy^2                     
        //      x^3  x^2y  xy^2  y^3                                        
        //                                                                        

        //Memory Allocation - allocate in blocks to avoid a bunch of mallocs
        double *** PascalsPyramid = new double** [extrapOrder+1];
        double ** DubPtrPtr = new double *[CompletePolynomialGetNumMonomials(2,extrapOrder)];
        int PtrPtrIndex = 0;
        int PtrIndex = 0;
        for(int i=0;i<extrapOrder+1;++i){
            PascalsPyramid[i] = &DubPtrPtr[PtrPtrIndex];
            PtrPtrIndex += extrapOrder+1-i;
            for(int j=0;j<extrapOrder+1-i;++j){
                PascalsPyramid[i][j] = &GVec[PtrIndex];
                PtrIndex += j+1;
            }
        }

        //Create the first Pascal's Triangle that has no contribution from the third coordinate
        PascalsPyramid[0][0][0]=1.0;
        for(int i=1;i<extrapOrder+1;++i){
            PascalsPyramid[0][i][0]=coords[0]*PascalsPyramid[0][i-1][0];
            for(int j=1;j<i;++j){
                PascalsPyramid[0][i][j]=PascalsPyramid[0][i-1][j-1]*PascalsPyramid[0][i-1][j];
            }
            PascalsPyramid[0][i][i] = coords[1]*PascalsPyramid[0][i-1][i-1];
        }

        //Now create the rest of the Pascal's Triangles...
        for(int i=1;i<extrapOrder+1;++i){
            for(int j=0;j<extrapOrder+1-i;++j){
                for(int k=0;k<j+1;++k){
                    PascalsPyramid[i][j][k]=coords[2]*PascalsPyramid[i-1][j][k];
                }
            }
        }

        //Now loop through the completed Pascal's Pyramid and populate Gvec.
        for(int i=0;i<extrapOrder+1;++i){
            for(int j=0;j<extrapOrder+1-i;++j){
                for(int k=0;k<j+1;++k){
                    //std::cout << PascalsPyramid[i][j][k] << "  ";
                }
                //std::cout << std::endl;
            }
            //std::cout << std::endl;
        }

        delete [] PascalsPyramid;
        delete [] DubPtrPtr;
    }	
}
//============================================================================
int TensorProductGetNumMonomials(const int numDims,
                                 const int extrapOrder)
{
    //Returns the number of monomials for the given dimensionality and order.
    int numTerms = 1;
    for(int i=0;i<numDims;++i){
        numTerms *= (extrapOrder+1);
    }
    return numTerms;
}
//============================================================================
void  TensorProductMonomialVector(const int &numDims, 
                                  const double * coords, 
                                  const int &extrapOrder,
                                  std::vector<double> &GVec)
{
	//returns a vector of monomials in the tensor product up to the input extrapOrder evaluated at the input coordinate
    int NumMonomials = TensorProductGetNumMonomials(numDims,extrapOrder);
    GVec.resize(NumMonomials);

    //1D is very straightforward...
	if(numDims == 1){
		GVec[0]=1.0;
        //std::cout << GVec.at(0) << "\t";
		for(int i=1;i<extrapOrder+1;++i){
			GVec[i] = coords[0]*GVec[i-1];
            //std::cout << GVec.at(i) << "\t";
		}
	}

	if(numDims == 2){
        vector<double> Dim1Terms(extrapOrder+1,1.0);
        vector<double> Dim2Terms(extrapOrder+1,1.0);
        for(int i=1;i<extrapOrder+1;++i){
            Dim1Terms[i] = coords[0]*Dim1Terms[i-1];
            Dim2Terms[i] = coords[1]*Dim2Terms[i-1];
        }
        for(int i=0;i<extrapOrder+1;++i){
            for(int j=0;j<extrapOrder+1;++j){
                GVec[i*(extrapOrder+1) + j] = Dim1Terms[i] * Dim2Terms[j];
            }
        }
    }

    if(numDims == 3){
        vector<double> Dim1Terms(extrapOrder+1,1.0);
        vector<double> Dim2Terms(extrapOrder+1,1.0);
        vector<double> Dim3Terms(extrapOrder+1,1.0);
        for(int i=1;i<extrapOrder+1;++i){
            Dim1Terms[i] = coords[0]*Dim1Terms[i-1];
            Dim2Terms[i] = coords[1]*Dim2Terms[i-1];
            Dim3Terms[i] = coords[2]*Dim3Terms[i-1];
        }
        for(int i=0;i<extrapOrder+1;++i){
            for(int j=0;j<extrapOrder+1;++j){
                for(int k=0;k<extrapOrder+1;++k){
                    GVec[(i*(extrapOrder+1) + j)*(extrapOrder+1) + k] = Dim1Terms[i] * Dim2Terms[j] * Dim3Terms[k];
                }
            }
        }
    }	
}