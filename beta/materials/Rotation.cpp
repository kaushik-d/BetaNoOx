#include"stdafx.h"

#include<cmath>
#include<vector>
#include"utility/utility.h"
#include"Rotation.hpp"

///////////////////////////// Rotation Methods ////////////////////////////////

Rotation::Rotation(const int &CoordinateAxis, const double &RotationAngle):
b0(cos( M_PI/180 * RotationAngle/2.0)),
b1(0.0),
b2(0.0),
b3(0.0)
{
    switch(CoordinateAxis){
    case 1: // x axis rotation
        b1 = sin(M_PI/180 * RotationAngle/2.0);
        break;
    case 2: // y axis rotation
        b2 = sin(M_PI/180 * RotationAngle/2.0);
        break;
    case 3: // z axis rotation
        b3 = sin(M_PI/180 * RotationAngle/2.0);
        break;
    default:
        cout << "ERROR - Rotation::Rotation - Invalid value for CoordinateAxis" << endl;
        exit(1);
    }
}

///////////////////////////////////////////////////////////////////////////////
Rotation::Rotation(const double *Axis, const double &RotationAngle):
b0(cos( M_PI/180.0 * RotationAngle/2.0))
{
    //Assumes RotationAngle is size 3
    double magnitude = sqrt(Axis[0]*Axis[0] + Axis[1]*Axis[1] + Axis[2]*Axis[2]);
    //Catch to prevent divide by zero
    if(magnitude == 0){
        cout << "ERROR - Rotation::Rotation - Input rotation axis has zero length" << endl;
        exit(1);
    }
    b1 = Axis[0]/magnitude * sin( M_PI/180.0 * RotationAngle/2.0);
    b2 = Axis[1]/magnitude * sin( M_PI/180.0 * RotationAngle/2.0);
    b3 = Axis[2]/magnitude * sin( M_PI/180.0 * RotationAngle/2.0);
}

///////////////////////////////////////////////////////////////////////////////
Rotation::Rotation(const InterpAndSolDerivatives* Interp, const Rotation* &NodalRotations):
b0(0.0),b1(0.0),b2(0.0),b3(0.0)
{
    interpolate_rotation(Interp,NodalRotations);
}

void Rotation::interpolate_rotation(const InterpAndSolDerivatives* Interp, const Rotation* NodalRotations)
{
    // This method does not give the same result as the old approach used in Beta
    // since in essence you are interpolating cos(theta/2) instead of theta itself.
    // The result of this interpolation is somewhat less accurate in that regard.

    // This method of interpolation _is_ frame indifferent (applying an additional
    // rotation to the interpolated rotation yields the same result as applying the 
    // additional rotation to all the nodal rotations and then interpolating those)

    // If some really smart graduate student comes along some day, they could
    // probably read up on interpolation of quaternions and improve this method
    // to avoid the creation of small errors associated with using Galerkin
    // interpolation on quaternions and then re-normalizing.  There are a few
    // papers on doing this for 1D elements, but at this point I haven't seen
    // any regarding fully 3D finite elements.
    b0=b1=b2=b3=0.0;

    for(int i=0;i<Interp->numberOfInterp;++i){
        b0 += Interp->S[i] * NodalRotations[i].getb0();
        b1 += Interp->S[i] * NodalRotations[i].getb1();
        b2 += Interp->S[i] * NodalRotations[i].getb2();
        b3 += Interp->S[i] * NodalRotations[i].getb3();
    }
    normalize_self();
}

///////////////////////////////////////////////////////////////////////////////
void Rotation::interpolate_old_beta(const InterpAndSolDerivatives* Interp, 
                                    const Rotation* NodalRotations,
                                    const Rotation &MaterialRotation)
{
    // This method is designed to give the same result as the old Beta when
    // only the tow undulation rotations are being interpolated.
    // !!!ONLY WORKS FOR UNDULATION+MATERIAL COMPOSITE ROTATION!!!
    // First, the material rotation is removed from the overall rotation.
    // Then, it uses a quaternion interpolation to calculate the axis of rotation,
    // then it interpolates the angle itself instead of interpolating
    // cos(theta/2) (b0 from the quaternion).  Then the material rotation is added
    // back onto the end of the rotation sequence.  The result is more accurate
    // but obviously more expensive as well.
    b0=b1=b2=b3=0.0;

    for(int i=0;i<Interp->numberOfInterp;++i){
        Rotation RotWithoutMatRotation = NodalRotations[i] - MaterialRotation;
        b0 += Interp->S[i] * RotWithoutMatRotation.getb0();
        b1 += Interp->S[i] * RotWithoutMatRotation.getb1();
        b2 += Interp->S[i] * RotWithoutMatRotation.getb2();
        b3 += Interp->S[i] * RotWithoutMatRotation.getb3();
    }
    normalize_self();

    // Now get the axis for the rotation
    // This has to be done this way because if the rotation for some node 
    // is zero, then the axis of rotation for that node is undefined.
    double SineAngleBy2 = sin(acos(b0));
    double axis[] = {1.0, 0.0, 0.0};
    if(!(SineAngleBy2 < 1.0e-10)){
        axis[0] = b1/SineAngleBy2;
        axis[1] = b2/SineAngleBy2;
        axis[2] = b3/SineAngleBy2;
    }

    //Now interpolate the angle iteslf using Galerkin interpolation
    double angle = 0.0;
    for(int i=0;i<Interp->numberOfInterp;++i){
        // Note, since only one of the parameters is needed here, efficiency could be improved by
        // explicitly coding the subtraction of the material rotation (instead of using the 
        // minus operator which calculates all 4 parameters)
        angle += Interp->S[i] * 2*acos((NodalRotations[i] - MaterialRotation).getb0());
    }
    //Now reconstruct the quaternion for this interpolated angle...
    SineAngleBy2 = sin(angle/2);
    b0 = cos(angle/2.0);
    b1 = axis[0]*SineAngleBy2;
    b2 = axis[1]*SineAngleBy2;
    b3 = axis[2]*SineAngleBy2;

    //Add the material rotation back onto the end.
    add_after_self(MaterialRotation);
}

///////////////////////////////////////////////////////////////////////////////
void Rotation::interpolate_rotation_by_angle(const InterpAndSolDerivatives* Interp, 
                                             const Rotation* NodalRotations)
{
    // This is a more accurate but _much_ more expensive way to interpolate the 
    // rotation at a quadrature point.  It is based on the fact that for most
    // cases in textiles, the material rotations in an element vary by only a 
    // single rotation about a single, principle axis.  The angle of this single 
    // rotation is then interpolated, rather than simply interpolating and 
    // re-normalizing the quaternion parameters.

    // Determine if there is a common axis among all the orientations.
    bool XAxisConstant = true;
    bool YAxisConstant = true;
    bool ZAxisConstant = true;
    // Each local axis expressed in the global axes forms the columns of the rotation's DCM
    Matrix FirstRotationDCM(3,3), OtherRotationDCM(3,3);
    NodalRotations[0].get_direction_cosine_matrix(FirstRotationDCM);
    for(int i=1;i<Interp->numberOfInterp;++i){
        NodalRotations[i].get_direction_cosine_matrix(OtherRotationDCM);
        for(int j=0;j<3;++j){
            if(fabs(OtherRotationDCM(j,0)-FirstRotationDCM(j,0)) > 1.0e-8) XAxisConstant = false;
            if(fabs(OtherRotationDCM(j,1)-FirstRotationDCM(j,1)) > 1.0e-8) YAxisConstant = false;
            if(fabs(OtherRotationDCM(j,2)-FirstRotationDCM(j,2)) > 1.0e-8) ZAxisConstant = false;
        }
    }
    if(XAxisConstant && YAxisConstant && ZAxisConstant){
        // All rotations are the same.
        *this = NodalRotations[0];
    }else if(XAxisConstant || YAxisConstant || ZAxisConstant){
        // All the rotated coordinate systems share one common axis.

        // Determine the rotation angle at the quadrature point relative to first rotation
        // Relative angle in radians can be determined directly from the quaternions...
        // The math depends on the axis that the rotation is about
        // Basically, you are calculating the b1, b2, or b3 component of the quaternion that will
        // rotate the local axes at node zero to the local axes at node i, and then taking
        // 2 * arcsin of that component.  This will give you the angle of rotation with the 
        // correct sign.  Look at Eq 3.98 in Schaub & Junkins (P.98).  Also recall that the
        // rotation class stores the rotation that takes the axes from the local coordinate
        // system to the global coordinate system.  Therefore, the LHS (unprimed beta) is the
        // inverse rotation at the ith node of the element.  The the single-primed beta is the
        // inverse rotation at node zero.  The result of rotating the primed beta by the
        // double-primed beta is the LHS (unprimed beta).  Therefore, we want to determine this
        // double primed from the single-primed and unprimed beta to calculate the rotation angle
        // between the two orientations.  Note that the double-primed beta will only have a single
        // component b1, b2, b3 that is non-zero since it represents a rotation about a principle
        // axis.

        // Note that a simpler but  more computationally expensive way to obtain the rotation
        // between NodalRotations[0] and NodalRotations[i] would be to take 
        // NodalRotations[0] - NodalRotations[i].  This would unnecessarily calculate all
        // values of the resulting rotation, most of which would be zero
        vector<double> AngleRadians(Interp->numberOfInterp,0.0);
        Rotation ZeroNodeInverseRotation(NodalRotations[0].get_reverse());
        Rotation iNodeInverseRotation;
        if(XAxisConstant){
            for(int i=1;i<Interp->numberOfInterp;++i){
                iNodeInverseRotation = NodalRotations[i].get_reverse();
                AngleRadians[i] = 2*asin(-ZeroNodeInverseRotation.b1 * iNodeInverseRotation.b0+
                                          ZeroNodeInverseRotation.b0 * iNodeInverseRotation.b1+
                                          ZeroNodeInverseRotation.b3 * iNodeInverseRotation.b2-
                                          ZeroNodeInverseRotation.b2 * iNodeInverseRotation.b3);
            }
        }else if(YAxisConstant){
            for(int i=1;i<Interp->numberOfInterp;++i){
                iNodeInverseRotation = NodalRotations[i].get_reverse();
                AngleRadians[i] = 2*asin(-ZeroNodeInverseRotation.b2 * iNodeInverseRotation.b0-
                                          ZeroNodeInverseRotation.b3 * iNodeInverseRotation.b1+
                                          ZeroNodeInverseRotation.b0 * iNodeInverseRotation.b2+
                                          ZeroNodeInverseRotation.b1 * iNodeInverseRotation.b3);
            }
        }else if(ZAxisConstant){
            for(int i=1;i<Interp->numberOfInterp;++i){
                iNodeInverseRotation = NodalRotations[i].get_reverse();
                AngleRadians[i] = 2*asin(-ZeroNodeInverseRotation.b3 * iNodeInverseRotation.b0+
                                          ZeroNodeInverseRotation.b2 * iNodeInverseRotation.b1-
                                          ZeroNodeInverseRotation.b1 * iNodeInverseRotation.b2+
                                          ZeroNodeInverseRotation.b0 * iNodeInverseRotation.b3);
            }
        }
        //Calculate the angle at the quadrature point
        double QuadPointAngleRadians = 0.0;
        for(int i=0;i<Interp->numberOfInterp;++i){
            QuadPointAngleRadians += AngleRadians[i]*Interp->S[i];
        }
        //Create the quaternion at the quad point by applying an incremental rotation to the
        //orientation at the first node.
        *this = NodalRotations[0].get_reverse();
        if(XAxisConstant){
            add_after_self(Rotation(1,QuadPointAngleRadians*180.0/M_PI));
        }else if(YAxisConstant){
            add_after_self(Rotation(2,QuadPointAngleRadians*180.0/M_PI));
        }else if(ZAxisConstant){
            add_after_self(Rotation(3,QuadPointAngleRadians*180.0/M_PI));
        }
        // The rotation was calculated from global to local, but it needs to be expressed
        // in terms of going from local to global (the inverse)
        *this = get_reverse();
    }else{
        // If there is no common axis between the nodes.
        // It might be a good idea to print a warning here.  Won't for now.
        interpolate_rotation(Interp,NodalRotations);
        return;
    }
}

///////////////////////////////////////////////////////////////////////////////
void Rotation::add_after_self(const Rotation &Arg)
{
    // This is a matrix-vector multiply given on P98 of Schaub & Junkins
    double newb0, newb1, newb2, newb3;
    newb0 = Arg.b0*b0 - Arg.b1*b1 - Arg.b2*b2 - Arg.b3*b3;
    newb1 = Arg.b1*b0 + Arg.b0*b1 + Arg.b3*b2 - Arg.b2*b3;
    newb2 = Arg.b2*b0 - Arg.b3*b1 + Arg.b0*b2 + Arg.b1*b3;
    newb3 = Arg.b3*b0 + Arg.b2*b1 - Arg.b1*b2 + Arg.b0*b3;

    b0 = newb0;
    b1 = newb1;
    b2 = newb2;
    b3 = newb3;
}

///////////////////////////////////////////////////////////////////////////////
void Rotation::add_before_self(const Rotation &Arg)
{
    double newb0, newb1, newb2, newb3;
    newb0 = b0*Arg.b0 - b1*Arg.b1 - b2*Arg.b2 - b3*Arg.b3;
    newb1 = b1*Arg.b0 + b0*Arg.b1 + b3*Arg.b2 - b2*Arg.b3;
    newb2 = b2*Arg.b0 - b3*Arg.b1 + b0*Arg.b2 + b1*Arg.b3;
    newb3 = b3*Arg.b0 + b2*Arg.b1 - b1*Arg.b2 + b0*Arg.b3;

    b0 = newb0;
    b1 = newb1;
    b2 = newb2;
    b3 = newb3;
}

///////////////////////////////////////////////////////////////////////////////
Rotation Rotation::operator + (const Rotation &Arg) const
{
    Rotation NewRotation(*this);
    NewRotation += Arg;
    return NewRotation;
}

///////////////////////////////////////////////////////////////////////////////
void Rotation::get_Euler_axis_and_angle(vector<double> &EulerAxis, double &Angle) const
{
    double AngleRad = 2*acos(b0);
    Angle = AngleRad * 180.0/M_PI;
    EulerAxis.resize(3);
    if(Angle <1.0e-8){
        EulerAxis[0] = 1.0;
        EulerAxis[1] = 0.0;
        EulerAxis[2] = 0.0;
    }else{
        EulerAxis[0] = b1/sin(AngleRad/2);
        EulerAxis[1] = b2/sin(AngleRad/2);
        EulerAxis[2] = b3/sin(AngleRad/2);
    }
}

///////////////////////////////////////////////////////////////////////////////
void Rotation::output_rotation(ostream* OS) const
{
    vector<double> EulerAxis(3,0.0);
    double Angle;
    get_Euler_axis_and_angle(EulerAxis,Angle);
    *OS << " Axis = <" << EulerAxis[0] << ", " << EulerAxis[1] << ", " << EulerAxis[2] << ">" << "\tAngle = " << Angle;
}

///////////////////////////////////////////////////////////////////////////////
void Rotation::get_direction_cosine_matrix(Matrix &DCM) const
{
    //Time is important here, so assume that DCM is already 3x3 (don't check it).
    //If it isn't, you'll probably get a seg fault or corrupt some
    //other memory in the problem and end up with a really nasty crash.

    //could possibly put in a check to see if b0 == 1.0, in which case rotation
    //angle is zero, in which case DCM should be identity

    //See Schaub and Junkins, P.96, Eq 3.92
    //Perform all multiplications first
    double b0sq = b0*b0;
    double b1sq = b1*b1;
    double b2sq = b2*b2;
    double b3sq = b3*b3;
    double b0b1 = b0*b1;
    double b0b2 = b0*b2;
    double b0b3 = b0*b3;
    double b1b2 = b1*b2;
    double b1b3 = b1*b3;
    double b2b3 = b2*b3;

    //Assemble DCM with the above terms
    DCM(0,0) = b0sq + b1sq - b2sq - b3sq;
    DCM(0,1) = 2*(b1b2 + b0b3);
    DCM(0,2) = 2*(b1b3 - b0b2);

    DCM(1,0) = 2*(b1b2 - b0b3);
    DCM(1,1) = b0sq - b1sq + b2sq - b3sq;
    DCM(1,2) = 2*(b2b3 + b0b1);

    DCM(2,0) = 2*(b1b3 + b0b2);
    DCM(2,1) = 2*(b2b3 - b0b1);
    DCM(2,2) = b0sq - b1sq - b2sq + b3sq;
}

///////////////////////////////////////////////////////////////////////////////
Matrix Rotation::return_direction_cosine_matrix() const
{
    Matrix DCM(3,3);
    get_direction_cosine_matrix(DCM);
    return DCM;
}

///////////////////////////////////////////////////////////////////////////////
void Rotation::get_Voigt_rotation_matrix(Matrix &VRM) const
{
    //See rotation_transformation/VoigtTransformation.mws for this derivation.
    //Assumes that VRM is 6x6, will probably crash if it isn't

    Matrix a(3,3); //The Direction Cosine Matrix

    //Get the DCM
    //Note that there doesn't currently appear to be no performance benefit to be had from
    //directly using the quaternions in the computation of the Voigt rotation matrix, but
    //the checks made so far haven't been extremely thorough.  Perhaps some speedup could
    //be obtained with such an approach and using temporary intermediate variables
    get_direction_cosine_matrix(a);

    VRM(0,0) = a(0,0) * a(0,0);
    VRM(0,1) = a(0,1) * a(0,1);
    VRM(0,2) = a(0,2) * a(0,2);
    VRM(0,3) = 2.0 * a(0,1) * a(0,0);
    VRM(0,4) = 2.0 * a(0,2) * a(0,1);
    VRM(0,5) = 2.0 * a(0,2) * a(0,0);

    VRM(1,0) = a(1,0) * a(1,0);
    VRM(1,1) = a(1,1) * a(1,1);
    VRM(1,2) = a(1,2) * a(1,2);
    VRM(1,3) = 2.0 * a(1,1) * a(1,0);
    VRM(1,4) = 2.0 * a(1,2) * a(1,1);
    VRM(1,5) = 2.0 * a(1,2) * a(1,0);

    VRM(2,0) = a(2,0)*a(2,0);
    VRM(2,1) = a(2,1)*a(2,1);
    VRM(2,2) = a(2,2)*a(2,2);
    VRM(2,3) = 2.0 * a(2,1) * a(2,0);
    VRM(2,4) = 2.0 * a(2,2) * a(2,1);
    VRM(2,5) = 2.0 * a(2,2) * a(2,0);

    VRM(3,0) = a(0,0) * a(1,0);
    VRM(3,1) = a(0,1) * a(1,1);
    VRM(3,2) = a(0,2) * a(1,2);
    VRM(3,3) = a(0,1) * a(1,0) + a(0,0) * a(1,1);
    VRM(3,4) = a(0,2) * a(1,1) + a(0,1) * a(1,2);
    VRM(3,5) = a(0,2) * a(1,0) + a(0,0) * a(1,2);

    VRM(4,0) = a(1,0) * a(2,0);
    VRM(4,1) = a(1,1) * a(2,1);
    VRM(4,2) = a(1,2) * a(2,2);
    VRM(4,3) = a(1,1) * a(2,0) + a(1,0) * a(2,1);
    VRM(4,4) = a(1,2) * a(2,1) + a(1,1) * a(2,2);
    VRM(4,5) = a(1,2) * a(2,0) + a(1,0) * a(2,2);

    VRM(5,0) = a(0,0) * a(2,0);
    VRM(5,1) = a(0,1) * a(2,1);
    VRM(5,2) = a(0,2) * a(2,2);
    VRM(5,3) = a(0,1) * a(2,0) + a(0,0) * a(2,1);
    VRM(5,4) = a(0,2) * a(2,1) + a(0,1) * a(2,2);
    VRM(5,5) = a(0,2) * a(2,0) + a(0,0) * a(2,2);    
}

///////////////////////////////////////////////////////////////////////////////
void Rotation::get_Voigt_engineering_rotation_matrix(Matrix &VESRM) const
{
    //See rotation_transformation/VoigtTransformation.mws for this derivation.
    //Assumes that VESRM is 6x6, will probably crash if it isn't

    Matrix a(3,3); //The Direction Cosine Matrix

    //Get the DCM
    //Note that there doesn't currently appear to be no performance benefit to be had from
    //directly using the quaternions in the computation of the Voigt rotation matrix, but
    //the checks made so far haven't been extremely thorough.  Perhaps some speedup could
    //be obtained with such an approach and using temporary intermediate variables
    get_direction_cosine_matrix(a);

    VESRM(0,0) = a(0,0) * a(0,0);
    VESRM(0,1) = a(0,1) * a(0,1);
    VESRM(0,2) = a(0,2) * a(0,2);
    VESRM(0,3) = a(0,1) * a(0,0);
    VESRM(0,4) = a(0,2) * a(0,1);
    VESRM(0,5) = a(0,2) * a(0,0);

    VESRM(1,0) = a(1,0) * a(1,0);
    VESRM(1,1) = a(1,1) * a(1,1);
    VESRM(1,2) = a(1,2) * a(1,2);
    VESRM(1,3) = a(1,1) * a(1,0);
    VESRM(1,4) = a(1,2) * a(1,1);
    VESRM(1,5) = a(1,2) * a(1,0);

    VESRM(2,0) = a(2,0) * a(2,0);
    VESRM(2,1) = a(2,1) * a(2,1);
    VESRM(2,2) = a(2,2) * a(2,2);
    VESRM(2,3) = a(2,1) * a(2,0);
    VESRM(2,4) = a(2,2) * a(2,1);
    VESRM(2,5) = a(2,2) * a(2,0);

    VESRM(3,0) = 2.0 * a(0,0) * a(1,0);
    VESRM(3,1) = 2.0 * a(0,1) * a(1,1);
    VESRM(3,2) = 2.0 * a(0,2) * a(1,2);
    VESRM(3,3) = a(0,1) * a(1,0) + a(0,0) * a(1,1);
    VESRM(3,4) = a(0,2) * a(1,1) + a(0,1) * a(1,2);
    VESRM(3,5) = a(0,2) * a(1,0) + a(0,0) * a(1,2);

    VESRM(4,0) = 2.0 * a(1,0) * a(2,0);
    VESRM(4,1) = 2.0 * a(1,1) * a(2,1);
    VESRM(4,2) = 2.0 * a(1,2) * a(2,2);
    VESRM(4,3) = a(1,1) * a(2,0) + a(1,0) * a(2,1);
    VESRM(4,4) = a(1,2) * a(2,1) + a(1,1) * a(2,2);
    VESRM(4,5) = a(1,2) * a(2,0) + a(1,0) * a(2,2);

    VESRM(5,0) = 2.0 * a(0,0) * a(2,0);
    VESRM(5,1) = 2.0 * a(0,1) * a(2,1);
    VESRM(5,2) = 2.0 * a(0,2) * a(2,2);
    VESRM(5,3) = a(0,1) * a(2,0) + a(0,0) * a(2,1);
    VESRM(5,4) = a(0,2) * a(2,1) + a(0,1) * a(2,2);
    VESRM(5,5) = a(0,2) * a(2,0) + a(0,0) * a(2,2);

}

///////////////////////////////////////////////////////////////////////////////
void Rotation::express_3D_2nd_Order_Tensor_in_primed_CSYS(const Matrix &UnprimedTensor,
                                                                Matrix &PrimedTensor) const
{
    if(no_rotation()){
        PrimedTensor = UnprimedTensor;
        return;
    }

    Matrix a(3,3);
    get_direction_cosine_matrix(a);

    PrimedTensor = (a * UnprimedTensor) * a.transposeCopy();

    //double *UnprimedTensorRow,*aTransposedRow,*aRow;
    //int i,j,k,l;
    //{
    //    register double sum;
    //    for(i=0;i<3;++i) {
    //        aTransposedRow = a[i];
    //        for(j=0;j<3;++j) {
    //            aRow = a[j];
    //            sum = 0;
    //            for(k=0;k<3;++k) {
    //                UnprimedTensorRow = UnprimedTensor[k];
    //                for(l=0;l<3;++l) {
    //                    sum += aTransposedRow[k] * UnprimedTensorRow[l] * aRow[l];
    //                }
    //            }
    //            PrimedTensor(i,j) = sum;
    //        }
    //    }
    //}
}

///////////////////////////////////////////////////////////////////////////////
void Rotation::express_3D_1st_Order_Tensor_in_primed_CSYS(const double* unprimedVec, double* primedVec)const
{
    if(no_rotation()){
        copyVector(unprimedVec,primedVec,3);
        return;
    }

    Matrix a(3,3);
    get_direction_cosine_matrix(a);
    int i,j;
    primedVec[0] = primedVec[1] = primedVec[2] = 0.0;
    for(i=0;i<3;++i){
        for(j=0;j<3;++j){
            primedVec[i] += a(i,j) * unprimedVec[j];
        }
    }
}

///////////////////////////////////////////////////////////////////////////////
void Rotation::express_Cmat_in_primed_CSYS(const Matrix &unprimedCMat, Matrix &primedCMat) const
{
    // bypass method if there is no rotation
    if(no_rotation()){
        primedCMat = unprimedCMat;
        return;
    }

    //Special Logic if the C matrix is for a 2D material
    if(unprimedCMat.Size() == 3){
        cout << "Logic for rotating the C Matrix of a 2D" << endl;
        cout << "elastic solid is currently not implemented." << endl;
        exit(1);
    }

    Matrix Q(6,6);
    get_Voigt_rotation_matrix(Q);

    //C'_ij = Q_ik Q_jl C_kl or C' = Q C Q^T

    //This definitely does the operation correctly, but will be slower because of copies that take place
    //GlobalCMat = (Q * unprimedCMat) * Q.transposeCopy();

    // Note that one could potentially exploit the fact that for orthotropic materials, the unprimed Cmat
    // will contain a large number of zeros to reduce the number of computations required here

    double *unprimedCRow,*QTransposedRow,*QRow;
    int i,j,k,l;
    {
        register double sum;
        for(i=0;i<6;++i) {
            QTransposedRow = Q[i];
            for(j=0;j<6;++j) {
                QRow = Q[j];
                sum = 0;
                for(k=0;k<6;++k) {
                    unprimedCRow = unprimedCMat[k];
                    for(l=0;l<6;++l) {
                        sum += QTransposedRow[k] * unprimedCRow[l] * QRow[l];
                    }
                }
                primedCMat(i,j) = sum;
            }
        }
    }

}

///////////////////////////////////////////////////////////////////////////////
void Rotation::express_Voigt_stress_in_primed_CSYS(const double* unprimedStress,
                                                         double* primedStress) const
{
    // bypass method if there is no rotation
    if(no_rotation()){
        copyVector(unprimedStress,primedStress,6);
        return;
    }

    Matrix Q(6,6);
    get_Voigt_rotation_matrix(Q);

    //sigma'_i = Q_ij sigma_j
    int i,j;
    for(i=0;i<6;++i) {
        primedStress[i] = 0;
        for(j=0;j<6;++j) {
            primedStress[i] += Q(i,j) * unprimedStress[j];
        }
    }
}

///////////////////////////////////////////////////////////////////////////////
void Rotation::express_Voigt_engineering_strain_in_primed_CSYS(const double* unprimedStrain,
                                                                     double* primedStrain) const
{
    // bypass method if there is no rotation
    if(no_rotation()){
        copyVector(unprimedStrain,primedStrain,6);
        return;
    }

    Matrix Q(6,6);
    get_Voigt_engineering_rotation_matrix(Q);

    //sigma'_i = Q_ij sigma_j
    int i,j;
    for(i=0;i<6;++i) {
        primedStrain[i] = 0;
        for(j=0;j<6;++j) {
            primedStrain[i] += Q(i,j) * unprimedStrain[j];
        }
    }
}

///////////////////////////////////////////////////////////////////////////////
void Rotation::express_2D_stress_in_primed_CSYS(const double* unprimedStrain, 
                                                      double* primedStrain) const
{
    double Angle = 2*acos(b0);
    //If the rotation angle is nonzero and about some axis other than z, complain
    if( ( (b1 > 1.0e-8) || (b2 > 1.0e-8) ) && (Angle > 1.0e-8)){
        cout << "You are trying to rotate a 2D material about an axis other than z!" << endl;
        cout << "This does not make sense" << endl;
        exit(1);
    }
    double t1 = (unprimedStrain[0]+unprimedStrain[1])/2;
    double t2 = (unprimedStrain[0]-unprimedStrain[1])/2*cos(2*Angle);
    double t3 = unprimedStrain[2]*sin(2*Angle);
    primedStrain[0] = t1+t2+t3;
    primedStrain[1] = t1-t2-t3;
    primedStrain[2] = 2 * ((unprimedStrain[0]-unprimedStrain[1])/2*sin(2*Angle) + unprimedStrain[2]*cos(2*Angle));
}

///////////////////////////////////////////////////////////////////////////////
void Rotation::express_2D_engineering_strain_in_primed_CSYS(const double* unprimedStrain, 
                                                                  double* primedStrain) const
{
    double Angle = 2*acos(b0);
    //If the rotation angle is nonzero and about some axis other than z, complain
    if( ( (b1 > 1.0e-8) || (b2 > 1.0e-8) ) && (Angle > 1.0e-8)){
        cout << "You are trying to rotate a 2D material about an axis other than z!" << endl;
        cout << "This does not make sense" << endl;
        exit(1);
    }
    double tensorShear = unprimedStrain[2] / 2.0;
    double t1 = (unprimedStrain[0]+unprimedStrain[1])/2;
    double t2 = (unprimedStrain[0]-unprimedStrain[1])/2*cos(2*Angle);
    double t3 = tensorShear*sin(2*Angle);
    primedStrain[0] = t1+t2+t3;
    primedStrain[1] = t1-t2-t3;
    primedStrain[2] = 2 * ((unprimedStrain[0]-unprimedStrain[1])/2*sin(2*Angle) + tensorShear*cos(2*Angle));
}