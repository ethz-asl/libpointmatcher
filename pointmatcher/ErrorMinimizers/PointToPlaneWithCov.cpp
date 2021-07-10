// kate: replace-tabs off; indent-width 4; indent-mode normal
// vim: ts=4:sw=4:noexpandtab
/*

Copyright (c) 2010--2012,
François Pomerleau and Stephane Magnenat, ASL, ETHZ, Switzerland
You can contact the authors at <f dot pomerleau at gmail dot com> and
<stephane at magnenat dot net>

All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:
    * Redistributions of source code must retain the above copyright
      notice, this list of conditions and the following disclaimer.
    * Redistributions in binary form must reproduce the above copyright
      notice, this list of conditions and the following disclaimer in the
      documentation and/or other materials provided with the distribution.
    * Neither the name of the <organization> nor the
      names of its contributors may be used to endorse or promote products
      derived from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL ETH-ASL BE LIABLE FOR ANY
DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

*/

#include <iostream>

#include "Eigen/SVD"

#include "ErrorMinimizersImpl.h"
#include "PointMatcherPrivate.h"
#include "Functions.h"

using namespace Eigen;
using namespace std;

typedef PointMatcherSupport::Parametrizable Parametrizable;
typedef PointMatcherSupport::Parametrizable P;
typedef Parametrizable::Parameters Parameters;
typedef Parametrizable::ParameterDoc ParameterDoc;
typedef Parametrizable::ParametersDoc ParametersDoc;

template<typename T>
PointToPlaneWithCovErrorMinimizer<T>::PointToPlaneWithCovErrorMinimizer(const Parameters& params):
	PointToPlaneErrorMinimizer<T>(PointToPlaneWithCovErrorMinimizer::availableParameters(), params),
        sensorStdDev(Parametrizable::get<T>("sensorStdDev"))
{
}

template<typename T>
typename PointMatcher<T>::TransformationParameters PointToPlaneWithCovErrorMinimizer<T>::compute(const ErrorElements& mPts_const)
{
	ErrorElements mPts = mPts_const;
    typename PointMatcher<T>::TransformationParameters out = PointToPlaneErrorMinimizer<T>::compute_in_place(mPts);

    if(!PointToPlaneErrorMinimizer<T>::force2D) {
        this->covMatrix = this->estimateCovariance3d(mPts, out);
    } else {
        this->covMatrix = this->estimateCovariance2d(mPts, out);
    }

    return out;
}

template<typename T>
typename PointMatcher<T>::Matrix
PointToPlaneWithCovErrorMinimizer<T>::estimateCovariance3d(const ErrorElements& mPts, const TransformationParameters& transformation)
{

    // Create our empty Jacobian matices
    const int max_nbr_point = mPts.reading.getNbPoints();
    Matrix d2J_dX2(Matrix::Zero(6,6));
    Matrix d2J_dZdX(Matrix::Zero(6,6*max_nbr_point));

    // Normals of the reference cloud (taken to be true)
    Matrix normals = mPts.reference.getDescriptorViewByName("normals");

    // Return if we do not have any normals in DataPoints
    if (normals.rows() < 3)
        return std::numeric_limits<T>::max() * Matrix::Identity(6,6);

    T beta = -asin(transformation(2,0));
    T alpha = atan2(transformation(2,1), transformation(2,2));
    T gamma = atan2(transformation(1,0), transformation(0,0));
    T t_x = transformation(0,3);
    T t_y = transformation(1,3);
    T t_z = transformation(2,3);

    // According to the rotation matrix I used and after verification
    // It is Yaw Pitch ROLL = [a,b,c]== [R] matrix used in the MatLab
    double x, y, z, a, b, c;
    x = t_x; y = t_y; z = t_z;
    a = gamma; b = beta; c = alpha;

    // Loop through all points and calculate the state and partial Jacobians
    for(int i = 0; i < max_nbr_point; ++i)
    {

        double pix = mPts.reading.features(0,i);
        double piy = mPts.reading.features(1,i);
        double piz = mPts.reading.features(2,i);
        double qix = mPts.reference.features(0,i);
        double qiy = mPts.reference.features(1,i);
        double qiz = mPts.reference.features(2,i);
        double nix = normals(0,i);
        double niy = normals(1,i);
        double niz = normals(2,i);

        /***********************************************************************
        d2J_dX2 -- X is the [R|T] in the form of [x, y, z, a, b, c]
        x, y, z is the translation part
        a, b, c is the rotation part in Euler format
        [x, y, z, a, b, c] is acquired from the Transformation Matrix returned by ICP.
        Now d2J_dX2 is a 6x6 matrix of the form
        d2J_dx2
        d2J_dxdy    d2J_dy2
        d2J_dxdz    d2J_dydz    d2J_dz2
        d2J_dxda    d2J_dyda    d2J_dzda   d2J_da2
        d2J_dxdb    d2J_dydb    d2J_dzdb   d2J_dadb   d2J_db2
        d2J_dxdc    d2J_dydc    d2J_dzdc   d2J_dadc   d2J_dbdc   d2J_dc2
        *************************************************************************/

        double 	d2J_dx2,     d2J_dydx,	  d2J_dzdx,   d2J_dadx,   d2J_dbdx,     d2J_dcdx,
                d2J_dxdy,    d2J_dy2,	  d2J_dzdy,   d2J_dady,   d2J_dbdy,     d2J_dcdy,
                d2J_dxdz,    d2J_dydz,    d2J_dz2,    d2J_dadz,   d2J_dbdz,     d2J_dcdz,
                d2J_dxda,    d2J_dyda,    d2J_dzda,   d2J_da2,	  d2J_dbda,     d2J_dcda,
                d2J_dxdb,    d2J_dydb,    d2J_dzdb,   d2J_dadb,   d2J_db2,      d2J_dcdb,
                d2J_dxdc,    d2J_dydc,    d2J_dzdc,   d2J_dadc,   d2J_dbdc,     d2J_dc2;

        // These terms are generated from the provided Matlab scipts. We just have to copy
        // the expressions from the matlab output with two very simple changes.
        // The first one being the the sqaure of a number 'a' is shown as a^2 in matlab,
        // which is converted to pow(a,2) in the below expressions.
        // The second change is to add ';' at the end of each expression :)
        // In this way, matlab can be used to generate these terms for various objective functions of ICP
        // and they can simply be copied to the C++ files and with appropriate changes to ICP estimation,
        // its covariance can be easily estimated.

        d2J_dx2 = 2*pow(nix,2);
        d2J_dy2 = 2*pow(niy,2);
        d2J_dz2 = 2*pow(niz,2);
        d2J_dydx =2*nix*niy;
        d2J_dxdy = 2*nix*niy;
        d2J_dzdx = 2*nix*niz;
        d2J_dxdz = 2*nix*niz;
        d2J_dydz = 2*niy*niz;
        d2J_dzdy = 2*niy*niz;
        d2J_da2 = (niy*(piz*(sin(a)*sin(c) + cos(a)*cos(c)*sin(b)) - piy*(cos(c)*sin(a) - cos(a)*sin(b)*sin(c)) + pix*cos(a)*cos(b)) - nix*(piy*(cos(a)*cos(c) + sin(a)*sin(b)*sin(c)) - piz*(cos(a)*sin(c) - cos(c)*sin(a)*sin(b)) + pix*cos(b)*sin(a)))*(2*niy*(piz*(sin(a)*sin(c) + cos(a)*cos(c)*sin(b)) - piy*(cos(c)*sin(a) - cos(a)*sin(b)*sin(c)) + pix*cos(a)*cos(b)) - 2*nix*(piy*(cos(a)*cos(c) + sin(a)*sin(b)*sin(c)) - piz*(cos(a)*sin(c) - cos(c)*sin(a)*sin(b)) + pix*cos(b)*sin(a))) - (2*nix*(piz*(sin(a)*sin(c) + cos(a)*cos(c)*sin(b)) - piy*(cos(c)*sin(a) - cos(a)*sin(b)*sin(c)) + pix*cos(a)*cos(b)) + 2*niy*(piy*(cos(a)*cos(c) + sin(a)*sin(b)*sin(c)) - piz*(cos(a)*sin(c) - cos(c)*sin(a)*sin(b)) + pix*cos(b)*sin(a)))*(nix*(x - qix - piy*(cos(c)*sin(a) - cos(a)*sin(b)*sin(c)) + piz*(sin(a)*sin(c) + cos(a)*cos(c)*sin(b)) + pix*cos(a)*cos(b)) + niy*(y - qiy + piy*(cos(a)*cos(c) + sin(a)*sin(b)*sin(c)) - piz*(cos(a)*sin(c) - cos(c)*sin(a)*sin(b)) + pix*cos(b)*sin(a)) + niz*(z - qiz - pix*sin(b) + piz*cos(b)*cos(c) + piy*cos(b)*sin(c)));
        d2J_db2 = (niy*(piz*cos(b)*cos(c)*sin(a) - pix*sin(a)*sin(b) + piy*cos(b)*sin(a)*sin(c)) - niz*(pix*cos(b) + piz*cos(c)*sin(b) + piy*sin(b)*sin(c)) + nix*(piz*cos(a)*cos(b)*cos(c) - pix*cos(a)*sin(b) + piy*cos(a)*cos(b)*sin(c)))*(2*niy*(piz*cos(b)*cos(c)*sin(a) - pix*sin(a)*sin(b) + piy*cos(b)*sin(a)*sin(c)) - 2*niz*(pix*cos(b) + piz*cos(c)*sin(b) + piy*sin(b)*sin(c)) + 2*nix*(piz*cos(a)*cos(b)*cos(c) - pix*cos(a)*sin(b) + piy*cos(a)*cos(b)*sin(c))) - (2*niy*(pix*cos(b)*sin(a) + piz*cos(c)*sin(a)*sin(b) + piy*sin(a)*sin(b)*sin(c)) + 2*niz*(piz*cos(b)*cos(c) - pix*sin(b) + piy*cos(b)*sin(c)) + 2*nix*(pix*cos(a)*cos(b) + piz*cos(a)*cos(c)*sin(b) + piy*cos(a)*sin(b)*sin(c)))*(nix*(x - qix - piy*(cos(c)*sin(a) - cos(a)*sin(b)*sin(c)) + piz*(sin(a)*sin(c) + cos(a)*cos(c)*sin(b)) + pix*cos(a)*cos(b)) + niy*(y - qiy + piy*(cos(a)*cos(c) + sin(a)*sin(b)*sin(c)) - piz*(cos(a)*sin(c) - cos(c)*sin(a)*sin(b)) + pix*cos(b)*sin(a)) + niz*(z - qiz - pix*sin(b) + piz*cos(b)*cos(c) + piy*cos(b)*sin(c)));
        d2J_dc2 = (nix*(piy*(sin(a)*sin(c) + cos(a)*cos(c)*sin(b)) + piz*(cos(c)*sin(a) - cos(a)*sin(b)*sin(c))) - niy*(piy*(cos(a)*sin(c) - cos(c)*sin(a)*sin(b)) + piz*(cos(a)*cos(c) + sin(a)*sin(b)*sin(c))) + niz*(piy*cos(b)*cos(c) - piz*cos(b)*sin(c)))*(2*nix*(piy*(sin(a)*sin(c) + cos(a)*cos(c)*sin(b)) + piz*(cos(c)*sin(a) - cos(a)*sin(b)*sin(c))) - 2*niy*(piy*(cos(a)*sin(c) - cos(c)*sin(a)*sin(b)) + piz*(cos(a)*cos(c) + sin(a)*sin(b)*sin(c))) + 2*niz*(piy*cos(b)*cos(c) - piz*cos(b)*sin(c))) - (2*niy*(piy*(cos(a)*cos(c) + sin(a)*sin(b)*sin(c)) - piz*(cos(a)*sin(c) - cos(c)*sin(a)*sin(b))) - 2*nix*(piy*(cos(c)*sin(a) - cos(a)*sin(b)*sin(c)) - piz*(sin(a)*sin(c) + cos(a)*cos(c)*sin(b))) + 2*niz*(piz*cos(b)*cos(c) + piy*cos(b)*sin(c)))*(nix*(x - qix - piy*(cos(c)*sin(a) - cos(a)*sin(b)*sin(c)) + piz*(sin(a)*sin(c) + cos(a)*cos(c)*sin(b)) + pix*cos(a)*cos(b)) + niy*(y - qiy + piy*(cos(a)*cos(c) + sin(a)*sin(b)*sin(c)) - piz*(cos(a)*sin(c) - cos(c)*sin(a)*sin(b)) + pix*cos(b)*sin(a)) + niz*(z - qiz - pix*sin(b) + piz*cos(b)*cos(c) + piy*cos(b)*sin(c)));
        d2J_dxda = nix*(2*niy*(piz*(sin(a)*sin(c) + cos(a)*cos(c)*sin(b)) - piy*(cos(c)*sin(a) - cos(a)*sin(b)*sin(c)) + pix*cos(a)*cos(b)) - 2*nix*(piy*(cos(a)*cos(c) + sin(a)*sin(b)*sin(c)) - piz*(cos(a)*sin(c) - cos(c)*sin(a)*sin(b)) + pix*cos(b)*sin(a)));
        d2J_dadx = 2*nix*(niy*(piz*(sin(a)*sin(c) + cos(a)*cos(c)*sin(b)) - piy*(cos(c)*sin(a) - cos(a)*sin(b)*sin(c)) + pix*cos(a)*cos(b)) - nix*(piy*(cos(a)*cos(c) + sin(a)*sin(b)*sin(c)) - piz*(cos(a)*sin(c) - cos(c)*sin(a)*sin(b)) + pix*cos(b)*sin(a)));
        d2J_dyda = niy*(2*niy*(piz*(sin(a)*sin(c) + cos(a)*cos(c)*sin(b)) - piy*(cos(c)*sin(a) - cos(a)*sin(b)*sin(c)) + pix*cos(a)*cos(b)) - 2*nix*(piy*(cos(a)*cos(c) + sin(a)*sin(b)*sin(c)) - piz*(cos(a)*sin(c) - cos(c)*sin(a)*sin(b)) + pix*cos(b)*sin(a)));
        d2J_dady = 2*niy*(niy*(piz*(sin(a)*sin(c) + cos(a)*cos(c)*sin(b)) - piy*(cos(c)*sin(a) - cos(a)*sin(b)*sin(c)) + pix*cos(a)*cos(b)) - nix*(piy*(cos(a)*cos(c) + sin(a)*sin(b)*sin(c)) - piz*(cos(a)*sin(c) - cos(c)*sin(a)*sin(b)) + pix*cos(b)*sin(a)));
        d2J_dzda = niz*(2*niy*(piz*(sin(a)*sin(c) + cos(a)*cos(c)*sin(b)) - piy*(cos(c)*sin(a) - cos(a)*sin(b)*sin(c)) + pix*cos(a)*cos(b)) - 2*nix*(piy*(cos(a)*cos(c) + sin(a)*sin(b)*sin(c)) - piz*(cos(a)*sin(c) - cos(c)*sin(a)*sin(b)) + pix*cos(b)*sin(a)));
        d2J_dadz = 2*niz*(niy*(piz*(sin(a)*sin(c) + cos(a)*cos(c)*sin(b)) - piy*(cos(c)*sin(a) - cos(a)*sin(b)*sin(c)) + pix*cos(a)*cos(b)) - nix*(piy*(cos(a)*cos(c) + sin(a)*sin(b)*sin(c)) - piz*(cos(a)*sin(c) - cos(c)*sin(a)*sin(b)) + pix*cos(b)*sin(a)));
        d2J_dxdb = nix*(2*niy*(piz*cos(b)*cos(c)*sin(a) - pix*sin(a)*sin(b) + piy*cos(b)*sin(a)*sin(c)) - 2*niz*(pix*cos(b) + piz*cos(c)*sin(b) + piy*sin(b)*sin(c)) + 2*nix*(piz*cos(a)*cos(b)*cos(c) - pix*cos(a)*sin(b) + piy*cos(a)*cos(b)*sin(c)));
        d2J_dbdx = 2*nix*(niy*(piz*cos(b)*cos(c)*sin(a) - pix*sin(a)*sin(b) + piy*cos(b)*sin(a)*sin(c)) - niz*(pix*cos(b) + piz*cos(c)*sin(b) + piy*sin(b)*sin(c)) + nix*(piz*cos(a)*cos(b)*cos(c) - pix*cos(a)*sin(b) + piy*cos(a)*cos(b)*sin(c)));
        d2J_dydb = niy*(2*niy*(piz*cos(b)*cos(c)*sin(a) - pix*sin(a)*sin(b) + piy*cos(b)*sin(a)*sin(c)) - 2*niz*(pix*cos(b) + piz*cos(c)*sin(b) + piy*sin(b)*sin(c)) + 2*nix*(piz*cos(a)*cos(b)*cos(c) - pix*cos(a)*sin(b) + piy*cos(a)*cos(b)*sin(c)));
        d2J_dbdy = 2*niy*(niy*(piz*cos(b)*cos(c)*sin(a) - pix*sin(a)*sin(b) + piy*cos(b)*sin(a)*sin(c)) - niz*(pix*cos(b) + piz*cos(c)*sin(b) + piy*sin(b)*sin(c)) + nix*(piz*cos(a)*cos(b)*cos(c) - pix*cos(a)*sin(b) + piy*cos(a)*cos(b)*sin(c)));
        d2J_dzdb = niz*(2*niy*(piz*cos(b)*cos(c)*sin(a) - pix*sin(a)*sin(b) + piy*cos(b)*sin(a)*sin(c)) - 2*niz*(pix*cos(b) + piz*cos(c)*sin(b) + piy*sin(b)*sin(c)) + 2*nix*(piz*cos(a)*cos(b)*cos(c) - pix*cos(a)*sin(b) + piy*cos(a)*cos(b)*sin(c)));
        d2J_dbdz = 2*niz*(niy*(piz*cos(b)*cos(c)*sin(a) - pix*sin(a)*sin(b) + piy*cos(b)*sin(a)*sin(c)) - niz*(pix*cos(b) + piz*cos(c)*sin(b) + piy*sin(b)*sin(c)) + nix*(piz*cos(a)*cos(b)*cos(c) - pix*cos(a)*sin(b) + piy*cos(a)*cos(b)*sin(c)));
        d2J_dxdc = nix*(2*nix*(piy*(sin(a)*sin(c) + cos(a)*cos(c)*sin(b)) + piz*(cos(c)*sin(a) - cos(a)*sin(b)*sin(c))) - 2*niy*(piy*(cos(a)*sin(c) - cos(c)*sin(a)*sin(b)) + piz*(cos(a)*cos(c) + sin(a)*sin(b)*sin(c))) + 2*niz*(piy*cos(b)*cos(c) - piz*cos(b)*sin(c)));
        d2J_dcdx = 2*nix*(nix*(piy*(sin(a)*sin(c) + cos(a)*cos(c)*sin(b)) + piz*(cos(c)*sin(a) - cos(a)*sin(b)*sin(c))) - niy*(piy*(cos(a)*sin(c) - cos(c)*sin(a)*sin(b)) + piz*(cos(a)*cos(c) + sin(a)*sin(b)*sin(c))) + niz*(piy*cos(b)*cos(c) - piz*cos(b)*sin(c)));
        d2J_dydc = niy*(2*nix*(piy*(sin(a)*sin(c) + cos(a)*cos(c)*sin(b)) + piz*(cos(c)*sin(a) - cos(a)*sin(b)*sin(c))) - 2*niy*(piy*(cos(a)*sin(c) - cos(c)*sin(a)*sin(b)) + piz*(cos(a)*cos(c) + sin(a)*sin(b)*sin(c))) + 2*niz*(piy*cos(b)*cos(c) - piz*cos(b)*sin(c)));
        d2J_dcdy = 2*niy*(nix*(piy*(sin(a)*sin(c) + cos(a)*cos(c)*sin(b)) + piz*(cos(c)*sin(a) - cos(a)*sin(b)*sin(c))) - niy*(piy*(cos(a)*sin(c) - cos(c)*sin(a)*sin(b)) + piz*(cos(a)*cos(c) + sin(a)*sin(b)*sin(c))) + niz*(piy*cos(b)*cos(c) - piz*cos(b)*sin(c)));
        d2J_dzdc = niz*(2*nix*(piy*(sin(a)*sin(c) + cos(a)*cos(c)*sin(b)) + piz*(cos(c)*sin(a) - cos(a)*sin(b)*sin(c))) - 2*niy*(piy*(cos(a)*sin(c) - cos(c)*sin(a)*sin(b)) + piz*(cos(a)*cos(c) + sin(a)*sin(b)*sin(c))) + 2*niz*(piy*cos(b)*cos(c) - piz*cos(b)*sin(c)));
        d2J_dcdz = 2*niz*(nix*(piy*(sin(a)*sin(c) + cos(a)*cos(c)*sin(b)) + piz*(cos(c)*sin(a) - cos(a)*sin(b)*sin(c))) - niy*(piy*(cos(a)*sin(c) - cos(c)*sin(a)*sin(b)) + piz*(cos(a)*cos(c) + sin(a)*sin(b)*sin(c))) + niz*(piy*cos(b)*cos(c) - piz*cos(b)*sin(c)));
        d2J_dadb = (niy*(piz*(sin(a)*sin(c) + cos(a)*cos(c)*sin(b)) - piy*(cos(c)*sin(a) - cos(a)*sin(b)*sin(c)) + pix*cos(a)*cos(b)) - nix*(piy*(cos(a)*cos(c) + sin(a)*sin(b)*sin(c)) - piz*(cos(a)*sin(c) - cos(c)*sin(a)*sin(b)) + pix*cos(b)*sin(a)))*(2*niy*(piz*cos(b)*cos(c)*sin(a) - pix*sin(a)*sin(b) + piy*cos(b)*sin(a)*sin(c)) - 2*niz*(pix*cos(b) + piz*cos(c)*sin(b) + piy*sin(b)*sin(c)) + 2*nix*(piz*cos(a)*cos(b)*cos(c) - pix*cos(a)*sin(b) + piy*cos(a)*cos(b)*sin(c))) - (2*nix*(piz*cos(b)*cos(c)*sin(a) - pix*sin(a)*sin(b) + piy*cos(b)*sin(a)*sin(c)) - 2*niy*(piz*cos(a)*cos(b)*cos(c) - pix*cos(a)*sin(b) + piy*cos(a)*cos(b)*sin(c)))*(nix*(x - qix - piy*(cos(c)*sin(a) - cos(a)*sin(b)*sin(c)) + piz*(sin(a)*sin(c) + cos(a)*cos(c)*sin(b)) + pix*cos(a)*cos(b)) + niy*(y - qiy + piy*(cos(a)*cos(c) + sin(a)*sin(b)*sin(c)) - piz*(cos(a)*sin(c) - cos(c)*sin(a)*sin(b)) + pix*cos(b)*sin(a)) + niz*(z - qiz - pix*sin(b) + piz*cos(b)*cos(c) + piy*cos(b)*sin(c)));
        d2J_dbda = (2*niy*(piz*(sin(a)*sin(c) + cos(a)*cos(c)*sin(b)) - piy*(cos(c)*sin(a) - cos(a)*sin(b)*sin(c)) + pix*cos(a)*cos(b)) - 2*nix*(piy*(cos(a)*cos(c) + sin(a)*sin(b)*sin(c)) - piz*(cos(a)*sin(c) - cos(c)*sin(a)*sin(b)) + pix*cos(b)*sin(a)))*(niy*(piz*cos(b)*cos(c)*sin(a) - pix*sin(a)*sin(b) + piy*cos(b)*sin(a)*sin(c)) - niz*(pix*cos(b) + piz*cos(c)*sin(b) + piy*sin(b)*sin(c)) + nix*(piz*cos(a)*cos(b)*cos(c) - pix*cos(a)*sin(b) + piy*cos(a)*cos(b)*sin(c))) - (2*nix*(piz*cos(b)*cos(c)*sin(a) - pix*sin(a)*sin(b) + piy*cos(b)*sin(a)*sin(c)) - 2*niy*(piz*cos(a)*cos(b)*cos(c) - pix*cos(a)*sin(b) + piy*cos(a)*cos(b)*sin(c)))*(nix*(x - qix - piy*(cos(c)*sin(a) - cos(a)*sin(b)*sin(c)) + piz*(sin(a)*sin(c) + cos(a)*cos(c)*sin(b)) + pix*cos(a)*cos(b)) + niy*(y - qiy + piy*(cos(a)*cos(c) + sin(a)*sin(b)*sin(c)) - piz*(cos(a)*sin(c) - cos(c)*sin(a)*sin(b)) + pix*cos(b)*sin(a)) + niz*(z - qiz - pix*sin(b) + piz*cos(b)*cos(c) + piy*cos(b)*sin(c)));
        d2J_dbdc = (2*nix*(piy*cos(a)*cos(b)*cos(c) - piz*cos(a)*cos(b)*sin(c)) - 2*niz*(piy*cos(c)*sin(b) - piz*sin(b)*sin(c)) + 2*niy*(piy*cos(b)*cos(c)*sin(a) - piz*cos(b)*sin(a)*sin(c)))*(nix*(x - qix - piy*(cos(c)*sin(a) - cos(a)*sin(b)*sin(c)) + piz*(sin(a)*sin(c) + cos(a)*cos(c)*sin(b)) + pix*cos(a)*cos(b)) + niy*(y - qiy + piy*(cos(a)*cos(c) + sin(a)*sin(b)*sin(c)) - piz*(cos(a)*sin(c) - cos(c)*sin(a)*sin(b)) + pix*cos(b)*sin(a)) + niz*(z - qiz - pix*sin(b) + piz*cos(b)*cos(c) + piy*cos(b)*sin(c))) + (niy*(piz*cos(b)*cos(c)*sin(a) - pix*sin(a)*sin(b) + piy*cos(b)*sin(a)*sin(c)) - niz*(pix*cos(b) + piz*cos(c)*sin(b) + piy*sin(b)*sin(c)) + nix*(piz*cos(a)*cos(b)*cos(c) - pix*cos(a)*sin(b) + piy*cos(a)*cos(b)*sin(c)))*(2*nix*(piy*(sin(a)*sin(c) + cos(a)*cos(c)*sin(b)) + piz*(cos(c)*sin(a) - cos(a)*sin(b)*sin(c))) - 2*niy*(piy*(cos(a)*sin(c) - cos(c)*sin(a)*sin(b)) + piz*(cos(a)*cos(c) + sin(a)*sin(b)*sin(c))) + 2*niz*(piy*cos(b)*cos(c) - piz*cos(b)*sin(c)));
        d2J_dcdb = (2*nix*(piy*cos(a)*cos(b)*cos(c) - piz*cos(a)*cos(b)*sin(c)) - 2*niz*(piy*cos(c)*sin(b) - piz*sin(b)*sin(c)) + 2*niy*(piy*cos(b)*cos(c)*sin(a) - piz*cos(b)*sin(a)*sin(c)))*(nix*(x - qix - piy*(cos(c)*sin(a) - cos(a)*sin(b)*sin(c)) + piz*(sin(a)*sin(c) + cos(a)*cos(c)*sin(b)) + pix*cos(a)*cos(b)) + niy*(y - qiy + piy*(cos(a)*cos(c) + sin(a)*sin(b)*sin(c)) - piz*(cos(a)*sin(c) - cos(c)*sin(a)*sin(b)) + pix*cos(b)*sin(a)) + niz*(z - qiz - pix*sin(b) + piz*cos(b)*cos(c) + piy*cos(b)*sin(c))) + (2*niy*(piz*cos(b)*cos(c)*sin(a) - pix*sin(a)*sin(b) + piy*cos(b)*sin(a)*sin(c)) - 2*niz*(pix*cos(b) + piz*cos(c)*sin(b) + piy*sin(b)*sin(c)) + 2*nix*(piz*cos(a)*cos(b)*cos(c) - pix*cos(a)*sin(b) + piy*cos(a)*cos(b)*sin(c)))*(nix*(piy*(sin(a)*sin(c) + cos(a)*cos(c)*sin(b)) + piz*(cos(c)*sin(a) - cos(a)*sin(b)*sin(c))) - niy*(piy*(cos(a)*sin(c) - cos(c)*sin(a)*sin(b)) + piz*(cos(a)*cos(c) + sin(a)*sin(b)*sin(c))) + niz*(piy*cos(b)*cos(c) - piz*cos(b)*sin(c)));
        d2J_dcda = (2*niy*(piz*(sin(a)*sin(c) + cos(a)*cos(c)*sin(b)) - piy*(cos(c)*sin(a) - cos(a)*sin(b)*sin(c)) + pix*cos(a)*cos(b)) - 2*nix*(piy*(cos(a)*cos(c) + sin(a)*sin(b)*sin(c)) - piz*(cos(a)*sin(c) - cos(c)*sin(a)*sin(b)) + pix*cos(b)*sin(a)))*(nix*(piy*(sin(a)*sin(c) + cos(a)*cos(c)*sin(b)) + piz*(cos(c)*sin(a) - cos(a)*sin(b)*sin(c))) - niy*(piy*(cos(a)*sin(c) - cos(c)*sin(a)*sin(b)) + piz*(cos(a)*cos(c) + sin(a)*sin(b)*sin(c))) + niz*(piy*cos(b)*cos(c) - piz*cos(b)*sin(c))) + (2*nix*(piy*(cos(a)*sin(c) - cos(c)*sin(a)*sin(b)) + piz*(cos(a)*cos(c) + sin(a)*sin(b)*sin(c))) + 2*niy*(piy*(sin(a)*sin(c) + cos(a)*cos(c)*sin(b)) + piz*(cos(c)*sin(a) - cos(a)*sin(b)*sin(c))))*(nix*(x - qix - piy*(cos(c)*sin(a) - cos(a)*sin(b)*sin(c)) + piz*(sin(a)*sin(c) + cos(a)*cos(c)*sin(b)) + pix*cos(a)*cos(b)) + niy*(y - qiy + piy*(cos(a)*cos(c) + sin(a)*sin(b)*sin(c)) - piz*(cos(a)*sin(c) - cos(c)*sin(a)*sin(b)) + pix*cos(b)*sin(a)) + niz*(z - qiz - pix*sin(b) + piz*cos(b)*cos(c) + piy*cos(b)*sin(c)));
        d2J_dadc = (niy*(piz*(sin(a)*sin(c) + cos(a)*cos(c)*sin(b)) - piy*(cos(c)*sin(a) - cos(a)*sin(b)*sin(c)) + pix*cos(a)*cos(b)) - nix*(piy*(cos(a)*cos(c) + sin(a)*sin(b)*sin(c)) - piz*(cos(a)*sin(c) - cos(c)*sin(a)*sin(b)) + pix*cos(b)*sin(a)))*(2*nix*(piy*(sin(a)*sin(c) + cos(a)*cos(c)*sin(b)) + piz*(cos(c)*sin(a) - cos(a)*sin(b)*sin(c))) - 2*niy*(piy*(cos(a)*sin(c) - cos(c)*sin(a)*sin(b)) + piz*(cos(a)*cos(c) + sin(a)*sin(b)*sin(c))) + 2*niz*(piy*cos(b)*cos(c) - piz*cos(b)*sin(c))) + (2*nix*(piy*(cos(a)*sin(c) - cos(c)*sin(a)*sin(b)) + piz*(cos(a)*cos(c) + sin(a)*sin(b)*sin(c))) + 2*niy*(piy*(sin(a)*sin(c) + cos(a)*cos(c)*sin(b)) + piz*(cos(c)*sin(a) - cos(a)*sin(b)*sin(c))))*(nix*(x - qix - piy*(cos(c)*sin(a) - cos(a)*sin(b)*sin(c)) + piz*(sin(a)*sin(c) + cos(a)*cos(c)*sin(b)) + pix*cos(a)*cos(b)) + niy*(y - qiy + piy*(cos(a)*cos(c) + sin(a)*sin(b)*sin(c)) - piz*(cos(a)*sin(c) - cos(c)*sin(a)*sin(b)) + pix*cos(b)*sin(a)) + niz*(z - qiz - pix*sin(b) + piz*cos(b)*cos(c) + piy*cos(b)*sin(c)));

        // Finally reconstruct the matrix
        Matrix d2J_dX2_temp(6,6);
        d2J_dX2_temp << d2J_dx2,     d2J_dydx,	  d2J_dzdx,   d2J_dadx,   d2J_dbdx,     d2J_dcdx,
                d2J_dxdy,    d2J_dy2,	  d2J_dzdy,   d2J_dady,   d2J_dbdy,     d2J_dcdy,
                d2J_dxdz,    d2J_dydz,    d2J_dz2,    d2J_dadz,   d2J_dbdz,     d2J_dcdz,
                d2J_dxda,    d2J_dyda,    d2J_dzda,   d2J_da2,	  d2J_dbda,     d2J_dcda,
                d2J_dxdb,    d2J_dydb,    d2J_dzdb,   d2J_dadb,   d2J_db2,      d2J_dcdb,
                d2J_dxdc,    d2J_dydc,    d2J_dzdc,   d2J_dadc,   d2J_dbdc,     d2J_dc2;
        d2J_dX2 += d2J_dX2_temp;

        /***********************************************************************
        d2J_dZdX -- X is the [R|T] in the form of [x, y, z, a, b, c]
        x, y, z is the translation part
        a, b, c is the rotation part in Euler format
        [x, y, z, a, b, c] is acquired from the Transformation Matrix returned by ICP.
        *************************************************************************/

        double 	d2J_dpix_dx,    d2J_dpiy_dx,	d2J_dpiz_dx,  	   d2J_dqix_dx,    d2J_dqiy_dx,	   d2J_dqiz_dx,
                d2J_dpix_dy,    d2J_dpiy_dy,	d2J_dpiz_dy,   	   d2J_dqix_dy,    d2J_dqiy_dy,	   d2J_dqiz_dy,
                d2J_dpix_dz,    d2J_dpiy_dz,    d2J_dpiz_dz,       d2J_dqix_dz,    d2J_dqiy_dz,    d2J_dqiz_dz,
                d2J_dpix_da,    d2J_dpiy_da,    d2J_dpiz_da,       d2J_dqix_da,    d2J_dqiy_da,    d2J_dqiz_da,
                d2J_dpix_db,    d2J_dpiy_db,    d2J_dpiz_db,       d2J_dqix_db,    d2J_dqiy_db,    d2J_dqiz_db,
                d2J_dpix_dc,    d2J_dpiy_dc,    d2J_dpiz_dc,       d2J_dqix_dc,    d2J_dqiy_dc,    d2J_dqiz_dc;

        d2J_dpix_dx = 2*nix*(nix*cos(a)*cos(b) - niz*sin(b) + niy*cos(b)*sin(a));
        d2J_dpix_dy = 2*niy*(nix*cos(a)*cos(b) - niz*sin(b) + niy*cos(b)*sin(a));
        d2J_dpix_dz = 2*niz*(nix*cos(a)*cos(b) - niz*sin(b) + niy*cos(b)*sin(a));
        d2J_dpix_da = (2*niy*cos(a)*cos(b) - 2*nix*cos(b)*sin(a))*(nix*(x - qix - piy*(cos(c)*sin(a) - cos(a)*sin(b)*sin(c)) + piz*(sin(a)*sin(c) + cos(a)*cos(c)*sin(b)) + pix*cos(a)*cos(b)) + niy*(y - qiy + piy*(cos(a)*cos(c) + sin(a)*sin(b)*sin(c)) - piz*(cos(a)*sin(c) - cos(c)*sin(a)*sin(b)) + pix*cos(b)*sin(a)) + niz*(z - qiz - pix*sin(b) + piz*cos(b)*cos(c) + piy*cos(b)*sin(c))) + (2*niy*(piz*(sin(a)*sin(c) + cos(a)*cos(c)*sin(b)) - piy*(cos(c)*sin(a) - cos(a)*sin(b)*sin(c)) + pix*cos(a)*cos(b)) - 2*nix*(piy*(cos(a)*cos(c) + sin(a)*sin(b)*sin(c)) - piz*(cos(a)*sin(c) - cos(c)*sin(a)*sin(b)) + pix*cos(b)*sin(a)))*(nix*cos(a)*cos(b) - niz*sin(b) + niy*cos(b)*sin(a));
        d2J_dpix_db = (2*niy*(piz*cos(b)*cos(c)*sin(a) - pix*sin(a)*sin(b) + piy*cos(b)*sin(a)*sin(c)) - 2*niz*(pix*cos(b) + piz*cos(c)*sin(b) + piy*sin(b)*sin(c)) + 2*nix*(piz*cos(a)*cos(b)*cos(c) - pix*cos(a)*sin(b) + piy*cos(a)*cos(b)*sin(c)))*(nix*cos(a)*cos(b) - niz*sin(b) + niy*cos(b)*sin(a)) - (2*niz*cos(b) + 2*nix*cos(a)*sin(b) + 2*niy*sin(a)*sin(b))*(nix*(x - qix - piy*(cos(c)*sin(a) - cos(a)*sin(b)*sin(c)) + piz*(sin(a)*sin(c) + cos(a)*cos(c)*sin(b)) + pix*cos(a)*cos(b)) + niy*(y - qiy + piy*(cos(a)*cos(c) + sin(a)*sin(b)*sin(c)) - piz*(cos(a)*sin(c) - cos(c)*sin(a)*sin(b)) + pix*cos(b)*sin(a)) + niz*(z - qiz - pix*sin(b) + piz*cos(b)*cos(c) + piy*cos(b)*sin(c)));
        d2J_dpix_dc = (2*nix*(piy*(sin(a)*sin(c) + cos(a)*cos(c)*sin(b)) + piz*(cos(c)*sin(a) - cos(a)*sin(b)*sin(c))) - 2*niy*(piy*(cos(a)*sin(c) - cos(c)*sin(a)*sin(b)) + piz*(cos(a)*cos(c) + sin(a)*sin(b)*sin(c))) + 2*niz*(piy*cos(b)*cos(c) - piz*cos(b)*sin(c)))*(nix*cos(a)*cos(b) - niz*sin(b) + niy*cos(b)*sin(a));
        d2J_dpiy_dx = 2*nix*(niy*(cos(a)*cos(c) + sin(a)*sin(b)*sin(c)) - nix*(cos(c)*sin(a) - cos(a)*sin(b)*sin(c)) + niz*cos(b)*sin(c));
        d2J_dpiy_dy = 2*niy*(niy*(cos(a)*cos(c) + sin(a)*sin(b)*sin(c)) - nix*(cos(c)*sin(a) - cos(a)*sin(b)*sin(c)) + niz*cos(b)*sin(c));
        d2J_dpiy_dz = 2*niz*(niy*(cos(a)*cos(c) + sin(a)*sin(b)*sin(c)) - nix*(cos(c)*sin(a) - cos(a)*sin(b)*sin(c)) + niz*cos(b)*sin(c));
        d2J_dpiy_da = (2*niy*(piz*(sin(a)*sin(c) + cos(a)*cos(c)*sin(b)) - piy*(cos(c)*sin(a) - cos(a)*sin(b)*sin(c)) + pix*cos(a)*cos(b)) - 2*nix*(piy*(cos(a)*cos(c) + sin(a)*sin(b)*sin(c)) - piz*(cos(a)*sin(c) - cos(c)*sin(a)*sin(b)) + pix*cos(b)*sin(a)))*(niy*(cos(a)*cos(c) + sin(a)*sin(b)*sin(c)) - nix*(cos(c)*sin(a) - cos(a)*sin(b)*sin(c)) + niz*cos(b)*sin(c)) - (2*nix*(cos(a)*cos(c) + sin(a)*sin(b)*sin(c)) + 2*niy*(cos(c)*sin(a) - cos(a)*sin(b)*sin(c)))*(nix*(x - qix - piy*(cos(c)*sin(a) - cos(a)*sin(b)*sin(c)) + piz*(sin(a)*sin(c) + cos(a)*cos(c)*sin(b)) + pix*cos(a)*cos(b)) + niy*(y - qiy + piy*(cos(a)*cos(c) + sin(a)*sin(b)*sin(c)) - piz*(cos(a)*sin(c) - cos(c)*sin(a)*sin(b)) + pix*cos(b)*sin(a)) + niz*(z - qiz - pix*sin(b) + piz*cos(b)*cos(c) + piy*cos(b)*sin(c)));
        d2J_dpiy_db = (2*niy*(piz*cos(b)*cos(c)*sin(a) - pix*sin(a)*sin(b) + piy*cos(b)*sin(a)*sin(c)) - 2*niz*(pix*cos(b) + piz*cos(c)*sin(b) + piy*sin(b)*sin(c)) + 2*nix*(piz*cos(a)*cos(b)*cos(c) - pix*cos(a)*sin(b) + piy*cos(a)*cos(b)*sin(c)))*(niy*(cos(a)*cos(c) + sin(a)*sin(b)*sin(c)) - nix*(cos(c)*sin(a) - cos(a)*sin(b)*sin(c)) + niz*cos(b)*sin(c)) + (2*nix*cos(a)*cos(b)*sin(c) - 2*niz*sin(b)*sin(c) + 2*niy*cos(b)*sin(a)*sin(c))*(nix*(x - qix - piy*(cos(c)*sin(a) - cos(a)*sin(b)*sin(c)) + piz*(sin(a)*sin(c) + cos(a)*cos(c)*sin(b)) + pix*cos(a)*cos(b)) + niy*(y - qiy + piy*(cos(a)*cos(c) + sin(a)*sin(b)*sin(c)) - piz*(cos(a)*sin(c) - cos(c)*sin(a)*sin(b)) + pix*cos(b)*sin(a)) + niz*(z - qiz - pix*sin(b) + piz*cos(b)*cos(c) + piy*cos(b)*sin(c)));
        d2J_dpiy_dc = (2*nix*(piy*(sin(a)*sin(c) + cos(a)*cos(c)*sin(b)) + piz*(cos(c)*sin(a) - cos(a)*sin(b)*sin(c))) - 2*niy*(piy*(cos(a)*sin(c) - cos(c)*sin(a)*sin(b)) + piz*(cos(a)*cos(c) + sin(a)*sin(b)*sin(c))) + 2*niz*(piy*cos(b)*cos(c) - piz*cos(b)*sin(c)))*(niy*(cos(a)*cos(c) + sin(a)*sin(b)*sin(c)) - nix*(cos(c)*sin(a) - cos(a)*sin(b)*sin(c)) + niz*cos(b)*sin(c)) + (2*nix*(sin(a)*sin(c) + cos(a)*cos(c)*sin(b)) - 2*niy*(cos(a)*sin(c) - cos(c)*sin(a)*sin(b)) + 2*niz*cos(b)*cos(c))*(nix*(x - qix - piy*(cos(c)*sin(a) - cos(a)*sin(b)*sin(c)) + piz*(sin(a)*sin(c) + cos(a)*cos(c)*sin(b)) + pix*cos(a)*cos(b)) + niy*(y - qiy + piy*(cos(a)*cos(c) + sin(a)*sin(b)*sin(c)) - piz*(cos(a)*sin(c) - cos(c)*sin(a)*sin(b)) + pix*cos(b)*sin(a)) + niz*(z - qiz - pix*sin(b) + piz*cos(b)*cos(c) + piy*cos(b)*sin(c)));
        d2J_dpiz_dx = 2*nix*(nix*(sin(a)*sin(c) + cos(a)*cos(c)*sin(b)) - niy*(cos(a)*sin(c) - cos(c)*sin(a)*sin(b)) + niz*cos(b)*cos(c));
        d2J_dpiz_dy = 2*niy*(nix*(sin(a)*sin(c) + cos(a)*cos(c)*sin(b)) - niy*(cos(a)*sin(c) - cos(c)*sin(a)*sin(b)) + niz*cos(b)*cos(c));
        d2J_dpiz_dz = 2*niz*(nix*(sin(a)*sin(c) + cos(a)*cos(c)*sin(b)) - niy*(cos(a)*sin(c) - cos(c)*sin(a)*sin(b)) + niz*cos(b)*cos(c));
        d2J_dpiz_da = (2*niy*(piz*(sin(a)*sin(c) + cos(a)*cos(c)*sin(b)) - piy*(cos(c)*sin(a) - cos(a)*sin(b)*sin(c)) + pix*cos(a)*cos(b)) - 2*nix*(piy*(cos(a)*cos(c) + sin(a)*sin(b)*sin(c)) - piz*(cos(a)*sin(c) - cos(c)*sin(a)*sin(b)) + pix*cos(b)*sin(a)))*(nix*(sin(a)*sin(c) + cos(a)*cos(c)*sin(b)) - niy*(cos(a)*sin(c) - cos(c)*sin(a)*sin(b)) + niz*cos(b)*cos(c)) + (2*nix*(cos(a)*sin(c) - cos(c)*sin(a)*sin(b)) + 2*niy*(sin(a)*sin(c) + cos(a)*cos(c)*sin(b)))*(nix*(x - qix - piy*(cos(c)*sin(a) - cos(a)*sin(b)*sin(c)) + piz*(sin(a)*sin(c) + cos(a)*cos(c)*sin(b)) + pix*cos(a)*cos(b)) + niy*(y - qiy + piy*(cos(a)*cos(c) + sin(a)*sin(b)*sin(c)) - piz*(cos(a)*sin(c) - cos(c)*sin(a)*sin(b)) + pix*cos(b)*sin(a)) + niz*(z - qiz - pix*sin(b) + piz*cos(b)*cos(c) + piy*cos(b)*sin(c)));
        d2J_dpiz_db = (2*niy*(piz*cos(b)*cos(c)*sin(a) - pix*sin(a)*sin(b) + piy*cos(b)*sin(a)*sin(c)) - 2*niz*(pix*cos(b) + piz*cos(c)*sin(b) + piy*sin(b)*sin(c)) + 2*nix*(piz*cos(a)*cos(b)*cos(c) - pix*cos(a)*sin(b) + piy*cos(a)*cos(b)*sin(c)))*(nix*(sin(a)*sin(c) + cos(a)*cos(c)*sin(b)) - niy*(cos(a)*sin(c) - cos(c)*sin(a)*sin(b)) + niz*cos(b)*cos(c)) + (2*nix*cos(a)*cos(b)*cos(c) - 2*niz*cos(c)*sin(b) + 2*niy*cos(b)*cos(c)*sin(a))*(nix*(x - qix - piy*(cos(c)*sin(a) - cos(a)*sin(b)*sin(c)) + piz*(sin(a)*sin(c) + cos(a)*cos(c)*sin(b)) + pix*cos(a)*cos(b)) + niy*(y - qiy + piy*(cos(a)*cos(c) + sin(a)*sin(b)*sin(c)) - piz*(cos(a)*sin(c) - cos(c)*sin(a)*sin(b)) + pix*cos(b)*sin(a)) + niz*(z - qiz - pix*sin(b) + piz*cos(b)*cos(c) + piy*cos(b)*sin(c)));
        d2J_dpiz_dc = (2*nix*(piy*(sin(a)*sin(c) + cos(a)*cos(c)*sin(b)) + piz*(cos(c)*sin(a) - cos(a)*sin(b)*sin(c))) - 2*niy*(piy*(cos(a)*sin(c) - cos(c)*sin(a)*sin(b)) + piz*(cos(a)*cos(c) + sin(a)*sin(b)*sin(c))) + 2*niz*(piy*cos(b)*cos(c) - piz*cos(b)*sin(c)))*(nix*(sin(a)*sin(c) + cos(a)*cos(c)*sin(b)) - niy*(cos(a)*sin(c) - cos(c)*sin(a)*sin(b)) + niz*cos(b)*cos(c)) - (2*niy*(cos(a)*cos(c) + sin(a)*sin(b)*sin(c)) - 2*nix*(cos(c)*sin(a) - cos(a)*sin(b)*sin(c)) + 2*niz*cos(b)*sin(c))*(nix*(x - qix - piy*(cos(c)*sin(a) - cos(a)*sin(b)*sin(c)) + piz*(sin(a)*sin(c) + cos(a)*cos(c)*sin(b)) + pix*cos(a)*cos(b)) + niy*(y - qiy + piy*(cos(a)*cos(c) + sin(a)*sin(b)*sin(c)) - piz*(cos(a)*sin(c) - cos(c)*sin(a)*sin(b)) + pix*cos(b)*sin(a)) + niz*(z - qiz - pix*sin(b) + piz*cos(b)*cos(c) + piy*cos(b)*sin(c)));
        d2J_dqix_dx = -2*pow(nix,2);
        d2J_dqix_dy = -2*nix*niy;
        d2J_dqix_dz = -2*nix*niz;
        d2J_dqix_da = -nix*(2*niy*(piz*(sin(a)*sin(c) + cos(a)*cos(c)*sin(b)) - piy*(cos(c)*sin(a) - cos(a)*sin(b)*sin(c)) + pix*cos(a)*cos(b)) - 2*nix*(piy*(cos(a)*cos(c) + sin(a)*sin(b)*sin(c)) - piz*(cos(a)*sin(c) - cos(c)*sin(a)*sin(b)) + pix*cos(b)*sin(a)));
        d2J_dqix_db = -nix*(2*niy*(piz*cos(b)*cos(c)*sin(a) - pix*sin(a)*sin(b) + piy*cos(b)*sin(a)*sin(c)) - 2*niz*(pix*cos(b) + piz*cos(c)*sin(b) + piy*sin(b)*sin(c)) + 2*nix*(piz*cos(a)*cos(b)*cos(c) - pix*cos(a)*sin(b) + piy*cos(a)*cos(b)*sin(c)));
        d2J_dqix_dc = -nix*(2*nix*(piy*(sin(a)*sin(c) + cos(a)*cos(c)*sin(b)) + piz*(cos(c)*sin(a) - cos(a)*sin(b)*sin(c))) - 2*niy*(piy*(cos(a)*sin(c) - cos(c)*sin(a)*sin(b)) + piz*(cos(a)*cos(c) + sin(a)*sin(b)*sin(c))) + 2*niz*(piy*cos(b)*cos(c) - piz*cos(b)*sin(c)));
        d2J_dqiy_dx = -2*nix*niy;
        d2J_dqiy_dy = -2*pow(niy,2);
        d2J_dqiy_dz = -2*niy*niz;
        d2J_dqiy_da = -niy*(2*niy*(piz*(sin(a)*sin(c) + cos(a)*cos(c)*sin(b)) - piy*(cos(c)*sin(a) - cos(a)*sin(b)*sin(c)) + pix*cos(a)*cos(b)) - 2*nix*(piy*(cos(a)*cos(c) + sin(a)*sin(b)*sin(c)) - piz*(cos(a)*sin(c) - cos(c)*sin(a)*sin(b)) + pix*cos(b)*sin(a)));
        d2J_dqiy_db = -niy*(2*niy*(piz*cos(b)*cos(c)*sin(a) - pix*sin(a)*sin(b) + piy*cos(b)*sin(a)*sin(c)) - 2*niz*(pix*cos(b) + piz*cos(c)*sin(b) + piy*sin(b)*sin(c)) + 2*nix*(piz*cos(a)*cos(b)*cos(c) - pix*cos(a)*sin(b) + piy*cos(a)*cos(b)*sin(c)));
        d2J_dqiy_dc = -niy*(2*nix*(piy*(sin(a)*sin(c) + cos(a)*cos(c)*sin(b)) + piz*(cos(c)*sin(a) - cos(a)*sin(b)*sin(c))) - 2*niy*(piy*(cos(a)*sin(c) - cos(c)*sin(a)*sin(b)) + piz*(cos(a)*cos(c) + sin(a)*sin(b)*sin(c))) + 2*niz*(piy*cos(b)*cos(c) - piz*cos(b)*sin(c)));
        d2J_dqiz_dx = -2*nix*niz;
        d2J_dqiz_dy = -2*niy*niz;
        d2J_dqiz_dz = -2*pow(niz,2);
        d2J_dqiz_da = -niz*(2*niy*(piz*(sin(a)*sin(c) + cos(a)*cos(c)*sin(b)) - piy*(cos(c)*sin(a) - cos(a)*sin(b)*sin(c)) + pix*cos(a)*cos(b)) - 2*nix*(piy*(cos(a)*cos(c) + sin(a)*sin(b)*sin(c)) - piz*(cos(a)*sin(c) - cos(c)*sin(a)*sin(b)) + pix*cos(b)*sin(a)));
        d2J_dqiz_db = -niz*(2*niy*(piz*cos(b)*cos(c)*sin(a) - pix*sin(a)*sin(b) + piy*cos(b)*sin(a)*sin(c)) - 2*niz*(pix*cos(b) + piz*cos(c)*sin(b) + piy*sin(b)*sin(c)) + 2*nix*(piz*cos(a)*cos(b)*cos(c) - pix*cos(a)*sin(b) + piy*cos(a)*cos(b)*sin(c)));
        d2J_dqiz_dc = -niz*(2*nix*(piy*(sin(a)*sin(c) + cos(a)*cos(c)*sin(b)) + piz*(cos(c)*sin(a) - cos(a)*sin(b)*sin(c))) - 2*niy*(piy*(cos(a)*sin(c) - cos(c)*sin(a)*sin(b)) + piz*(cos(a)*cos(c) + sin(a)*sin(b)*sin(c))) + 2*niz*(piy*cos(b)*cos(c) - piz*cos(b)*sin(c)));

        // Finally reconstruct the matrix
        Matrix d2J_dZdX_temp(6,6);
        d2J_dZdX_temp <<    d2J_dpix_dx,    d2J_dpiy_dx,	d2J_dpiz_dx,  	   d2J_dqix_dx,    d2J_dqiy_dx,	   d2J_dqiz_dx,
                d2J_dpix_dy,    d2J_dpiy_dy,	d2J_dpiz_dy,   	   d2J_dqix_dy,    d2J_dqiy_dy,	   d2J_dqiz_dy,
                d2J_dpix_dz,    d2J_dpiy_dz,        d2J_dpiz_dz,       d2J_dqix_dz,    d2J_dqiy_dz,    d2J_dqiz_dz,
                d2J_dpix_da,    d2J_dpiy_da,        d2J_dpiz_da,       d2J_dqix_da,    d2J_dqiy_da,    d2J_dqiz_da,
                d2J_dpix_db,    d2J_dpiy_db,        d2J_dpiz_db,       d2J_dqix_db,    d2J_dqiy_db,    d2J_dqiz_db,
                d2J_dpix_dc,    d2J_dpiy_dc,        d2J_dpiz_dc,       d2J_dqix_dc,    d2J_dqiy_dc,    d2J_dqiz_dc;
        d2J_dZdX.block(0,6*i,6,6) = d2J_dZdX_temp;

    }

    // Inverse of our state Hessian
    Matrix d2J_dX2_inv = d2J_dX2.inverse();

    // Calculate covariance
    Matrix covariance = d2J_dZdX * d2J_dZdX.transpose();
    covariance = d2J_dX2_inv * covariance * d2J_dX2_inv;

    // Multiply with our measurement noise
    return (sensorStdDev * sensorStdDev) * covariance;
}


template<typename T>
typename PointMatcher<T>::Matrix
PointToPlaneWithCovErrorMinimizer<T>::estimateCovariance2d(const ErrorElements& mPts, const TransformationParameters& transformation)
{

    // Create our empty Jacobian matices
    const int max_nbr_point = mPts.reading.getNbPoints();
    Matrix d2J_dX2(Matrix::Zero(3,3));
    Matrix d2J_dZdX(Matrix::Zero(3,4*max_nbr_point));

    // Normals of the reference cloud (taken to be true)
    Matrix normals = mPts.reference.getDescriptorViewByName("normals");

    // Return if we do not have any normals in DataPoints
    if (normals.rows() < 2)
        return std::numeric_limits<T>::max() * Matrix::Identity(3,3);

    T gamma = asin(transformation(1,0));
    T t_x = transformation(0,2);
    T t_y = transformation(1,2);

    // According to the rotation matrix I used and after verification
    // It is Yaw == [R] matrix used in the MatLab
    double x, y, a;
    x = t_x; y = t_y;
    a = gamma;

    // Loop through all points and calculate the state and partial Jacobians
    for(int i = 0; i < max_nbr_point; ++i)
    {

        double pix = mPts.reading.features(0,i);
        double piy = mPts.reading.features(1,i);
        double qix = mPts.reference.features(0,i);
        double qiy = mPts.reference.features(1,i);
        double nix = normals(0,i);
        double niy = normals(1,i);

        /***********************************************************************
        d2J_dX2 -- X is the [R|T] in the form of [x, y, a]
        x, y is the translation part
        a is the rotation part in Euler format (yaw)
        [x, y, a] is acquired from the Transformation Matrix returned by ICP.
        Now d2J_dX2 is a 3x3 matrix of the form
        d2J_dx2
        d2J_dxdy    d2J_dy2
        d2J_dxda    d2J_dyda    d2J_da2
        *************************************************************************/

        double 	d2J_dx2,     d2J_dydx,	  d2J_dadx,
                d2J_dxdy,    d2J_dy2,	  d2J_dady,
                d2J_dxda,    d2J_dyda,    d2J_da2;

        // These terms are generated from the provided Matlab scipts. We just have to copy
        // the expressions from the matlab output with two very simple changes.
        // The first one being the the sqaure of a number 'a' is shown as a^2 in matlab,
        // which is converted to pow(a,2) in the below expressions.
        // The second change is to add ';' at the end of each expression :)
        // In this way, matlab can be used to generate these terms for various objective functions of ICP
        // and they can simply be copied to the C++ files and with appropriate changes to ICP estimation,
        // its covariance can be easily estimated.

        d2J_dx2 = 2*pow(nix,2);
        d2J_dy2 = 2*pow(niy,2);
        d2J_dydx = 2*nix*niy;
        d2J_dxdy = 2*nix*niy;
        d2J_da2 = (2*nix*(pix*cos(a) - piy*sin(a)) + 2*niy*(piy*cos(a) + pix*sin(a)))*(nix*(qix - x - pix*cos(a) + piy*sin(a)) - niy*(y - qiy + piy*cos(a) + pix*sin(a))) + (nix*(piy*cos(a) + pix*sin(a)) - niy*(pix*cos(a) - piy*sin(a)))*(2*nix*(piy*cos(a) + pix*sin(a)) - 2*niy*(pix*cos(a) - piy*sin(a)));
        d2J_dxda = -nix*(2*nix*(piy*cos(a) + pix*sin(a)) - 2*niy*(pix*cos(a) - piy*sin(a)));
        d2J_dadx = -2*nix*(nix*(piy*cos(a) + pix*sin(a)) - niy*(pix*cos(a) - piy*sin(a)));
        d2J_dyda = -niy*(2*nix*(piy*cos(a) + pix*sin(a)) - 2*niy*(pix*cos(a) - piy*sin(a)));
        d2J_dady = -2*niy*(nix*(piy*cos(a) + pix*sin(a)) - niy*(pix*cos(a) - piy*sin(a)));

        // Finally reconstruct the matrix
        Matrix d2J_dX2_temp(3,3);
        d2J_dX2_temp << d2J_dx2,     d2J_dydx,	  d2J_dadx,
                d2J_dxdy,    d2J_dy2,	  d2J_dady,
                d2J_dxda,    d2J_dyda,    d2J_da2;
        d2J_dX2 += d2J_dX2_temp;

        /***********************************************************************
        d2J_dZdX -- X is the [R|T] in the form of [x, y, a]
        x, y is the translation part
        a is the rotation part in Euler format (yaw)
        *************************************************************************/

        double 	d2J_dpix_dx,    d2J_dpiy_dx,	d2J_dqix_dx,    d2J_dqiy_dx,
                d2J_dpix_dy,    d2J_dpiy_dy,	d2J_dqix_dy,    d2J_dqiy_dy,
                d2J_dpix_da,    d2J_dpiy_da,    d2J_dqix_da,    d2J_dqiy_da;

        d2J_dpix_dx = 2*nix*(nix*cos(a) + niy*sin(a));
        d2J_dpix_dy = 2*niy*(nix*cos(a) + niy*sin(a));
        d2J_dpix_da = - (2*nix*(piy*cos(a) + pix*sin(a)) - 2*niy*(pix*cos(a) - piy*sin(a)))*(nix*cos(a) + niy*sin(a)) - (2*niy*cos(a) - 2*nix*sin(a))*(nix*(qix - x - pix*cos(a) + piy*sin(a)) - niy*(y - qiy + piy*cos(a) + pix*sin(a)));
        d2J_dpiy_dx = 2*nix*(niy*cos(a) - nix*sin(a));
        d2J_dpiy_dy = 2*niy*(niy*cos(a) - nix*sin(a));
        d2J_dpiy_da = (2*nix*cos(a) + 2*niy*sin(a))*(nix*(qix - x - pix*cos(a) + piy*sin(a)) - niy*(y - qiy + piy*cos(a) + pix*sin(a))) - (2*nix*(piy*cos(a) + pix*sin(a)) - 2*niy*(pix*cos(a) - piy*sin(a)))*(niy*cos(a) - nix*sin(a));
        d2J_dqix_dx = -2*pow(nix,2);
        d2J_dqix_dy = -2*nix*niy;
        d2J_dqix_da = nix*(2*nix*(piy*cos(a) + pix*sin(a)) - 2*niy*(pix*cos(a) - piy*sin(a)));
        d2J_dqiy_dx = 2*nix*niy;
        d2J_dqiy_dy = -2*pow(niy,2);
        d2J_dqiy_da = niy*(2*nix*(piy*cos(a) + pix*sin(a)) - 2*niy*(pix*cos(a) - piy*sin(a)));

        // Finally reconstruct the matrix
        Matrix d2J_dZdX_temp(3,4);
        d2J_dZdX_temp <<  d2J_dpix_dx,    d2J_dpiy_dx,	d2J_dqix_dx,    d2J_dqiy_dx,
                d2J_dpix_dy,    d2J_dpiy_dy,	d2J_dqix_dy,    d2J_dqiy_dy,
                d2J_dpix_da,    d2J_dpiy_da,    d2J_dqix_da,    d2J_dqiy_da;
        d2J_dZdX.block(0,4*i,3,4) = d2J_dZdX_temp;

    }

    // Inverse of our state Hessian
    Matrix d2J_dX2_inv = d2J_dX2.inverse();

    // Calculate covariance
    Matrix covariance = d2J_dZdX * d2J_dZdX.transpose();
    covariance = d2J_dX2_inv * covariance * d2J_dX2_inv;

    // Multiply with our measurement noise
    return (sensorStdDev * sensorStdDev) * covariance;
}


template<typename T>
typename PointMatcher<T>::Matrix PointToPlaneWithCovErrorMinimizer<T>::getCovariance() const
{
  return covMatrix;
}

template struct PointToPlaneWithCovErrorMinimizer<float>;
template struct PointToPlaneWithCovErrorMinimizer<double>;
