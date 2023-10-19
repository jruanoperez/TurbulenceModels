/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2022 Upstream CFD GmbH
    Copyright (C) 2015-2022 OpenCFD Ltd.
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "mySpalartAllmarasDDES.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace LESModels
{

// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

template<class BasicTurbulenceModel>
tmp<volScalarField> mySpalartAllmarasDDES<BasicTurbulenceModel>::fD
(
    const volScalarField& rd
) const
{
    return
        1  - tanh(pow(this->Cd1_*rd, Cd2_));
}

template<class BasicTurbulenceModel>
tmp<volScalarField> mySpalartAllmarasDDES<BasicTurbulenceModel>::fd
(
    const volScalarField& magGradU
) const
{

    const volScalarField rval(this->r(this->nuEff(), magGradU,this->y_));
    
    if (useDeck_){
    
        const dimensionedScalar gradySMALL(this->y_.dimensions()/dimLength, SMALL);
        const dimensionedScalar gradNuTildaSMALL(this->nuTilda_.dimensions()/dimLength, SMALL);
        const dimensionedScalar magGradUSMALL(magGradU.dimensions(), SMALL);
        const dimensionedScalar magGradUDenomSMALL(dimensionSet(0,0,-3,0,0,0,0), SMALL);
	
	const volScalarField vorticityMag(mag(fvc::curl(this->U_)));

        const volVectorField n(fvc::grad(this->y_)/max(gradySMALL, mag(fvc::grad(this->y_))));
	const volScalarField Gr
        (
            C3_*max(gradNuTildaSMALL, -1.0*fvc::grad(this->nuTilda_)&n )/(max(magGradUSMALL, magGradU)*this->kappa_*this->y_)
        );

        const volScalarField GOmega
        (
            (fvc::grad(vorticityMag) & n)*sqrt(mag(this->nuTilda_)/max(magGradUDenomSMALL, pow(magGradU, 3)))
        );

        const dimensionedScalar invGOmega(GOmega.dimensions(), 1.0/200.0);

        const volScalarField alphaExp(max(-1.0,min(1.0,7.0-GOmega/invGOmega))); //està acotada entre -1 i 1
        const volScalarField exponent(max(-30.0,min(30.0,-6.0*alphaExp/max(1e-14, 1.0-pow(alphaExp,2.0)))));

	const volScalarField normVortGrad (fvc::grad(vorticityMag) & n);

        Info << "Deck shielding applied" << endl;
	
        return fD(rval)-(fD(rval)-fD(Gr)*fD(Beta_*rval))*1.0/(1.0+exp(exponent));
    }

    else return fD(rval);
}

template<class BasicTurbulenceModel>
tmp<volScalarField> mySpalartAllmarasDDES<BasicTurbulenceModel>::calc_QG
(
        const volTensorField& gradU
) const
{
        return this->calc_second_invariant(gradU);
}

template<class BasicTurbulenceModel>
tmp<volScalarField> mySpalartAllmarasDDES<BasicTurbulenceModel>::calc_QS
(
        const volTensorField& gradU
) const
{
        return this->calc_second_invariant(symm(gradU));
}

template<class BasicTurbulenceModel>
tmp<volScalarField> mySpalartAllmarasDDES<BasicTurbulenceModel>::calc_RG
(
        const volTensorField& gradU
) const
{
        return this->calc_third_invariant(gradU);
}

template<class BasicTurbulenceModel>
tmp<volScalarField> mySpalartAllmarasDDES<BasicTurbulenceModel>::calc_RS
(
        const volTensorField& gradU
) const
{
        return this->calc_third_invariant(symm(gradU));
}

template<class BasicTurbulenceModel>
tmp<volScalarField> mySpalartAllmarasDDES<BasicTurbulenceModel>::calc_V2
(
        const volTensorField& gradU
) const
{
        volScalarField trS202 = tr((symm(gradU) & symm(gradU)) & (skew(gradU) & skew(gradU)));
        volScalarField QS = this->calc_QS(gradU);
        volScalarField QO = this->calc_QG(gradU)-QS;

        return  (4.0*(trS202-2.0*QS*QO));
}

template<class BasicTurbulenceModel>
tmp<volScalarField> mySpalartAllmarasDDES<BasicTurbulenceModel>::calc_first_invariant
(
        const volTensorField& gradU
) const
{
        return tr(gradU);
}

template<class BasicTurbulenceModel>
tmp<volScalarField> mySpalartAllmarasDDES<BasicTurbulenceModel>::calc_second_invariant
(
        const volTensorField& gradU
) const
{
        return (0.5*(sqr(tr(gradU))-tr(gradU & gradU)));
}

template<class BasicTurbulenceModel>
tmp<volScalarField> mySpalartAllmarasDDES<BasicTurbulenceModel>::calc_second_invariant
(
        const volSymmTensorField& gradU
) const
{
        return (0.5*(sqr(tr(gradU))-tr(gradU & gradU)));
}

template<class BasicTurbulenceModel>
tmp<volScalarField> mySpalartAllmarasDDES<BasicTurbulenceModel>::calc_third_invariant
(
        const volTensorField& gradU
) const
{
        return det(gradU);
}

template<class BasicTurbulenceModel>
tmp<volScalarField> mySpalartAllmarasDDES<BasicTurbulenceModel>::calc_third_invariant
(
        const volSymmTensorField& gradU
) const
{
        return det(gradU);
}

template<class BasicTurbulenceModel>
vector mySpalartAllmarasDDES<BasicTurbulenceModel>::solve_cubic_eq
(
        scalar a,
        scalar b,
        scalar c,
        scalar d
) const
{
        vector x;

        scalar A,B,C;
    scalar R,Q;
        scalar sol1,sol0;
    scalar sqrtQ;
    scalar theta;

    /*normalitzem*/
    A=b/a;
    B=c/a;
    C=d/a;

    Q=(A*A-3.0*B)/9.0;
    R=(2.0*A*A*A - 9.0*A*B + 27*C)/54.0;

    /*In this case, we assume that A=0.0 and instead of solving the cubic equation
      x^3+Ax^2+Bx+C=0
      we solve the quadratic equation
      x^2+Ax+B=0*/
        if(fabs(C/A)<1e-11) {     /*before the threshold was 1e-13*/ /*XAVI95.x3*/
                sol1=0.5*(-A+sqrt(A*A-4.0*B));
                sol0=0.5*(-A-sqrt(A*A-4.0*B));
                if(sol0>0.0) {
                        x[0]=sol1; x[1]=sol0; x[2]=0.0;
                        return x;
                }

                if(sol1<0.0) {
                        x[0]=0.0; x[1]=sol1; x[2]=sol0;
                        return x;
                }

                x[0]=sol1; x[1]=0.0; x[2]=sol0;

                return x;
        }

        if(R*R>=Q*Q*Q)
                FatalErrorInFunction << "solve_cubic_eq\t does not have real roots  R="<<R<<" Q="<<Q<<" R^2-Q^3="<<R*R-Q*Q*Q<<" a="<<a<<" b="<<b<<" c="<<c<<" d="<<d<<" theta="<<acos(R/pow(Q,1.5)) << exit(FatalError) ;

        sqrtQ=Foam::sqrt(Q);

        theta=acos(R/pow(sqrtQ,3.0));

        //In this way the roots are already ordered x[0]>=x[1]>=x[2]
        x[2]=-2.0*sqrtQ*cos(theta/3.0      )-A/3.0;
        x[0]=-2.0*sqrtQ*cos((theta+6.28318530717958647692)/3.0)-A/3.0;
        x[1]=-2.0*sqrtQ*cos((theta-6.28318530717958647692)/3.0)-A/3.0;

        return x;
}

template<class BasicTurbulenceModel>
tmp<volScalarField> mySpalartAllmarasDDES<BasicTurbulenceModel>::SMG_model
(
        const volScalarField& Q
) const
{
        forAll(Q.internalField(), celli)
        {
                if(Q[celli]>1e-15)
                        Info << "QS(-): " << Q[celli] << endl;
        }

        //const cellList& cells = this->mesh().cells();
        //forAll(cells, celli)
        //{
        //      if(Q[celli]>=0.0)
        //              Info << "Positive Value: " << Q[celli] << endl;
        //}

        return Foam::sqrt(-4.0*Q);
}

template<class BasicTurbulenceModel>
tmp<volScalarField> mySpalartAllmarasDDES<BasicTurbulenceModel>::WALE_model
(
        const volScalarField& QS,
        const volScalarField& QG,
        const volScalarField& V2
) const
{
        volScalarField aux = 0.5*V2+2.0/3.0*QG*QG;

        forAll(aux.internalField(), celli)
        {
                if(fabs(aux[celli])<1e-15)
                        aux[celli] *= 0.0;
                //if(V2[celli]<1e-15)
                //      Info << "V2(-): " << V2[celli] << endl;
                //if(QG[celli]>1e-15)
                //      Info << "QG(-): " << QG[celli] << endl;
                //if(aux[celli]<1e-15)
                //      Info << "aux(-): " << aux[celli] << endl;
                //if(QS[celli]>1e-15)
                //      Info << "QS(-): " << QS[celli] << endl;
        }

        forAll(aux.boundaryField(), patchi)
        {
                fvPatchScalarField&     aux_patch = aux.boundaryFieldRef()[patchi];
                //const fvPatchScalarField&     V2_patch =  V2.boundaryFieldRef()[patchi];
                //const fvPatchScalarField&     QG_patch =  QG.boundaryFieldRef()[patchi];
                forAll(aux_patch, facei)
                {
                        if(fabs(aux_patch[facei])<1e-15)
                                aux_patch[facei] *= 0.0;
                        //if(V2_patch[facei]<1e-15)
                        //      Info << "V2_patch(-): " << V2_patch[facei] << endl;
                        //if(QG_patch[facei]<1e-15)
                        //      Info << "QG_patch(-): " << QG_patch[facei] << endl;
                        //if(aux_patch[facei]<1e-15)
                        //      Info << "aux_patch(-): " << aux_patch[facei] << endl;
                }
        }

        return  (pow(aux,3.0/2.0)/(pow(-2.0*QS,5.0/2.0)+pow(aux,5.0/4.0)));
}

template<class BasicTurbulenceModel>
tmp<volScalarField> mySpalartAllmarasDDES<BasicTurbulenceModel>::VREMAN_model
(
        const volScalarField& QS,
        const volScalarField& QG,
        const volScalarField& V2
) const
{
        volScalarField den = 2.0*(QG-2.0*QS);
        volScalarField diffOperator = (QG*QG+V2)/den;

        forAll(diffOperator.internalField(), celli)
        {
                if(fabs(QG[celli])<=1e-15 )
                        diffOperator[celli] *= 0.0;
        }
        forAll(diffOperator.boundaryField(), patchi)
        {
                fvPatchScalarField&     aux_patch = diffOperator.boundaryFieldRef()[patchi];
                forAll(aux_patch, facei)
                {
                        if(fabs(aux_patch[facei])<=1e-15 )
                                aux_patch[facei] *= 0.0;
                }
        }

        return  Foam::sqrt(diffOperator);
}

template<class BasicTurbulenceModel>
tmp<volScalarField> mySpalartAllmarasDDES<BasicTurbulenceModel>::SIGMA_model
(
        const volScalarField& QS,
        const volScalarField& QG,
        const volScalarField& RG,
        const volScalarField& V2
) const
{
        int checking = 1;

        volVectorField eigenval_scaled(
                IOobject
                (
                        "eigenval_scaled",
                        this->time().timeName(),
                        this->mesh(),
                        IOobject::NO_READ,
                        IOobject::NO_WRITE
                ),
                this->mesh(),
        dimensionedVector("zero", dimless, Zero)
        );

        volScalarField P = 2.0*(QG-2.0*QS);
        volScalarField Q = V2+QG*QG;
        volScalarField R = RG*RG;

        volScalarField scaling=1.0/P;

        volScalarField Pscaled=scaling*P;
        volScalarField Qscaled=scaling*scaling*Q;
        volScalarField Rscaled=scaling*scaling*scaling*R;
        forAll(Pscaled.internalField(), celli)
        {
                /*solving the cubic equation with the invariants of scaling*GG^t;
                  the results are the ordered eigenvalues of the scaling*GG^t tensor*/
                eigenval_scaled[celli] = solve_cubic_eq(1.0,-Pscaled[celli],Qscaled[celli],-Rscaled[celli]);
        }

        forAll(Pscaled.boundaryField(), patchi)
        {
                fvPatchVectorField&     aux_patch_ev = eigenval_scaled.boundaryFieldRef()[patchi];
                fvPatchScalarField&     aux_patch_p = Pscaled.boundaryFieldRef()[patchi];
                fvPatchScalarField&     aux_patch_q = Qscaled.boundaryFieldRef()[patchi];
                fvPatchScalarField&     aux_patch_r = Rscaled.boundaryFieldRef()[patchi];

                forAll(aux_patch_p, facei)
                {
                        /*solving the cubic equation with the invariants of scaling*GG^t;
                        the results are the ordered eigenvalues of the scaling*GG^t tensor*/
                        aux_patch_ev[facei] = solve_cubic_eq(1.0,-aux_patch_p[facei],aux_patch_q[facei],-aux_patch_r[facei]);
                }
        }

        /*computing the eigenvalues of GG^t*/
        volVectorField eigenval=eigenval_scaled/scaling;
        volVectorField sigma = eigenval/pow(mag(eigenval),0.5);
        sigma.component(0) =Foam::sqrt(eigenval.component(0));
        sigma.component(1) =Foam::sqrt(eigenval.component(1));
        sigma.component(2) =Foam::sqrt(eigenval.component(2));
        /*checkings*/
        if(checking) {
                volScalarField errorP = P-(eigenval.component(0)+eigenval.component(1)+eigenval.component(2));
                volScalarField errorQ = Q-(eigenval.component(0)*eigenval.component(1)+eigenval.component(0)*eigenval.component(2)+eigenval.component(1)*eigenval.component(2));
                volScalarField errorR = R-(eigenval.component(0)*eigenval.component(1)*eigenval.component(2));

        forAll(sigma.internalField(), celli)
                {
                forAll(sigma[celli], i){
                                if(i<2)
                                        if(sigma[celli][i]<sigma[celli][i+1])
                                                FatalErrorInFunction << "Alguna cosa a fallat amb les sigmes!!!" << exit(FatalError);
                        }
                        if(fabs(errorP[celli]/P[celli])>1e-7 && fabs(errorP[celli])>1e-9)
                                FatalErrorInFunction << "Something wrong with P!!!" << exit(FatalError);
                        if(fabs(errorQ[celli]/P[celli])>1e-7 && fabs(errorQ[celli])>1e-7)
                                FatalErrorInFunction << "Something wrong with Q!!!" << exit(FatalError);
                        if(fabs(errorR[celli]/P[celli])>1e-7 && fabs(errorR[celli])>1e-5)
                                FatalErrorInFunction << "Something wrong with R!!!" << exit(FatalError);
                }

                forAll(sigma.boundaryField(), patchi)
                {
                        fvPatchVectorField&     aux_patch_sigma         = sigma.boundaryFieldRef()[patchi];
                        fvPatchScalarField&     aux_patch_p             = P.boundaryFieldRef()[patchi];
                        fvPatchScalarField&     aux_patch_errorp        = errorP.boundaryFieldRef()[patchi];
                        fvPatchScalarField&     aux_patch_errorq        = errorQ.boundaryFieldRef()[patchi];
                        fvPatchScalarField&     aux_patch_errorr        = errorR.boundaryFieldRef()[patchi];

                        forAll(aux_patch_sigma, facei)
                        {
                        forAll(aux_patch_sigma[facei], i){
                                        if(i<2)
                                                if(aux_patch_sigma[facei][i]<aux_patch_sigma[facei][i+1])
                                                        FatalErrorInFunction << "Alguna cosa a fallat amb les sigmes (BOUNDARY)!!!" << exit(FatalError);
                                }
                                if(fabs(aux_patch_errorp[facei]/aux_patch_p[facei])>1e-7 && fabs(aux_patch_errorp[facei])>1e-9)
                                        FatalErrorInFunction << "Something wrong with P!!! (BOUNDARY)" << exit(FatalError);
                                if(fabs(aux_patch_errorq[facei]/aux_patch_p[facei])>1e-7 && fabs(aux_patch_errorq[facei])>1e-7)
                                        FatalErrorInFunction << "Something wrong with Q!!! (BOUNDARY)" << exit(FatalError);
                                if(fabs(aux_patch_errorr[facei]/aux_patch_p[facei])>1e-7 && fabs(aux_patch_errorr[facei])>1e-5)
                                        FatalErrorInFunction << "Something wrong with R!!! (BOUNDARY)" << exit(FatalError);
                        }
                }
        }

        return sigma.component(2)*(sigma.component(0)-sigma.component(1))*(sigma.component(1)-sigma.component(2))/(sigma.component(0)*sigma.component(0));
}

template<class BasicTurbulenceModel>
tmp<volScalarField> mySpalartAllmarasDDES<BasicTurbulenceModel>::SIGMA_Alex_model
(
        const volScalarField& QS,
        const volScalarField& QG,
        const volScalarField& RG,
        const volScalarField& V2
) const
{
        //int checking = 1;

        scalar C1_3 = 0.333333333333333333333333333333;
        scalar C1_6 = 0.166666666666666666666666666666;
        scalar C1_9 = 0.111111111111111111111111111111;
        scalar C1_27 = C1_9 * C1_3;

        scalar tiny = 1e-16;
        scalar tinyflt = 1e-8;

        volScalarField filter(
                IOobject
                (
                        "filter",
                        this->time().timeName(),
                        this->mesh(),
                        IOobject::NO_READ,
                        IOobject::NO_WRITE
                ),
                this->mesh(),
        dimensionedScalar("zero", dimless, Zero)
        );
        filter = 1.0;
        volScalarField P = 2.0*(QG-2.0*QS);
        volScalarField Q = V2+QG*QG;
        volScalarField R = RG*RG;

        volScalarField alpha1 = C1_9 * P*P - C1_3 * Q;
        volScalarField alpha2 = C1_27 * P*P*P- C1_6 * P*Q + 0.5*R;

        forAll(alpha1.internalField(), celli)
        {
                if( alpha1[celli] < tinyflt ) {
                        alpha1[celli] = 1.0;
                        filter[celli] = 0.0;
                }
        }

        forAll(alpha1.boundaryField(), patchi)
        {
                fvPatchScalarField&     aux_patch_alpha1 = alpha1.boundaryFieldRef()[patchi];
                fvPatchScalarField&     aux_patch_filter = filter.boundaryFieldRef()[patchi];

                forAll(aux_patch_alpha1,facei)
                {
                        /*solving the cubic equation with the invariants of scaling*GG^t;
                        the results are the ordered eigenvalues of the scaling*GG^t tensor*/
                        if( aux_patch_alpha1[facei] < tinyflt ) {
                                aux_patch_alpha1[facei] = 1.0;
                                aux_patch_filter[facei] = 0.0;
                        }
                }
        }
    volScalarField alpha1_sqrt = Foam::sqrt(alpha1);
        volScalarField alpha3      = alpha2 / (alpha1*alpha1_sqrt);

        forAll(alpha3.internalField(), celli)
        {
                if( alpha3[celli] < -1 || alpha3[celli] > 1 ) {
                        alpha3[celli] = 0.0;
                        filter[celli] = 0.0;
                } else
                        alpha3[celli] = C1_3 * acos(alpha3[celli]);
        }

        forAll(alpha3.boundaryField(), patchi)
        {
                fvPatchScalarField&     aux_patch_alpha3 = alpha3.boundaryFieldRef()[patchi];
                fvPatchScalarField&     aux_patch_filter = filter.boundaryFieldRef()[patchi];

                forAll(aux_patch_alpha3,facei)
                {
                        /*solving the cubic equation with the invariants of scaling*GG^t;
                        the results are the ordered eigenvalues of the scaling*GG^t tensor*/
                        if( aux_patch_alpha3[facei] < -1 || aux_patch_alpha3[facei] > 1 ) {
                                aux_patch_alpha3[facei] = 0.0;
                                aux_patch_filter[facei] = 0.0;
                        }
                }
        }

        volVectorField sigma(
                IOobject
                (
                        "sigma",
                        this->time().timeName(),
                        this->mesh(),
                        IOobject::NO_READ,
                        IOobject::NO_WRITE
                ),
                this->mesh(),
        dimensionedVector("zero", dimensionSet(0, 0, -1, 0, 0, 0 ,0), Zero)
        );

        forAll(sigma,celli)
        {
                sigma[celli][0] = C1_3*P[celli] + 2.0*alpha1_sqrt[celli]*cos(alpha3[celli]);
                sigma[celli][1] = C1_3*P[celli] - 2.0*alpha1_sqrt[celli]*cos(C1_3*M_PI + alpha3[celli]);
                sigma[celli][2] = C1_3*P[celli] - 2.0*alpha1_sqrt[celli]*cos(C1_3*M_PI - alpha3[celli]);
        }

        forAll(sigma.internalField(), celli)
        {
                if( sigma[celli][0] < tiny  || sigma[celli][1] < tiny || sigma[celli][2] < tiny || !(sigma[celli][0] > sigma[celli][1] && sigma[celli][1] > sigma[celli][2])){
                        //Info << "SHIT[" << filter[celli] << "]: " << alpha1[celli] << "\t"<< alpha3[celli] << "\t" << sigma[celli][0] << "\t" << sigma[celli][1] << "\t" << sigma[celli][2] << endl;
                        sigma[celli] = vector(1.0,0.0,0.0);
                        //sigma[celli][0] =  0.0*sigma[celli][0] + 1.0;
                        //sigma[celli][1] =  0.0*sigma[celli][1];
                        //sigma[celli][2] =  0.0*sigma[celli][2];
                        filter[celli] = 0.0;
                }
    }

        forAll(sigma.boundaryField(), patchi)
        {
                fvPatchVectorField&     aux_patch_sigma  = sigma.boundaryFieldRef()[patchi];
                fvPatchScalarField&     aux_patch_filter = filter.boundaryFieldRef()[patchi];

                forAll(aux_patch_sigma,facei)
                {
                        /*solving the cubic equation with the invariants of scaling*GG^t;
                        the results are the ordered eigenvalues of the scaling*GG^t tensor*/
                        if( aux_patch_sigma[facei][0] < tiny  || aux_patch_sigma[facei][1] < tiny || aux_patch_sigma[facei][2] < tiny || !(aux_patch_sigma[facei][0] > aux_patch_sigma[facei][1] && aux_patch_sigma[facei][1] > aux_patch_sigma[facei][2])){
                                aux_patch_sigma[facei] = vector(1.0,0.0,0.0);
                                //aux_patch_sigma[facei][0]  = 0.0*aux_patch_sigma[facei][0] + 1.0;
                                //aux_patch_sigma[facei][1]  = 0.0*aux_patch_sigma[facei][1];
                                //aux_patch_sigma[facei][2]  = 0.0*aux_patch_sigma[facei][2];
                                aux_patch_filter[facei] = 0.0;
                        }
                }
        }

        sigma.component(0) = Foam::sqrt(sigma.component(0));
        sigma.component(1) = Foam::sqrt(sigma.component(1));
        sigma.component(2) = Foam::sqrt(sigma.component(2));

        return filter*sigma.component(2)*(sigma.component(0)-sigma.component(1))*(sigma.component(1)-sigma.component(2))/(sigma.component(0)*sigma.component(0));
        //return sigma.component(2)*(sigma.component(0)-sigma.component(1))*(sigma.component(1)-sigma.component(2))/(sigma.component(0)*sigma.component(0));
}

template<class BasicTurbulenceModel>
tmp<volScalarField> mySpalartAllmarasDDES<BasicTurbulenceModel>::SIGMA_Fuchs_model
(
        const volScalarField& QS,
        const volScalarField& QG,
        const volScalarField& RG,
        const volScalarField& V2
) const
{
        volScalarField P = 2.0*(QG-2.0*QS);
        volScalarField Q = V2+QG*QG;
        volScalarField R = RG*RG;

        const volScalarField alpha1 = (sqr(P)/9.0)-(Q/3.0);
        const volScalarField alpha2 = (pow(P, 3.0)/27.0)-(P*Q/6.0)+(R/2.0);
        const volScalarField alpha3 = (1.0/3.0)*acos(
                        max(
                                scalar(-1) + dimensionedScalar("SMALL", dimless, SMALL),
                                min(
                                        scalar(1) - dimensionedScalar("SMALL", dimless, SMALL),
                                        alpha2/pow(alpha1, 3.0/2.0)
                                   )
                           ));

        const volScalarField sigma1 = sqrt(max(((P/3.0)+(2.0*sqrt(alpha1)*cos(alpha3))), dimensionedScalar("SMALL", dimless/sqr(dimTime), SMALL)));
        const volScalarField sigma2 = sqrt(max(((P/3.0)-(2.0*sqrt(alpha1)*cos((constant::mathematical::pi/3.0)+alpha3))), dimensionedScalar("SMALL", dimless/sqr(dimTime), SMALL)));
        const volScalarField sigma3 = sqrt(max(((P/3.0)-(2.0*sqrt(alpha1)*cos((constant::mathematical::pi/3.0)-alpha3))), dimensionedScalar("SMALL", dimless/sqr(dimTime), SMALL)));

        return sigma3*(sigma1-sigma2)*(sigma2-sigma3)/max(sqr(sigma1), dimensionedScalar("SMALL", dimless/sqr(dimTime), SMALL));
}

template<class BasicTurbulenceModel>
tmp<volScalarField> mySpalartAllmarasDDES<BasicTurbulenceModel>::S3PQR_model
(
        const volScalarField& QS,
        const volScalarField& QG,
        const volScalarField& RG,
        const volScalarField& V2,
        const scalar p,
        const scalar q,
        const scalar r,
        const bool p_sqrt2,
        const bool q_sqrt2,
        const bool r_sqrt2
) const
{
        volScalarField P = 2.0*(QG-2.0*QS);
        volScalarField Q = V2+QG*QG;
        volScalarField R = RG*RG;

        forAll(QG.internalField(), celli)
        {
                if( fabs(QG[celli]) < 1e-15  ){
                        P[celli] *= 0.0;
                        Q[celli] *= 0.0;
                        R[celli] *= 0.0;
                        Info << "RECTIFYING: "<< P[celli] << Q[celli] << R[celli] << endl ;
                }
        }
        //forAll(P.boundaryField(), patchi)
        //{
        //      fvPatchScalarField      aux_patch_qg = QG.boundaryFieldRef()[patchi];
        //      fvPatchScalarField&     aux_patch_p = P.boundaryFieldRef()[patchi];
        //      fvPatchScalarField&     aux_patch_q = Q.boundaryFieldRef()[patchi];
        //      fvPatchScalarField&     aux_patch_r = R.boundaryFieldRef()[patchi];

        //      forAll(aux_patch_qg, facei)
        //      {
        //              if( fabs(aux_patch_qg[facei]) < 1e-15 ){
        //                      aux_patch_p[facei] *= 0.0;
        //                      aux_patch_q[facei] *= 0.0;
        //                      aux_patch_r[facei] *= 0.0;
        //              }
        //      }
        //}

        return pow(P,p)*pow(Q,q)*pow(R,r);
}

template<class BasicTurbulenceModel>
tmp<volScalarField> mySpalartAllmarasDDES<BasicTurbulenceModel>::S3PQ_model
(
        const volScalarField& QS,
        const volScalarField& QG,
        const volScalarField& RG,
        const volScalarField& V2
) const
{
        //scalar p=-5.0/2.0, q=3.0/2.0, r=0.0;
        scalar p=-5.0/2.0, q=3.0/2.0, r=0.0;
        bool   p_sqrt2 = true, q_sqrt2 = true, r_sqrt2 = false;

        return S3PQR_model(QS,QG,RG,V2,p,q,r,p_sqrt2, q_sqrt2, r_sqrt2 );
}

template<class BasicTurbulenceModel>
tmp<volScalarField> mySpalartAllmarasDDES<BasicTurbulenceModel>::S3PR_model
(
        const volScalarField& QS,
        const volScalarField& QG,
        const volScalarField& RG,
        const volScalarField& V2
) const
{
        scalar p=-1.0, q=0.0, r=1.0/2.0;
        bool   p_sqrt2 = false, q_sqrt2 = false, r_sqrt2 = true;

        return S3PQR_model(QS,QG,RG,V2,p,q,r,p_sqrt2, q_sqrt2, r_sqrt2 );
}

template<class BasicTurbulenceModel>
tmp<volScalarField> mySpalartAllmarasDDES<BasicTurbulenceModel>::S3QR_model
(
        const volScalarField& QS,
        const volScalarField& QG,
        const volScalarField& RG,
        const volScalarField& V2
) const
{
        scalar p=0.0, q=-1.0, r=5.0/6.0;
        bool   p_sqrt2 = false, q_sqrt2 = false, r_sqrt2 = true;

        return S3PQR_model(QS,QG,RG,V2,p,q,r,p_sqrt2, q_sqrt2, r_sqrt2 );
}

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

template<class BasicTurbulenceModel>
tmp<volScalarField> mySpalartAllmarasDDES<BasicTurbulenceModel>::Stilda
(
    const volScalarField& chi,
    const volScalarField& fv1,
    const volTensorField& gradU,
    const volScalarField& dTilda
) const
{
        const volScalarField& lRAS(this->y_);
        const volScalarField fv2(this->fv2(chi, fv1));
        const volScalarField lLES(this->lengthScaleLES(chi, fv1));
        const volScalarField Omega(this->Omega(gradU));

	volScalarField S_LES(
                IOobject
                (
                        "S_LES",
                        this->time().timeName(),
                        this->mesh(),
                        IOobject::NO_READ,
                        IOobject::NO_WRITE,
                        false
                ),
                Omega
        );

        Info << "Assessing Invariants" << endl;
        scalar B_LES_v = B_LES.value();
        const volTensorField gradUdev = dev(gradU);

        if(LES_MODEL=="SMG"){
                if(B_LES_v==-1) B_LES_v = 0.7225; // Validated for DHIT mesh64x64x64

		Info << "I2" << endl;
		volScalarField QS(this->calc_QS(gradUdev));

		S_LES = this->SMG_model(QS);
        } else if(LES_MODEL=="WALE"){
                if(B_LES_v==-1) B_LES_v = 5.5;  // Validated for DHIT mesh64x64x64

        	Info << "I1" << endl;
        	volScalarField QG(this->calc_QG(gradUdev));
                Info << "I2" << endl;
                volScalarField QS(this->calc_QS(gradUdev));
	        Info << "I5" << endl;
		volScalarField V2(this->calc_V2(gradUdev));


		S_LES = this->WALE_model(QS,QG,V2);
        } else if(LES_MODEL=="VREMAN"){
                if(B_LES_v==-1) B_LES_v = 1.75; // Validated for DHIT mesh64x64x64

                Info << "I1" << endl;
                volScalarField QG(this->calc_QG(gradUdev));
                Info << "I2" << endl;
                volScalarField QS(this->calc_QS(gradUdev));
                Info << "I5" << endl;
                volScalarField V2(this->calc_V2(gradUdev));

		S_LES = this->VREMAN_model(QS,QG,V2);
        } else if(LES_MODEL=="SIGMA"){
                if(B_LES_v==-1) B_LES_v = 270; // Validated for DHIT mesh64x64x64

                Info << "I1" << endl;
                volScalarField QG(this->calc_QG(gradUdev));
                Info << "I2" << endl;
                volScalarField QS(this->calc_QS(gradUdev));
	        Info << "I3" << endl;
		volScalarField RG(this->calc_RG(gradUdev));
		Info << "I5" << endl;
                volScalarField V2(this->calc_V2(gradUdev));

		S_LES = this->SIGMA_model(QS,QG,RG,V2);
        } else if(LES_MODEL=="SIGMA_Alex"){
                if(B_LES_v==-1) B_LES_v = 65;
         
                Info << "I1" << endl;
                volScalarField QG(this->calc_QG(gradUdev));
                Info << "I2" << endl;
                volScalarField QS(this->calc_QS(gradUdev));
                Info << "I3" << endl;
                volScalarField RG(this->calc_RG(gradUdev));
                Info << "I5" << endl;
                volScalarField V2(this->calc_V2(gradUdev));

		S_LES = this->SIGMA_Alex_model(QS,QG,RG,V2); // Validated for DHIT mesh64x64x64
        } else if(LES_MODEL=="SIGMA_Fuchs"){
                if(B_LES_v==-1) B_LES_v = 70;
         
                Info << "I1" << endl;
                volScalarField QG(this->calc_QG(gradUdev));
                Info << "I2" << endl;
                volScalarField QS(this->calc_QS(gradUdev));
                Info << "I3" << endl;
                volScalarField RG(this->calc_RG(gradUdev));
                Info << "I5" << endl;
                volScalarField V2(this->calc_V2(gradUdev));

                S_LES = this->SIGMA_Fuchs_model(QS,QG,RG,V2);
        } else if(LES_MODEL=="S3PQ"){
                if(B_LES_v==-1) B_LES_v = 10.5; // Validated for DHIT mesh64x64x64
         
                Info << "I1" << endl;
                volScalarField QG(this->calc_QG(gradUdev));
                Info << "I2" << endl;
                volScalarField QS(this->calc_QS(gradUdev));
                Info << "I3" << endl;
                volScalarField RG(this->calc_RG(gradUdev));
                Info << "I5" << endl;
                volScalarField V2(this->calc_V2(gradUdev));

                S_LES = this->S3PQ_model(QS,QG,RG,V2);
        } else if(LES_MODEL=="S3PR"){
                if(B_LES_v==-1) B_LES_v = 17.0; // Validated for DHIT mesh64x64x64
         
                Info << "I1" << endl;
                volScalarField QG(this->calc_QG(gradUdev));
                Info << "I2" << endl;
                volScalarField QS(this->calc_QS(gradUdev));
                Info << "I3" << endl;
                volScalarField RG(this->calc_RG(gradUdev));
                Info << "I5" << endl;
                volScalarField V2(this->calc_V2(gradUdev));

                S_LES = this->S3PR_model(QS,QG,RG,V2);
        } else if(LES_MODEL=="S3QR"){
                if(B_LES_v==-1) B_LES_v = 20.0; // Validated for DHIT mesh64x64x64
         
                Info << "I1" << endl;
                volScalarField QG(this->calc_QG(gradUdev));
                Info << "I2" << endl;
                volScalarField QS(this->calc_QS(gradUdev));
                Info << "I3" << endl;
                volScalarField RG(this->calc_RG(gradUdev));
                Info << "I5" << endl;
                volScalarField V2(this->calc_V2(gradUdev));

                S_LES = this->S3QR_model(QS,QG,RG,V2);
        } else {
                FatalErrorInFunction
                        << "LES_MODEL[WTF?]: " << LES_MODEL
                        << exit(FatalError);
        }

        Info << "Invariants assessed" << endl;
        Info << "Differential Operator Assessed" << endl;

        const volScalarField SDDES
        (
            Omega - fd(mag(gradU))*pos(lRAS - lLES)*(Omega - B_LES_v*S_LES)
        );

        return
            max
            (
                SDDES + fv2*this->nuTilda_/sqr(this->kappa_*dTilda),
                this->Cs_*SDDES
            );

}


template<class BasicTurbulenceModel>
tmp<volScalarField> mySpalartAllmarasDDES<BasicTurbulenceModel>::dTilda
(
    const volScalarField& chi,
    const volScalarField& fv1,
    const volTensorField& gradU
) const
{
    const volScalarField& lRAS(this->y_);
    const volScalarField lLES(this->lengthScaleLES(chi, fv1));
    const dimensionedScalar l0(dimLength, Zero);

    return max
    (
        lRAS - fd(mag(gradU))*max(lRAS - lLES, l0),
        dimensionedScalar("small", dimLength, SMALL)
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class BasicTurbulenceModel>
mySpalartAllmarasDDES<BasicTurbulenceModel>::mySpalartAllmarasDDES
(
    const alphaField& alpha,
    const rhoField& rho,
    const volVectorField& U,
    const surfaceScalarField& alphaRhoPhi,
    const surfaceScalarField& phi,
    const transportModel& transport,
    const word& propertiesName,
    const word& type
)
:
    SpalartAllmarasDES<BasicTurbulenceModel>
    (
        alpha,
        rho,
        U,
        alphaRhoPhi,
        phi,
        transport,
        propertiesName,
        type
    ),

    Cd1_
    (
        this->useSigma_
          ? dimensioned<scalar>::getOrAddToDict
            (
                "Cd1Sigma",
                this->coeffDict_,
                10
            )
          : dimensioned<scalar>::getOrAddToDict
            (
                "Cd1",
                this->coeffDict_,
                8
            )
    ),
    Cd2_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "Cd2",
            this->coeffDict_,
            3
        )
    ),
    C3_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "C3",
            this->coeffDict_,
            25
        )
    ),
    C4_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "C4",
            this->coeffDict_,
            0.03
        )
    ),
    Beta_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "Beta",
            this->coeffDict_,
            2.5
        )
    ),
    useDeck_
    (
    	Switch::getOrAddToDict
    	(
    	    "useDeck",
    	    this->coeffDict_,
    	    false
    	)
    ),
    LES_MODEL(
            this->coeffDict_.lookup("LES_MODEL")
    ),
    B_LES
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "B_LES",
            this->coeffDict_,
            -1
        )
    )

{
    if (type == typeName)
    {
        this->printCoeffs(type);
    }
    if(LES_MODEL != "SMG"){
        Info << "WARNING: Cd1_ from 8 to 10 [Mockett 2015]" << endl;
        Cd1_ = 10.0;
    }
    if(B_LES.value() != -1){
        Info << "WARNING: B_LES value has been modified." << endl;
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class BasicTurbulenceModel>
bool mySpalartAllmarasDDES<BasicTurbulenceModel>::read()
{
    if (SpalartAllmarasDES<BasicTurbulenceModel>::read())
    {
        Cd1_.readIfPresent(this->coeffDict());
        Cd2_.readIfPresent(this->coeffDict());
        C3_.readIfPresent(this->coeffDict());
        useDeck_.readIfPresent("useDeck", this->coeffDict());
        LES_MODEL = word(this->coeffDict().lookup("LES_MODEL"));
        B_LES.readIfPresent(this->coeffDict());

        return true;
    }

    return false;
}


template<class BasicTurbulenceModel>
tmp<volScalarField> mySpalartAllmarasDDES<BasicTurbulenceModel>::fd() const
{
    return fd(mag(fvc::grad(this->U_)));
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace LESModels
} // End namespace Foam

// ************************************************************************* //
