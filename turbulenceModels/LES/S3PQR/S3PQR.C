/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2022 Upstream CFD GmbH
    Copyright (C) 2022-2023 OpenCFD Ltd.
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

#include "S3PQR.H"
#include "fvOptions.H"
#include "DESModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace LESModels
{
// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

template<class BasicTurbulenceModel>
tmp<volScalarField> S3PQR<BasicTurbulenceModel>::calc_QG
(
        const volTensorField& gradU
) const
{
        return this->calc_second_invariant(gradU);
}

template<class BasicTurbulenceModel>
tmp<volScalarField> S3PQR<BasicTurbulenceModel>::calc_QS
(
        const volTensorField& gradU
) const
{
        return this->calc_second_invariant(symm(gradU));
}

template<class BasicTurbulenceModel>
tmp<volScalarField> S3PQR<BasicTurbulenceModel>::calc_RG
(
        const volTensorField& gradU
) const
{
        return this->calc_third_invariant(gradU);
}

template<class BasicTurbulenceModel>
tmp<volScalarField> S3PQR<BasicTurbulenceModel>::calc_RS
(
        const volTensorField& gradU
) const
{
        return this->calc_third_invariant(symm(gradU));
}

template<class BasicTurbulenceModel>
tmp<volScalarField> S3PQR<BasicTurbulenceModel>::calc_V2
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
tmp<volScalarField> S3PQR<BasicTurbulenceModel>::calc_first_invariant
(
        const volTensorField& gradU
) const
{
        return tr(gradU);
}

template<class BasicTurbulenceModel>
tmp<volScalarField> S3PQR<BasicTurbulenceModel>::calc_second_invariant
(
        const volTensorField& gradU
) const
{
        return (0.5*(sqr(tr(gradU))-tr(gradU & gradU)));
}

template<class BasicTurbulenceModel>
tmp<volScalarField> S3PQR<BasicTurbulenceModel>::calc_second_invariant
(
        const volSymmTensorField& gradU
) const
{
        return (0.5*(sqr(tr(gradU))-tr(gradU & gradU)));
}

template<class BasicTurbulenceModel>
tmp<volScalarField> S3PQR<BasicTurbulenceModel>::calc_third_invariant
(
        const volTensorField& gradU
) const
{
        return det(gradU);
}

template<class BasicTurbulenceModel>
tmp<volScalarField> S3PQR<BasicTurbulenceModel>::calc_third_invariant
(
        const volSymmTensorField& gradU
) const
{
        return det(gradU);
}

template<class BasicTurbulenceModel>
vector S3PQR<BasicTurbulenceModel>::solve_cubic_eq
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

    //normalitzem
    A=b/a;
    B=c/a;
    C=d/a;

    Q=(A*A-3.0*B)/9.0;
    R=(2.0*A*A*A - 9.0*A*B + 27*C)/54.0;

    //In this case, we assume that A=0.0 and instead of solving the cubic equation
//      x^3+Ax^2+Bx+C=0
//      we solve the quadratic equation
//      x^2+Ax+B=0
        if(fabs(C/A)<1e-11) {     //before the threshold was 1e-13 //XAVI95.x3
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
tmp<volScalarField> S3PQR<BasicTurbulenceModel>::S3PQR_model
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
        return pow(P,p)*pow(Q,q)*pow(R,r);
}

template<class BasicTurbulenceModel>
tmp<volScalarField> S3PQR<BasicTurbulenceModel>::S3PQ_model
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
tmp<volScalarField> S3PQR<BasicTurbulenceModel>::S3PR_model
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
tmp<volScalarField> S3PQR<BasicTurbulenceModel>::S3QR_model
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
void S3PQR<BasicTurbulenceModel>::correctNut()
{

    const volTensorField gradUdev = dev(fvc::grad(this->U_));

    Info << "I1" << endl;
    volScalarField QG(this->calc_QG(gradUdev));
    Info << "I2" << endl;
    volScalarField QS(this->calc_QS(gradUdev));
    Info << "I3" << endl;
    volScalarField RG(this->calc_RG(gradUdev));
    Info << "I5" << endl;
    volScalarField V2(this->calc_V2(gradUdev));

    scalar C_S3PQR_v = C_S3PQR_.value();

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
            mag(gradUdev)
    );

        if(S3PQR_MODEL=="S3PQ"){
                if(C_S3PQR_v==-1) C_S3PQR_v = 0.572; // Validated for DHIT mesh64x64x64
                S_LES = this->S3PQ_model(QS,QG,RG,V2);
        } else if(S3PQR_MODEL=="S3PR"){
		if(C_S3PQR_v==-1) C_S3PQR_v = 0.709; // Validated for DHIT mesh64x64x64
                S_LES = this->S3PR_model(QS,QG,RG,V2);
        } else if(S3PQR_MODEL=="S3QR"){
		if(C_S3PQR_v==-1) C_S3PQR_v = 0.762; // Validated for DHIT mesh64x64x64
                S_LES = this->S3QR_model(QS,QG,RG,V2);
        } else {
                FatalErrorInFunction
                        << "S3PQR_MODEL[WTF?]: " << S3PQR_MODEL
                        << exit(FatalError);
        }


    this->nut_ =
        sqr(C_S3PQR_v*this->delta())
       *S_LES;

    this->nut_.correctBoundaryConditions();
    fv::options::New(this->mesh_).correct(this->nut_);

    BasicTurbulenceModel::correctNut();
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class BasicTurbulenceModel> S3PQR<BasicTurbulenceModel>::S3PQR
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
    LESeddyViscosity<BasicTurbulenceModel>
    (
        type,
        alpha,
        rho,
        U,
        alphaRhoPhi,
        phi,
        transport,
        propertiesName
    ),

    Ck_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Ck",
            this->coeffDict_,
            0.094
        )
    ),

    Cw_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Cw",
            this->coeffDict_,
            0.325
        )
    ),
    S3PQR_MODEL(
            this->coeffDict_.lookup("S3PQR_MODEL")
    ),
    C_S3PQR_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "CS3PQR",
            this->coeffDict_,
            -1
        )
    )
{
    if (type == typeName)
    {
        this->printCoeffs(type);
    }
    if(C_S3PQR_.value() != -1){
        Info << "WARNING: B_LES value has been modified." << endl;
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class BasicTurbulenceModel>
bool S3PQR<BasicTurbulenceModel>::read()
{
    if (LESeddyViscosity<BasicTurbulenceModel>::read())
    {
        Ck_.readIfPresent(this->coeffDict());
        Cw_.readIfPresent(this->coeffDict());
        S3PQR_MODEL = word(this->coeffDict().lookup("S3PQR_MODEL"));
        C_S3PQR_.readIfPresent(this->coeffDict());

        return true;
    }

    return false;
}


template<class BasicTurbulenceModel>
tmp<volScalarField> S3PQR<BasicTurbulenceModel>::k() const
{
    return tmp<volScalarField>::New
    (
        IOobject::groupName("k", this->U_.group()),
        (2.0*Ck_/this->Ce_)
       *sqr(this->delta())
       *magSqr(devSymm(fvc::grad(this->U_)))
    );
}


template<class BasicTurbulenceModel>
tmp<volScalarField> S3PQR<BasicTurbulenceModel>::epsilon() const
{
    return tmp<volScalarField>::New
    (
        IOobject::groupName("epsilon", this->U_.group()),
        this->Ce_*k()*sqrt(k())/this->delta()
    );
}


template<class BasicTurbulenceModel>
void S3PQR<BasicTurbulenceModel>::correct()
{
    if (!this->turbulence_)
    {
        return;
    }

    LESeddyViscosity<BasicTurbulenceModel>::correct();
    correctNut();
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace LESModels
} // End namespace Foam

// ************************************************************************* //
