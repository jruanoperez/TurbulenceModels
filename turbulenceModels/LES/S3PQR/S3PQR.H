/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2022 Upstream CFD GmbH
    Copyright (C) 2022 OpenCFD Ltd.
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

Class
    Foam::LESModels::S3PQR

Group
    grpLESTurbulence

Description
    The S3PQR SGS model.

    Reference:
    \verbatim
	Trias, F. X., Folch, D., Gorobets, A., & Oliva, A. (2015). 
	Building proper invariants for eddy-viscosity subgrid-scale models. 
	Physics of Fluids, 27(6).
        DOI:10.1063/1.4921817
    \endverbatim

Note
    The default value of the C_S3PQR constant implemented was calibrated for
    OpenFOAM using decaying isotropic turbulence and is slightly higher than
    the value suggested in the reference publication.

SourceFiles
    S3PQR.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_LESModels_S3PQR_H
#define Foam_LESModels_S3PQR_H

#include "LESModel.H"
#include "LESeddyViscosity.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace LESModels
{

/*---------------------------------------------------------------------------*\
                           Class S3PQR Declaration
\*---------------------------------------------------------------------------*/

template<class BasicTurbulenceModel>
class S3PQR
:
    public LESeddyViscosity<BasicTurbulenceModel>
{
    // Private Member Functions

	tmp<volScalarField> calc_QG(const volTensorField& gradU) const;
        tmp<volScalarField> calc_QS(const volTensorField& gradU) const;
        tmp<volScalarField> calc_RG(const volTensorField& gradU) const;
        tmp<volScalarField> calc_RS(const volTensorField& gradU) const;
        tmp<volScalarField> calc_V2(const volTensorField& gradU) const;
        tmp<volScalarField> calc_first_invariant(const volTensorField& gradU) const;
        tmp<volScalarField> calc_second_invariant(const volTensorField& gradU) const;
        tmp<volScalarField> calc_second_invariant(const volSymmTensorField& gradU) const;
        tmp<volScalarField> calc_third_invariant(const volTensorField& gradU) const;
        tmp<volScalarField> calc_third_invariant(const volSymmTensorField& gradU) const;
        vector solve_cubic_eq(scalar a, scalar b, scalar c, scalar d ) const;

	tmp<volScalarField> S3PQR_model(const volScalarField& QS, const volScalarField& QG, const volScalarField& RG, const volScalarField& V2, const scalar p, const scalar q, const scalar r, const bool p_sqrt2, const bool q_sqrt2, const bool r_sqrt2 ) const;
        tmp<volScalarField> S3PQ_model(const volScalarField& QS, const volScalarField& QG, const volScalarField& RG, const volScalarField& V2) const;
        tmp<volScalarField> S3PR_model(const volScalarField& QS, const volScalarField& QG, const volScalarField& RG, const volScalarField& V2) const;
        tmp<volScalarField> S3QR_model(const volScalarField& QS, const volScalarField& QG, const volScalarField& RG, const volScalarField& V2) const;

        //- No copy construct
        S3PQR(const S3PQR&) = delete;

        //- No copy assignment
        void operator=(const S3PQR&) = delete;


protected:

    // Protected Data

	dimensionedScalar Cw_;
	dimensionedScalar Ck_;
        dimensionedScalar C_S3PQR_;
	word S3PQR_MODEL;

    // Protected Member Functions

        //- Update the SGS eddy-viscosity
        virtual void correctNut();


public:

    typedef typename BasicTurbulenceModel::alphaField alphaField;
    typedef typename BasicTurbulenceModel::rhoField rhoField;
    typedef typename BasicTurbulenceModel::transportModel transportModel;


    //- Runtime type information
    TypeName("S3PQR");


    // Constructors

        //- Construct from components
        S3PQR
        (
            const alphaField& alpha,
            const rhoField& rho,
            const volVectorField& U,
            const surfaceScalarField& alphaRhoPhi,
            const surfaceScalarField& phi,
            const transportModel& transport,
            const word& propertiesName = turbulenceModel::propertiesName,
            const word& type = typeName
        );


    //- Destructor
    virtual ~S3PQR()
    {}


    // Member Functions

        //- Read model coefficients if they have changed
        virtual bool read();

        //- Return SGS kinetic energy
        virtual tmp<volScalarField> k() const;

        //- Return SGS disipation rate
        virtual tmp<volScalarField> epsilon() const;

        //- Correct Eddy-Viscosity and related properties
        virtual void correct();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace LESModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "S3PQR.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
