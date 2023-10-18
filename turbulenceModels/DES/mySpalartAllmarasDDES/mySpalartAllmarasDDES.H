/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2022 Upstream CFD GmbH
    Copyright (C) 2019-2022 OpenCFD Ltd.
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
    Foam::LESModels::mySpalartAllmarasDDES

Group
    grpDESTurbulence

Description
    SpalartAllmaras DDES turbulence model for incompressible and compressible
    flows.

    Reference:
    \verbatim
        Spalart, P. R., Deck, S., Shur, M. L., Squires,
        K. D., Strelets, M. K., & Travin, A. (2006).
        A new version of detached-eddy simulation,
        resistant to ambiguous grid densities.
        Theoretical and computational fluid dynamics, 20(3), 181-195.
        DOI:10.1007/s00162-006-0015-0
    \endverbatim

SourceFiles
    mySpalartAllmarasDDES.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_mySpalartAllmarasDDES_H
#define Foam_mySpalartAllmarasDDES_H

#include "SpalartAllmarasBase.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace LESModels
{

/*---------------------------------------------------------------------------*\
                     Class mySpalartAllmarasDDES Declaration
\*---------------------------------------------------------------------------*/

template<class BasicTurbulenceModel>
class mySpalartAllmarasDDES
:
    public SpalartAllmarasDES<BasicTurbulenceModel>
{
    // Private Member Functions

        //- Return the shielding function
        tmp<volScalarField> fd(const volScalarField& magGradU) const;
        tmp<volScalarField> fD(const volScalarField& rd) const;

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

        tmp<volScalarField> SMG_model(const volScalarField& Q) const;
        tmp<volScalarField> WALE_model(const volScalarField& QS, const volScalarField& QG, const volScalarField& V2) const;
        tmp<volScalarField> VREMAN_model(const volScalarField& QS, const volScalarField& QG, const volScalarField& V2) const;
        tmp<volScalarField> SIGMA_model(const volScalarField& QS, const volScalarField& QG, const volScalarField& RG, const volScalarField& V2) const;
        tmp<volScalarField> SIGMA_Alex_model(const volScalarField& QS, const volScalarField& QG, const volScalarField& RG, const volScalarField& V2) const;
        tmp<volScalarField> SIGMA_Fuchs_model(const volScalarField& QS, const volScalarField& QG, const volScalarField& RG, const volScalarField& V2) const;
        tmp<volScalarField> S3PQR_model(const volScalarField& QS, const volScalarField& QG, const volScalarField& RG, const volScalarField& V2, const scalar p, const scalar q, const scalar r, const bool p_sqrt2, const bool q_sqrt2, const bool r_sqrt2 ) const;
        tmp<volScalarField> S3PQ_model(const volScalarField& QS, const volScalarField& QG, const volScalarField& RG, const volScalarField& V2) const;
        tmp<volScalarField> S3PR_model(const volScalarField& QS, const volScalarField& QG, const volScalarField& RG, const volScalarField& V2) const;
        tmp<volScalarField> S3QR_model(const volScalarField& QS, const volScalarField& QG, const volScalarField& RG, const volScalarField& V2) const;

        //- No copy construct
        mySpalartAllmarasDDES(const mySpalartAllmarasDDES&) = delete;

        //- No copy assignment
        void operator=(const mySpalartAllmarasDDES&) = delete;

protected:

    // Protected Data
    
    	
	    
        // Model coefficients

            dimensionedScalar Cd1_;
            dimensionedScalar Cd2_;
            dimensionedScalar C3_;
            dimensionedScalar C4_;
            dimensionedScalar Beta_;
	    dimensionedScalar B_LES;	

        // Model booleans
    	
    	    Switch useDeck_;    
	    word LES_MODEL;

    // Protected Member Functions

        //- Return the production term
        virtual tmp<volScalarField> Stilda
        (
            const volScalarField& chi,
            const volScalarField& fv1,
            const volTensorField& gradU,
            const volScalarField& dTilda
        ) const;

        //- Return the length scale
        virtual tmp<volScalarField> dTilda
        (
            const volScalarField& chi,
            const volScalarField& fv1,
            const volTensorField& gradU
        ) const;
        
        

public:

    typedef typename BasicTurbulenceModel::alphaField alphaField;
    typedef typename BasicTurbulenceModel::rhoField rhoField;
    typedef typename BasicTurbulenceModel::transportModel transportModel;


    //- Runtime type information
    TypeName("mySpalartAllmarasDDES");


    // Constructors

        //- Construct from components
        mySpalartAllmarasDDES
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
    virtual ~mySpalartAllmarasDDES() = default;


    // Member Functions

        //- Read from dictionary
        virtual bool read();

        //- Return the shielding function
        virtual tmp<volScalarField> fd() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace LESModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "mySpalartAllmarasDDES.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
