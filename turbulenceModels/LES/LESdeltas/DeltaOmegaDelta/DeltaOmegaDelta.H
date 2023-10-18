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
    Foam::LESModels::DeltaOmegaDelta

Description
    Delta calculated by taking the maximum distance between the cell centre
    and any face centre.  For a regular hex cell, the computed delta will
    equate to half of the cell width; accordingly, the deltaCoeff model
    coefficient should be set to 2 for this case.

SourceFiles
    DeltaOmegaDelta.C

\*---------------------------------------------------------------------------*/

#ifndef LESModels_DeltaOmegaDelta_H
#define LESModels_DeltaOmegaDelta_H

#include "LESdelta.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace LESModels
{

/*---------------------------------------------------------------------------*\
                    Class DeltaOmegaDelta Declaration
\*---------------------------------------------------------------------------*/

class DeltaOmegaDelta
:
    public LESdelta
{
    // Private Data

        //- Run-time selectable delta for hmax
        //  Defaults to the maxDeltaxyz model if not supplied
        autoPtr<LESdelta> hmaxPtr_;

        //- Model coefficient
        scalar deltaCoeff_;
        label calcInterval_;

        //- Flag to indicate whether hmax requires updating
        bool requireUpdate_;


    // Private Member Functions

        //- Calculate the delta values
        void calcDelta();

        //- No copy construct
        DeltaOmegaDelta(const DeltaOmegaDelta&) = delete;

        //- No copy assignment
        void operator=(const DeltaOmegaDelta&) = delete;


public:

    //- Runtime type information
    TypeName("DeltaOmega");


    // Constructors

        //- Construct from name, turbulenceModel and dictionary
        DeltaOmegaDelta
        (
            const word& name,
            const turbulenceModel& turbulence,
            const dictionary&
        );


    //- Destructor
    virtual ~DeltaOmegaDelta() = default;


    // Member Functions

        //- Read the LESdelta dictionary
        virtual void read(const dictionary&);

        //- Return the hmax delta field
        const volScalarField& hmax() const
        {
            return hmaxPtr_();
        }

        // Correct values
        void correct();

	inline const label& getCalcInterval() const
        {
                        return calcInterval_;
        }

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace LESModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
