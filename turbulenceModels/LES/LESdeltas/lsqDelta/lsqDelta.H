/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2015 OpenFOAM Foundation
    Copyright (C) 2016 OpenCFD Ltd.
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
    Foam::lsqDelta

Description
    Delta formulation based on the velocity gradient that minimizes mesh anisotropies effects.
    lsqDelta is able to detect 2D flow regions; in these regions Delta value is strongly reduced to accelerate the transition from RANS to LES in hybrid RANS/LES simulations.

    Reference:
    \verbatim
	Trias, F. X., Gorobets, A., Silvis, M. H., Verstappen, R. W. C. P., & Oliva, A. (2017). 
	A new subgrid characteristic length for turbulence simulations on anisotropic grids. 
	Physics of Fluids, 29(11).
	DOI:10.1063/1.5012546
    \endverbatim

SourceFiles
    lsqDelta.C

\*---------------------------------------------------------------------------*/

#ifndef lsqDeltaDelta_H
#define lsqDeltaDelta_H

#include "LESdelta.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace LESModels
{

/*---------------------------------------------------------------------------*\
                           Class cubeRootVolDelta Declaration
\*---------------------------------------------------------------------------*/

class lsqDelta
:
    public LESdelta
{
    // Private data

        //- Model coefficient
        volSymmTensorField jacobi;

        label calcInterval_;


    // Private Member Functions

        //- No copy construct
        lsqDelta(const lsqDelta&) = delete;

        //- No copy assignment
        void operator=(const lsqDelta&) = delete;


public:

    //- Runtime type information
    TypeName("lsqDelta");


    // Constructors

        //- Construct from name, turbulenceModel and dictionary
        lsqDelta
        (
            const word& name,
            const turbulenceModel& turbulence,
            const dictionary&
        );


    //- Destructor
    virtual ~lsqDelta() = default;


    // Member Functions

        //- Calculate the delta values
        void calcDelta();
        void calcJacobi();

        //- Read the LESdelta dictionary
        virtual void read(const dictionary&);

        // Correct values
        virtual void correct();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace LESModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
