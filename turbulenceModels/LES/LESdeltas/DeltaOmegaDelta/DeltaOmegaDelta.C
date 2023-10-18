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

\*---------------------------------------------------------------------------*/

#include "DeltaOmegaDelta.H"
#include "fvcCurl.H"
#include "maxDeltaxyz.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace LESModels
{
    defineTypeNameAndDebug(DeltaOmegaDelta, 0);
    addToRunTimeSelectionTable(LESdelta, DeltaOmegaDelta, dictionary);
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::LESModels::DeltaOmegaDelta::calcDelta()
{
    const fvMesh& mesh = turbulenceModel_.mesh();

    label nD = mesh.nGeometricD();

    const volVectorField& U = turbulenceModel_.U();
    const volVectorField curlU = fvc::curl(U);
    const volVectorField n_w(curlU/mag(curlU));

    const cellList& cells = mesh.cells();
    const vectorField& cellC = mesh.cellCentres();
    const pointField& pts = mesh.points();

    scalarField hmax(cells.size());

    forAll(cells, celli)
    {
        scalar deltaMaxTmp = 0.0;
                const labelList& cPoints = mesh.cellPoints()[celli];
        const point& cc = cellC[celli];

                List<vector> I(cPoints.size());
       forAll(cPoints, i)
        {
            label pointi = cPoints[i];
            const point& v = pts[pointi];

            I[i] = ( n_w[celli] ^ (v - cc));
        }

                forAll(I, n)
                {
                        forAll(I, m) {
                                if (m>n){
                                        scalar tmp = mag(I[n]-I[m]);
                        if (tmp > deltaMaxTmp)
                        {
                            deltaMaxTmp = tmp;
                        }
                                }
                        }
                }

        hmax[celli] = deltaCoeff_*Foam::sqrt(1.0/3.0)*deltaMaxTmp;
    }

    if (nD == 3)
    {
        delta_.primitiveFieldRef() = hmax;
    }
    else if (nD == 2)
    {
        WarningInFunction
            << "Case is 2D, LES is not strictly applicable" << nl
            << endl;

        delta_.primitiveFieldRef() = hmax;
    }
    else
    {
        FatalErrorInFunction
            << "Case is not 3D or 2D, LES is not applicable"
            << exit(FatalError);
    }

    // Handle coupled boundaries
    delta_.correctBoundaryConditions();


}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::LESModels::DeltaOmegaDelta::DeltaOmegaDelta
(
    const word& name,
    const turbulenceModel& turbulence,
    const dictionary& dict
)
:
    LESdelta(name, turbulence),
    hmaxPtr_(nullptr),
    deltaCoeff_
    (
        dict.optionalSubDict(type() + "Coeffs").getOrDefault<scalar>
        (
            "deltaCoeff",
            1
        )
    ),
    requireUpdate_
    (
        dict.optionalSubDict(type() + "Coeffs").getOrDefault<bool>
        (
            "requireUpdate",
            true
        )
    ),
    calcInterval_
    (
        dict.optionalSubDict(type() + "Coeffs").lookupOrDefault<label>
        (
            "calcInterval",
            1
        )
    )
{
    if (dict.optionalSubDict(type() + "Coeffs").found("hmax"))
    {
        // User-defined hmax
        hmaxPtr_ =
            LESdelta::New
            (
                IOobject::groupName("hmax", turbulence.U().group()),
                turbulence,
                dict.optionalSubDict("hmaxCoeffs"),
                "hmax"
            );
    }
    else
    {
        Info<< "Employing " << maxDeltaxyz::typeName << " for hmax" << endl;

        hmaxPtr_.reset
        (
            new maxDeltaxyz
            (
                IOobject::groupName("hmax", turbulence.U().group()),
                turbulence,
                dict.optionalSubDict("hmaxCoeffs")
            )
        );
    }

    calcDelta();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::LESModels::DeltaOmegaDelta::read(const dictionary& dict)
{
    const dictionary& coeffsDict = dict.optionalSubDict(type() + "Coeffs");

    coeffsDict.readIfPresent<scalar>("deltaCoeff", deltaCoeff_);
    coeffsDict.readIfPresent<bool>("requireUpdate", requireUpdate_);
    coeffsDict.readIfPresent<label>("calcInterval", calcInterval_);

    calcDelta();
}


void Foam::LESModels::DeltaOmegaDelta::correct()
{
    if (turbulenceModel_.mesh().changing() && requireUpdate_)
    {
        hmaxPtr_->correct();
    }

    scalar writeLESdelta_ = turbulenceModel_.mesh().time().controlDict().lookupOrDefault("writeLESdelta", false);
    if(writeLESdelta_ && turbulenceModel_.mesh().time().outputTime() )
    {                                                             
        delta_.write();
    }

    if (turbulenceModel_.mesh().time().timeIndex() % calcInterval_ == 0)
    {
        calcDelta();
    }

}


// ************************************************************************* //
