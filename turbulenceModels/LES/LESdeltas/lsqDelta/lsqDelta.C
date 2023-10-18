/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2016-2020 OpenCFD Ltd.
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

#include "lsqDelta.H"
#include "addToRunTimeSelectionTable.H"
#include "fvc.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace LESModels
{
    defineTypeNameAndDebug(lsqDelta, 0);
    addToRunTimeSelectionTable(LESdelta, lsqDelta, dictionary);
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::LESModels::lsqDelta::calcJacobi()
{
    const fvMesh& mesh = turbulenceModel_.mesh();

    const cellList& cells = mesh.cells();
    //const vectorField& cellC = mesh.cellCentres();
        const labelUList& owner = mesh.owner();                    
        const labelUList& neighbour = mesh.neighbour();            
        const vectorField& Sf = mesh.Sf();

        surfaceScalarField w = mesh.surfaceInterpolation::weights();

        vectorField gradUnstrCoeff = vectorField(mesh.nCells());

        scalar count = 0;
        forAll ( cells, celli)
        {
                count += 1 ;
                //Info << count << ") New Cell" << cellC[celli] << endl;

        const labelList& cFaces = cells[celli];
                gradUnstrCoeff[celli] = vector(0.0,0.0,0.0);

        forAll(cFaces, cFacei)
        {
            label facei = cFaces[cFacei];
                        label facePatch = mesh.boundaryMesh().whichPatch(facei);
                        vector value(0.0,0.0,0.0);

                        if(facePatch==-1){
                                if(celli==owner[facei])
                                        value = (scalar(1.0)-w[facei])*Sf[facei];
                                else if(celli==neighbour[facei])
                                        value = w[facei]*Sf[facei];
                                //Info << "Face[" << facePatch << "]: " << facei << ", " << w[facei] <<", " << Sf[facei] << ", " << value << endl;
                        }

                        forAll(value,i){
                                if(value[i]<0.0)
                                        value[i]*=-1.0;
                        }
                        //Info << "Face, " << w[facei] <<", " << Sf[facei] << ", " << value << ", " << 1.0/(1e-13+deltaCoeffs[facei])<< endl;			
                        gradUnstrCoeff[celli] += value;
        }
        }

        forAll(mesh.boundary(), patchi)
    {
        const labelUList& pFaceCells = mesh.boundary()[patchi].faceCells();
        const vectorField& pSf = mesh.Sf().boundaryField()[patchi];
                const scalarField& pwb = w.boundaryField()[patchi];

                vector value(0.0,0.0,0.0);
        forAll(mesh.boundary()[patchi], facei)
        {
                        value = pSf[facei]*pwb[facei];
                        forAll(value,i){
                                if(value[i]<0.0)
                                        value[i]*=-1.0;
                        }

            gradUnstrCoeff[pFaceCells[facei]] += value;
                        //Info << "Boundary[" << patchi << "], " << facei << ": " << pwb[facei] <<", " << pSf[facei] << ", " << value << endl;
        }
        }
        gradUnstrCoeff /= mesh.V();

        forAll ( cells, celli)
        {
                scalar xDim = 1.0/gradUnstrCoeff[celli][0];
                scalar yDim = 1.0/gradUnstrCoeff[celli][1];
                scalar zDim = 1.0/gradUnstrCoeff[celli][2];
                //Info << "JACOBI: " << cellC[celli][0] << "\t" << cellC[celli][1] << "\t" << cellC[celli][2] << "\t"  << xDim <<"\t" << yDim << "\t" << zDim << endl;

                jacobi[celli] = symmTensor(xDim,0,0,yDim,0,zDim);
        }
}


void Foam::LESModels::lsqDelta::calcDelta()
{

const fvMesh& mesh = turbulenceModel_.mesh();

    label nD = mesh.nGeometricD();

    const volVectorField& U = turbulenceModel_.U();

        tmp<volTensorField> tgradU = fvc::grad(U);
        const volTensorField& gradUT = tgradU();

        const volTensorField jgradUTgradU((jacobi&gradUT) & gradUT.T() );
        const volTensorField gradUTgradU(gradUT & gradUT.T() );

        //We should consider the possibility of having (gradUTgradU && gradUTgradU) == 0. It have not been implemented yet.
    if (nD == 3)
    {
                delta_.primitiveFieldRef() = Foam::sqrt((jgradUTgradU && jgradUTgradU) / (gradUTgradU && gradUTgradU));
    }
    else if (nD == 2)
    {
        WarningInFunction
            << "Case is 2D, LES is not strictly applicable" << nl
            << endl;

                delta_.primitiveFieldRef() = Foam::sqrt((jgradUTgradU && jgradUTgradU) / (gradUTgradU && gradUTgradU));
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

Foam::LESModels::lsqDelta::lsqDelta
(
    const word& name,
    const turbulenceModel& turbulence,
    const dictionary& dict
)
:
    LESdelta(name, turbulence),
        jacobi
    (
        IOobject
        (
            "jacobi",
            turbulence.mesh().time().timeName(),
            turbulence.mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        turbulence.mesh(),
        dimensionedSymmTensor("jacobi", dimLength, Zero)
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
    calcJacobi();
    calcDelta();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::LESModels::lsqDelta::read(const dictionary& dict)
{
    dict.optionalSubDict(type() + "Coeffs").readIfPresent<label>
    (
        "calcInterval",
        calcInterval_
    );

    calcJacobi();
    calcDelta();
}


void Foam::LESModels::lsqDelta::correct()
{
    if (turbulenceModel_.mesh().changing())
    {
        calcJacobi();
    }
    if (turbulenceModel_.mesh().time().timeIndex() % calcInterval_ == 0)
    {
        calcDelta();
    }
}


// ************************************************************************* //
