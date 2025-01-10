/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
     \\/     M anipulation  |
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
    
InClass
    Foam::uniformFrameInvariantAnisotropicStressTensor

\*---------------------------------------------------------------------------*/

#include "uniformFrameInvariantAnisotropicStressTensor.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace anisotropicStressTensorModels
{
    defineTypeNameAndDebug(uniformFrameInvariantAnisotropicStressTensor, 0);

    addToRunTimeSelectionTable
    (
        anisotropicStressTensorModel,
        uniformFrameInvariantAnisotropicStressTensor,
        dictionary
    );
}
}
// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::anisotropicStressTensorModels::
uniformFrameInvariantAnisotropicStressTensor::
uniformFrameInvariantAnisotropicStressTensor
(
    const fvMesh& mesh,
    const volVectorField& U,
    const volScalarField& c,
    const dictionary& dict
)
:
    anisotropicStressTensorModel(mesh, U, c, dict),
    uniformFrameInvariantAnisotropicStressTensorDict_
    (
        dict.subDict(typeName+"Coeffs")
    ),
    lambda1_
    (
        uniformFrameInvariantAnisotropicStressTensorDict_.lookup("lambda1")
    ),
    lambda2_
    (
        uniformFrameInvariantAnisotropicStressTensorDict_.lookup("lambda2")
    ),
    lambda3_
    (
        uniformFrameInvariantAnisotropicStressTensorDict_.lookup("lambda3")
    ),    
    dirU
    (
        IOobject
        (
            "dirU",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE//IOobject::NO_WRITE
        ),
        mesh,
        dimensionedVector
        (
            "dirU", dimensionSet(0,0,0,0,0,0,0), vector(0,0,0)
        ),
        "zeroGradient"
    ),
    dirGradU
    (
        IOobject
        (
            "dirGradU",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE//IOobject::NO_WRITE
        ),
        mesh,
        dimensionedVector
        (
            "dirGradU", dimensionSet(0,0,0,0,0,0,0), vector(0,0,0)
        ),
        "zeroGradient"
    ),
    dirCurlU
    (
        IOobject
        (
            "dirCurlU",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE//IOobject::NO_WRITE
        ),
        mesh,
        dimensionedVector
        (
            "dirCurlU", dimensionSet(0,0,0,0,0,0,0), vector(0,0,0)
        ),
        "zeroGradient"
    )
{}

// * * * * * * * * * * * * * * member functions  * * * * * * * * * * * * * * //

void Foam::anisotropicStressTensorModels::
uniformFrameInvariantAnisotropicStressTensor::updateAnisotropicStressTensor()
{
    // e_1
    dirU = U_/(max(mag(U_), dimensionedScalar("V", dimVelocity, 1e-12)));
    dirU.correctBoundaryConditions();

    // e_3
    dirCurlU =
    (
        fvc::curl(U_)/max
        (
            mag(fvc::curl(U_)),
            dimensionedScalar("curlV", dimVelocity/dimLength, 1e-12)
        )
    );
    dirCurlU.correctBoundaryConditions();

    // e_2    
    dirGradU =
    (
        dirU ^ dirCurlU/max
        (
            mag(dirU ^ dirCurlU),
            dimensionedScalar("curlV", dimless, 1e-12)
        )
    );
    dirGradU.correctBoundaryConditions();
    
    volTensorField Qrot =
    (
        lambda1_*(dirU*dirU)
        + lambda2_*(dirGradU*dirGradU)
        + lambda3_*(dirCurlU*dirCurlU)
    );

    Q_ = Qrot;
    Q_.correctBoundaryConditions();
}

// -------------------------------------------------------------------------//
