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
    Foam::frameDependentAnisotropicStressTensor

\*---------------------------------------------------------------------------*/

#include "frameDependentAnisotropicStressTensor.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace anisotropicStressTensorModels
{
    defineTypeNameAndDebug(frameDependentAnisotropicStressTensor, 0);

    addToRunTimeSelectionTable
    (
        anisotropicStressTensorModel,
        frameDependentAnisotropicStressTensor,
        dictionary
    );
}
}
// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::anisotropicStressTensorModels::
frameDependentAnisotropicStressTensor::frameDependentAnisotropicStressTensor
(
    const fvMesh& mesh,
    const volVectorField& U,
    const volScalarField& c,
    const dictionary& dict
)
:
    anisotropicStressTensorModel(mesh, U, c, dict),
    frameDependentAnisotropicStressTensorDict_(dict.subDict(typeName+"Coeffs")),
    lambda1_
    (
        frameDependentAnisotropicStressTensorDict_.lookup("lambda1")
    ),
    lambda2_
    (
        frameDependentAnisotropicStressTensorDict_.lookup("lambda2")
    ),
    lambda3_
    (
        frameDependentAnisotropicStressTensorDict_.lookup("lambda3")
    )
{
    vector i (1,0,0);
    vector j (0,1,0);
    vector k (0,0,1);

    Q_ = lambda1_*i*i + lambda2_*j*j + lambda3_*k*k;
    Q_.correctBoundaryConditions();
}

// * * * * * * * * * * * * * * member functions  * * * * * * * * * * * * * * //

void Foam::anisotropicStressTensorModels::
frameDependentAnisotropicStressTensor::updateAnisotropicStressTensor()
{}

// -------------------------------------------------------------------------//
