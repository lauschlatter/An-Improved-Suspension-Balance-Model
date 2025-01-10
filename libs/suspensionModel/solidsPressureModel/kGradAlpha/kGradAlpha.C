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
    Foam::kGradAlpha

\*---------------------------------------------------------------------------*/

#include "kGradAlpha.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace solidsPressureModels
{
    defineTypeNameAndDebug(kGradAlpha, 0);

    addToRunTimeSelectionTable
    (
        solidsPressureModel,
        kGradAlpha,
        dictionary
    );
}
}
// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::solidsPressureModels::kGradAlpha::kGradAlpha
(
    const fvMesh& mesh,
    const volVectorField& U,
    const volScalarField& c,
    const dictionary& dict
)
:
    solidsPressureModel(mesh, U, c, dict),
    kGradAlphaDict_(dict.subDict(typeName+"Coeffs")),
    k_(dict.subDict(typeName+"Coeffs").lookup("k"))
{}

// * * * * * * * * * * * * * * member functions  * * * * * * * * * * * * * * //

void Foam::solidsPressureModels::kGradAlpha::updateSolidsPressure()
{
    Fp_ = -k_*fvc::grad(c_);
    Fp_.correctBoundaryConditions();

    Fpf_ = -k_*(fvc::snGrad(c_)*mesh_.magSf());
}

// -------------------------------------------------------------------------//
