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
    Foam::constantNonLocalShearRate

\*---------------------------------------------------------------------------*/

#include "constantNonLocalShearRate.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace nonLocalShearRateModels
{
    defineTypeNameAndDebug(constantNonLocalShearRate, 0);

    addToRunTimeSelectionTable
    (
        nonLocalShearRateModel,
        constantNonLocalShearRate,
        dictionary
    );
}
}
// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::nonLocalShearRateModels::
constantNonLocalShearRate::constantNonLocalShearRate
(
    const fvMesh& mesh,
    const volVectorField& U,
    const volScalarField& c,
    const dictionary& dict
)
:
    nonLocalShearRateModel(mesh, U, c, dict),
    constantNonLocalShearRateDict_(dict.subDict(typeName+"Coeffs")),
    epsilon_(constantNonLocalShearRateDict_.lookup("epsilon"))
{
    dimensionedScalar Umx =
        max(max(max(U_.component(0), U_.component(1)), U_.component(2)));

    gammaNL_ = epsilon_*Umx;
}

// * * * * * * * * * * * * * * member functions  * * * * * * * * * * * * * * //

void Foam::nonLocalShearRateModels::
constantNonLocalShearRate::updateNonLocalShearRate()
{
    dimensionedScalar Umx =
        max(max(max(U_.component(0), U_.component(1)), U_.component(2)));

    gammaNL_ = epsilon_*Umx;
}

// -------------------------------------------------------------------------//
