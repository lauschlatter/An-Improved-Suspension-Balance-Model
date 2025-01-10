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
    Foam::noNonLocalShearRate

\*---------------------------------------------------------------------------*/

#include "noNonLocalShearRate.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace nonLocalShearRateModels
{
    defineTypeNameAndDebug(noNonLocalShearRate, 0);

    addToRunTimeSelectionTable
    (
        nonLocalShearRateModel,
        noNonLocalShearRate,
        dictionary
    );
}
}
// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::nonLocalShearRateModels::noNonLocalShearRate::noNonLocalShearRate
(
    const fvMesh& mesh,
    const volVectorField& U,
    const volScalarField& c,
    const dictionary& dict
)
:
    nonLocalShearRateModel(mesh, U, c, dict),
    noNonLocalShearRateDict_(dict.subDict(typeName+"Coeffs"))
{}

// * * * * * * * * * * * * * * member functions  * * * * * * * * * * * * * * //

void Foam::nonLocalShearRateModels::
noNonLocalShearRate::updateNonLocalShearRate()
{}

// -------------------------------------------------------------------------//
