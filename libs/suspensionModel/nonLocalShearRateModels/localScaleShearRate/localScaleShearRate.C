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
    Foam::localScaleShearRate

\*---------------------------------------------------------------------------*/

#include "localScaleShearRate.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace nonLocalShearRateModels
{
    defineTypeNameAndDebug(localScaleShearRate, 0);

    addToRunTimeSelectionTable
    (
        nonLocalShearRateModel,
        localScaleShearRate,
        dictionary
    );
}
}
// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::nonLocalShearRateModels::localScaleShearRate::localScaleShearRate
(
    const fvMesh& mesh,
    const volVectorField& U,
    const volScalarField& c,
    const dictionary& dict
)
:
    nonLocalShearRateModel(mesh, U, c, dict),
    localScaleShearRateDict_(dict.subDict(typeName+"Coeffs")),
    a_(dict.lookup("a")),
    k_(dict.subDict(typeName+"Coeffs").lookup("k"))
{
    volScalarField magU = mag(U_);
    gammaNL_ = k_*magU/(2*a_);
}

// * * * * * * * * * * * * * * member functions  * * * * * * * * * * * * * * //

void Foam::nonLocalShearRateModels::
localScaleShearRate::updateNonLocalShearRate()
{
    volScalarField magU = mag(U_);
    gammaNL_ = k_*magU/(2*a_);
}

// -------------------------------------------------------------------------//
