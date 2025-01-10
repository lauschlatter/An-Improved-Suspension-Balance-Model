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
    Foam::morrisBoulayNormal

\*---------------------------------------------------------------------------*/

#include "morrisBoulayNormal.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace normalViscosityModels
{
    defineTypeNameAndDebug(morrisBoulayNormal, 0);

    addToRunTimeSelectionTable
    (
        normalViscosityModel,
        morrisBoulayNormal,
        dictionary
    );
}
}
// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::normalViscosityModels::morrisBoulayNormal::morrisBoulayNormal
(
    const fvMesh& mesh,
    const volVectorField& U,
    const volScalarField& c,
    const dictionary& dict
)
:
    normalViscosityModel(mesh, U, c, dict),
    morrisBoulayNormalDict_(dict.subDict(typeName+"Coeffs")),
    cMax_(dict.lookup("cMax")),
    nu0_(dict.lookup("nu0"))
{
    cMax_ = morrisBoulayNormalDict_.lookupOrDefault("cMax", cMax_);
    nu0_ = morrisBoulayNormalDict_.lookupOrDefault("nu0", nu0_);
}

// * * * * * * * * * * * * * * member functions  * * * * * * * * * * * * * * //

void Foam::normalViscosityModels::morrisBoulayNormal::updateNormalViscosity()
{
    nun_ = 0.75*nu0_*sqr(c_/cMax_)*pow(max((1.0 - (c_/cMax_)), 1e-04), -2.0);
    dNunDc_ = 1.5*nu0_*(c_/sqr(cMax_))*pow(max((1.0 - c_/cMax_), 1e-04), -3.0);
}

// -------------------------------------------------------------------------//
