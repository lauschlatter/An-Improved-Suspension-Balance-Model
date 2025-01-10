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
    Foam::morrisBoulay

\*---------------------------------------------------------------------------*/

#include "morrisBoulay.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace shearViscosityModels
{
    defineTypeNameAndDebug(morrisBoulay, 0);

    addToRunTimeSelectionTable
    (
        shearViscosityModel,
        morrisBoulay,
        dictionary
    );
}
}
// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::shearViscosityModels::morrisBoulay::morrisBoulay
(
    const fvMesh& mesh,
    const volVectorField& U,
    const volScalarField& c,
    const dictionary& dict
)
:
    shearViscosityModel(mesh, U, c, dict),
    morrisBoulayDict_(dict.subDict(typeName+"Coeffs")),
    cMax_(dict.lookup("cMax")),
    nu0_(dict.lookup("nu0")),
    aP_(morrisBoulayDict_.lookup("aP")),
    bP_(morrisBoulayDict_.lookup("bP")),
    cP_(morrisBoulayDict_.lookup("cP"))
{
    cMax_ = morrisBoulayDict_.lookupOrDefault("cMax", cMax_);
    nu0_ = morrisBoulayDict_.lookupOrDefault("nu0", nu0_);
}

// * * * * * * * * * * * * * * member functions  * * * * * * * * * * * * * * //

void Foam::shearViscosityModels::morrisBoulay::updateShearViscosity()
{
    nus_ =
    (
        nu0_*(aP_ + bP_*cMax_*pow(max((cMax_/c_ - 1), 1e-04), -1.0)
        + cP_*sqr(c_/cMax_)*pow(max((1.0 - (c_/cMax_)), 1e-04), - 2.0) + 1.)
    );

    nup_ = nus_ - nu0_;
}

// -------------------------------------------------------------------------//
