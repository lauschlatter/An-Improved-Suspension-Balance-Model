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
    Foam::normalViscosityModel

\*---------------------------------------------------------------------------*/

#include "normalViscosityModel.H"
#include "volFields.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(normalViscosityModel, 0);
    defineRunTimeSelectionTable(normalViscosityModel, dictionary);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::normalViscosityModel::normalViscosityModel
(
    const fvMesh& mesh,
    const volVectorField& U,
    const volScalarField& c,
    const dictionary& dict
)
:
    mesh_(mesh),
    U_(U),
    c_(c),
    nun_
    (
        IOobject
        (
            "nun",
            mesh.time().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("nun", dimensionSet(0,2,-1,0,0), 0.0),
        "calculated"
    ),
    dNunDc_
    (
        IOobject
        (
            "dNunDc",
            mesh.time().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("dNunDc", dimensionSet(0,2,-1,0,0), 0.0),
        "calculated"
    )
{}

// -------------------------------------------------------------------------//

// ************************************************************************* //
