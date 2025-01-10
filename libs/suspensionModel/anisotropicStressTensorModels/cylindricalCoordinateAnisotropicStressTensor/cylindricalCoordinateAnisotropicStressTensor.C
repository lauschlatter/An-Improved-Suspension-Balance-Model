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
    Foam::cylindricalCoordinateAnisotropicStressTensor

\*---------------------------------------------------------------------------*/

#include "cylindricalCoordinateAnisotropicStressTensor.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace anisotropicStressTensorModels
{
    defineTypeNameAndDebug(cylindricalCoordinateAnisotropicStressTensor, 0);

    addToRunTimeSelectionTable
    (
        anisotropicStressTensorModel,
        cylindricalCoordinateAnisotropicStressTensor,
        dictionary
    );
}
}
// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::anisotropicStressTensorModels::
cylindricalCoordinateAnisotropicStressTensor::
cylindricalCoordinateAnisotropicStressTensor
(
    const fvMesh& mesh,
    const volVectorField& U,
    const volScalarField& c,
    const dictionary& dict
)
:
    anisotropicStressTensorModel(mesh, U, c, dict),
    cylindricalCoordinateAnisotropicStressTensorDict_
    (
        dict.subDict(typeName+"Coeffs")
    ),
    lambda1_
    (
        cylindricalCoordinateAnisotropicStressTensorDict_.lookup("lambda1")
    ),
    lambda2_
    (
        cylindricalCoordinateAnisotropicStressTensorDict_.lookup("lambda2")
    ),
    lambda3_
    (
        cylindricalCoordinateAnisotropicStressTensorDict_.lookup("lambda3")
    ),
    arg
    (
        IOobject
        (
            "arg",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        0.0,
        "calculated"
    ),
    theta
    (
        IOobject
        (
            "theta",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        0.0,
        "calculated"
    ),
    er
    (
        IOobject
        (
            "er",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        vector(0,0,0),
        "calculated"
    ),
    etheta
    (
        IOobject
        (
            "etheta",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        vector(0,0,0),
        "calculated"
    ),
    Qrot
    (
        IOobject
        (
            "Qrot",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        tensor(0,0,0,0,0,0,0,0,0),
        "calculated"
    )
{
    // Cartesian unit vectors
    vector i (1,0,0);
    vector j (0,1,0);
    vector k (0,0,1);
    
    // Calculating er and etheta 
    arg = mesh_.C().component(vector::Y)/mesh_.C().component(vector::X);
    theta = atan(arg);
    er = cos(theta)*i + sin(theta)*j;
    etheta = cos(theta)*j - sin(theta)*i;
    
    // Using er and etheta to calculate Q
    Qrot = lambda1_*(k*k) + lambda2_*(er*er) + lambda3_*(etheta*etheta);
    Q_ = Qrot;
    Q_.correctBoundaryConditions();
}

// * * * * * * * * * * * * * * member functions  * * * * * * * * * * * * * * //

void Foam::anisotropicStressTensorModels::
cylindricalCoordinateAnisotropicStressTensor::updateAnisotropicStressTensor()
{}

// -------------------------------------------------------------------------//
