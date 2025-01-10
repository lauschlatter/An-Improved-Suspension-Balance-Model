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

Class
    Foam::uniformFrameInvariantAnisotropicStressTensor

SourceFiles
    uniformFrameInvariantAnisotropicStressTensor.C
    
Description
    Uses the frame-invariant formulation of the anisotropic stress tensor,
    in which Q = lambda1*e1*e1 + lambda2*e2*e2 + lambda3*e3*e3
    with e1, e2 and e3 determined based on the calculated velocity field.
    DOES NOT include the correction for regions with near-constant velocity.

\*---------------------------------------------------------------------------*/

#ifndef uniformFrameInvariantAnisotropicStressTensor_H
#define uniformFrameInvariantAnisotropicStressTensor_H

#include "volFields.H"
#include "fvCFD.H"
#include "anisotropicStressTensorModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace anisotropicStressTensorModels
{
/*---------------------------------------------------------------------------*\
        Class uniformFrameInvariantAnisotropicStressTensor Declaration
\*---------------------------------------------------------------------------*/

class uniformFrameInvariantAnisotropicStressTensor :
public anisotropicStressTensorModel
{

    protected:

    //- Protected data
    dictionary uniformFrameInvariantAnisotropicStressTensorDict_;
    dimensionedScalar lambda1_;
    dimensionedScalar lambda2_;
    dimensionedScalar lambda3_;
    volVectorField dirU;
    volVectorField dirGradU;
    volVectorField dirCurlU;


    //- Disallow copy construct
    uniformFrameInvariantAnisotropicStressTensor
    (
        const uniformFrameInvariantAnisotropicStressTensor&
    );


    public:

    //- Runtime type information
    TypeName("uniformFrameInvariantAnisotropicStressTensor");


    //- Construct from components
    uniformFrameInvariantAnisotropicStressTensor
    (
        const fvMesh& mesh,
        const volVectorField& U,
        const volScalarField& c,
        const dictionary& dict
    );


    //- Destructor
    virtual ~uniformFrameInvariantAnisotropicStressTensor()
    {}


    //- Member Functions
    virtual void updateAnisotropicStressTensor();

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
} // End namespace anisotropicStressTensorModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
