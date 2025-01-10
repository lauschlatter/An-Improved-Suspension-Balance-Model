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
    Foam::frameDependentAnisotropicStressTensor

SourceFiles
    frameDependentAnisotropicStressTensor.C
    
Description
    Uses the frame-dependent formulation of Morris and Boulay (1999)
    with Q = [(lambda1 0 0) (0 lambda2 0) (0 0 lambda3)]

\*---------------------------------------------------------------------------*/

#ifndef frameDependentAnisotropicStressTensor_H
#define frameDependentAnisotropicStressTensor_H

#include "volFields.H"
#include "fvCFD.H"
#include "anisotropicStressTensorModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace anisotropicStressTensorModels
{
/*---------------------------------------------------------------------------*\
            Class frameDependentAnisotropicStressTensor Declaration
\*---------------------------------------------------------------------------*/

class frameDependentAnisotropicStressTensor :
public anisotropicStressTensorModel
{

    protected:

    //- Protected data
    dictionary frameDependentAnisotropicStressTensorDict_;
    dimensionedScalar lambda1_;
    dimensionedScalar lambda2_;
    dimensionedScalar lambda3_;


    //- Disallow copy construct
    frameDependentAnisotropicStressTensor
    (
        const frameDependentAnisotropicStressTensor&
    );


    public:

    //- Runtime type information
    TypeName("frameDependentAnisotropicStressTensor");


    //- Construct from components
    frameDependentAnisotropicStressTensor
    (
        const fvMesh& mesh,
        const volVectorField& U,
        const volScalarField& c,
        const dictionary& dict
    );


    //- Destructor
    virtual ~frameDependentAnisotropicStressTensor()
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
