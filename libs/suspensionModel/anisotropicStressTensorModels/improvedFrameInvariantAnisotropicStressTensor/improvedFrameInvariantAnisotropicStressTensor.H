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
    Foam::improvedFrameInvariantAnisotropicStressTensor

SourceFiles
    improvedFrameInvariantAnisotropicStressTensor.C
    
Description
    Uses the frame-invariant formulation of the anisotropic stress tensor with a
    correction for regions with near-constant velocity.

\*---------------------------------------------------------------------------*/

#ifndef improvedFrameInvariantAnisotropicStressTensor_H
#define improvedFrameInvariantAnisotropicStressTensorr_H

#include "volFields.H"
#include "fvCFD.H"
#include "anisotropicStressTensorModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace anisotropicStressTensorModels
{
/*---------------------------------------------------------------------------*\
        Class improvedFrameInvariantAnisotropicStressTensor Declaration
\*---------------------------------------------------------------------------*/

class improvedFrameInvariantAnisotropicStressTensor :
public anisotropicStressTensorModel
{

    protected:

    //- Protected data
    dictionary improvedFrameInvariantAnisotropicStressTensorDict_;
    dimensionedScalar lambda1_;
    dimensionedScalar lambda2_;
    dimensionedScalar lambda3_;
    dimensionedScalar cut_;
    dimensionedScalar decay_;
    dimensionedScalar rTol_;
    dimensionedScalar gUCrit_;
    volScalarField epsilon_;
    volScalarField fx_;
    volScalarField y_;
    volVectorField dirU_;
    volVectorField dirGradU_;
    volVectorField dirCurlU_;
    volVectorField curlU_;
    volTensorField gradU_;
    volScalarField magGradU_;
    volScalarField magCurlU_;
    volScalarField magU_;
    volScalarField l1_;
    volScalarField l2_;
    volScalarField l3_;
    bool is2D_;
    vector e32D_;
    Switch useGradients_;


    //- Disallow copy construct
    improvedFrameInvariantAnisotropicStressTensor
    (
        const improvedFrameInvariantAnisotropicStressTensor&
    );


    public:

    //- Runtime type information
    TypeName("improvedFrameInvariantAnisotropicStressTensor");


    //- Construct from components
    improvedFrameInvariantAnisotropicStressTensor
    (
        const fvMesh& mesh,
        const volVectorField& U,
        const volScalarField& c,
        const dictionary& dict
    );


    //- Destructor
    virtual ~improvedFrameInvariantAnisotropicStressTensor()
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
