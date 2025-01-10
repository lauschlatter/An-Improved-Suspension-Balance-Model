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
    Foam::cylindricalCoordinateAnisotropicStressTensor

SourceFiles
    cylindricalCoordinateAnisotropicStressTensor.C
    
Description
    For simulation cases in which e1=etheta, e2=er and e3=ez. Calculates Q 
    using a coordinate system transformation.
    
\*---------------------------------------------------------------------------*/

#ifndef cylindricalCoordinateAnisotropicStressTensor_H
#define cylindricalCoordinateAnisotropicStressTensor_H

#include "volFields.H"
#include "fvCFD.H"
#include "anisotropicStressTensorModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace anisotropicStressTensorModels
{
/*---------------------------------------------------------------------------*\
        Class cylindricalCoordinateAnisotropicStressTensor Declaration
\*---------------------------------------------------------------------------*/

class cylindricalCoordinateAnisotropicStressTensor :
public anisotropicStressTensorModel
{

    protected:

    //- Protected data
    dictionary cylindricalCoordinateAnisotropicStressTensorDict_;
    dimensionedScalar lambda1_;
    dimensionedScalar lambda2_;
    dimensionedScalar lambda3_;
    volScalarField arg;
    volScalarField theta;
    volVectorField er;
    volVectorField etheta;
    volTensorField Qrot;
    

    //- Disallow copy construct
    cylindricalCoordinateAnisotropicStressTensor
    (
        const cylindricalCoordinateAnisotropicStressTensor&
    );


    public:

    //- Runtime type information
    TypeName("cylindricalCoordinateAnisotropicStressTensor");

    
    //- Construct from components
    cylindricalCoordinateAnisotropicStressTensor
    (
        const fvMesh& mesh,
        const volVectorField& U,
        const volScalarField& c,
        const dictionary& dict
    );


    //- Destructor
    virtual ~cylindricalCoordinateAnisotropicStressTensor()
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
