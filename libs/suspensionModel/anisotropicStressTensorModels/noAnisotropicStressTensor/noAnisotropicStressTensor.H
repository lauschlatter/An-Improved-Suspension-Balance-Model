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
    Foam::noAnisotropicStressTensor

SourceFiles
    noAnisotropicStressTensor.C
    
Description
    Foam::noAnisotropicStressTensor

\*---------------------------------------------------------------------------*/

#ifndef noAnisotropicStressTensor_H
#define noAnisotropicStressTensor_H

#include "volFields.H"
#include "anisotropicStressTensorModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace anisotropicStressTensorModels
{
/*---------------------------------------------------------------------------*\
                    Class noAnisotropicStressTensor Declaration
\*---------------------------------------------------------------------------*/

class noAnisotropicStressTensor : public anisotropicStressTensorModel
{

    protected:

    //- Protected data
    dictionary noAnisotropicStressTensorDict_;
 

    //- Disallow copy construct
    noAnisotropicStressTensor(const noAnisotropicStressTensor&);


    public:

    //- Runtime type information
    TypeName("noAnisotropicStressTensor");


    //- Construct from components
    noAnisotropicStressTensor
    (
        const fvMesh& mesh,
        const volVectorField& U,
        const volScalarField& c,
        const dictionary& dict
    );


    //- Destructor
    virtual ~noAnisotropicStressTensor()
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
