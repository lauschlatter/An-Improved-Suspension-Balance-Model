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
    Foam::morrisBoulayNormal

SourceFiles
    morrisBoulayNormal.C
    
Description
    Foam::morrisBoulayNormal

\*---------------------------------------------------------------------------*/

#ifndef morrisBoulayNormal_H
#define morrisBoulayNormal_H

#include "volFields.H"
#include "normalViscosityModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace normalViscosityModels
{
/*---------------------------------------------------------------------------*\
                  Class morrisBoulayNormal Declaration
\*---------------------------------------------------------------------------*/

class morrisBoulayNormal : public normalViscosityModel
{

    protected:

    //- Protected data
    dictionary morrisBoulayNormalDict_;
    dimensionedScalar cMax_;
    dimensionedScalar nu0_;


    //- Disallow copy construct
    morrisBoulayNormal(const morrisBoulayNormal&);


    public:

    //- Runtime type information
    TypeName("morrisBoulayNormal");


    //- Construct from components
    morrisBoulayNormal
    (
        const fvMesh& mesh,
        const volVectorField& U,
        const volScalarField& c,
        const dictionary& dict
    );


    //- Destructor
    virtual ~morrisBoulayNormal()
    {}


    //- Member Functions
    virtual void updateNormalViscosity();

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
} // End namespace normalViscosityModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
