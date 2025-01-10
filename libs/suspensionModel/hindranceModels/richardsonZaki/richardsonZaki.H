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
    Foam::richardsonZaki

SourceFiles
    richardsonZaki.C
    
Description
    Foam::richardsonZaki

\*---------------------------------------------------------------------------*/

#ifndef richardsonZaki_H
#define richardsonZaki_H

#include "volFields.H"
#include "hindranceModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace hindranceModels
{
/*---------------------------------------------------------------------------*\
                      Class richardsonZaki Declaration
\*---------------------------------------------------------------------------*/

class richardsonZaki : public hindranceModel
{

    protected:

    //- Protected data
    dictionary richardsonZakiDict_;
    dimensionedScalar alpha_;
    dimensionedScalar cMax_;


    //- Disallow copy construct
    richardsonZaki(const richardsonZaki&);


    public:

    //- Runtime type information
    TypeName("richardsonZaki");


    //- Construct from components
    richardsonZaki
    (
        const fvMesh& mesh,
        const volVectorField& U,
        const volScalarField& c,
        const dictionary& dict
    );


    //- Destructor
    virtual ~richardsonZaki()
    {}


    //- Member Functions
    virtual void updateHindrance();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
} // End namespace hindranceModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
