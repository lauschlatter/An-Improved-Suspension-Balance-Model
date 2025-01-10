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
    Foam::anisotropicStressTensorModel

Description
    anisotropicStressTensorModel module

\*---------------------------------------------------------------------------*/

#ifndef anisotropicStressTensorModel_H
#define anisotropicStressTensorModel_H

#include "dictionary.H"
#include "volFieldsFwd.H"
#include "volFields.H"
#include "runTimeSelectionTables.H"
#include "surfaceInterpolate.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class anisotropicStressTensorModel Declaration
\*---------------------------------------------------------------------------*/

class anisotropicStressTensorModel
{

    protected:

    //- Protected data
    const fvMesh& mesh_;
    const volVectorField& U_;
    const volScalarField& c_;
    volTensorField Q_;


    //- Disallow copy construct
    anisotropicStressTensorModel(const anisotropicStressTensorModel&);


    public:

    //- Runtime type information
    TypeName("anisotropicStressTensorModel");

    //- Declare run-time constructor selection table
    declareRunTimeSelectionTable
    (
        autoPtr,
        anisotropicStressTensorModel,
        dictionary,
        (
            const fvMesh& mesh,
            const volVectorField& U,
            const volScalarField& c,
            const dictionary& dict
        ),
        (mesh, U, c, dict)
    );


    //- Return a reference to the selected viscosity model
    static autoPtr<anisotropicStressTensorModel> New
    (
        const fvMesh& mesh,
        const volVectorField& U,
        const volScalarField& c,
        const dictionary& dict
    );


    //- Construct from components
    anisotropicStressTensorModel
    (
        const fvMesh& mesh,
        const volVectorField& U,
        const volScalarField& c,
        const dictionary& dict
    );


    //- Destructor
    virtual ~anisotropicStressTensorModel()
    {}


    //- Member Functions
    virtual void updateAnisotropicStressTensor() = 0;
    
    const volTensorField& Q() const
    {
        return Q_;
    }
    volTensorField& Q()
    {
        return Q_;
    }

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
