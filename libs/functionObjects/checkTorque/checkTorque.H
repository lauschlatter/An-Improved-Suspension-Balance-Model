/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright held by original author
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

Class
    checkTorque
    
Author
    LTFD.
    Laboratory of Thermo-Fluid Dynamics (LTFD)
    Chemical Engineering Program (PEQ)
    Alberto Luiz Coimbra Institute for Graduate Studies and Research in Engineering (COPPE)
    Federal University of Rio de Janeiro Rio de Janeiro, Brazil
    
Written by:
    Lauren Schlatter Fernandes
    laurenfernandes@peq.coppe.ufrj.br

    Gabriel Gonçalves da Silva Ferreira
    gferreira@peq.coppe.ufrj.br

SourceFiles
    checkTorque.C

Description
    FunctionObject that calculates and reports the torque on a given patch.
    To use, add the following to controlDict:    
    functions
    ( 
        history
        {
            type checkTorque;
            CofR (x y z);
            patches (patchName);
            functionObjectLibs
            (
                "libcheckTorque.so"
            );
        }
    );

\*---------------------------------------------------------------------------*/

#ifndef checkTorque_H
#define checkTorque_H

#include "functionObject.H"
#include "dictionary.H"
#include "fvMesh.H"
#include "OFstream.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                        Class checkTorque Declaration
\*---------------------------------------------------------------------------*/

class checkTorque : public functionObject
{
    // Private data
    const word name_;
    const Time& time_;
    word regionName_;
    vector CofR_;
    labelHashSet patchSet_;
    OFstream* historyFilePtr_;


// Disallow default bitwise copy construct
    checkTorque
    (
        const checkTorque&
    );


// Disallow default bitwise assignment
    void operator=(const checkTorque&);


    public:

// Runtime type information
    TypeName("checkTorque");


// Construct from components
    checkTorque
    (
        const word& name,
        const Time&,
        const dictionary&
    );


    // Member Functions
    virtual bool start();
    virtual bool execute();
    virtual bool write();
    virtual bool read(const dictionary& dict);

};

} // End namespace functionObject
} // End namespace Foam

#endif

// ************************************************************************* //
