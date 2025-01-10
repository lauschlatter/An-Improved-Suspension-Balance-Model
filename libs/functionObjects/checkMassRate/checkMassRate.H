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
    checkMassRate

Description
    FunctionObject reports the massRate for the domain entrance and exit
    
    Add the following to controlDict:
    
    functions
    ( 
       history
       {
           type checkMassRate;
           inletPatch   INLET;
           outletPatch  OUTLET;
           functionObjectLibs
           (
             "libcheckMassRate.so"
           );
       }
    );

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
    checkMassRate.C    

\*---------------------------------------------------------------------------*/

#ifndef checkMassRate_H
#define checkMassRate_H

#include "functionObject.H"
#include "dictionary.H"
#include "fvMesh.H"
#include "OFstream.H"
#include "suspensionModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                        Class checkMassRate Declaration
\*---------------------------------------------------------------------------*/

class checkMassRate
:
    public functionObject
{
// Private data
    const word name_;
    const Time& time_;
    word regionName_;
    label inletPatchID_;
    label outletPatchID_;
    OFstream* historyFilePtr_;


// Disallow default bitwise copy construct
    checkMassRate
    (
        const checkMassRate&
    );


// Disallow default bitwise assignment
    void operator=(const checkMassRate&);


    public:

// Runtime type information
    TypeName("checkMassRate");


// Construct from components
    checkMassRate
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
