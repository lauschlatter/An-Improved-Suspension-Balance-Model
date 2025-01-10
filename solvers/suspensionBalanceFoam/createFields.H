Info<< "Reading transportProperties\n" << endl;

IOdictionary transportProperties
(
    IOobject
    (
        "transportProperties",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    )
);

Info<< "Reading field p\n" << endl;
volScalarField p
(
    IOobject
    (
        "p",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info<< "Reading field U\n" << endl;
volVectorField U
(
    IOobject
    (
        "U",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info<< "Reading field c\n" << endl;
volScalarField c
(
    IOobject
    (
        "c",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

#include "newCreatePhi.H"

label pRefCell = 0;
scalar pRefValue = 0.0;
setRefCell(p, mesh.solutionDict().subDict("PISO"), pRefCell, pRefValue);

bool fixCAve =
    mesh.solutionDict().subDict("PISO").lookupOrDefault("fixCAve", false);
bool fixCflux =
    mesh.solutionDict().subDict("PISO").lookupOrDefault("fixCflux", false);
scalar cAve = mesh.solutionDict().subDict("PISO").lookupOrDefault("cAve", 0.0);
word patchName = "inlet";
mesh.solutionDict().subDict("PISO").readIfPresent("inletPatch", patchName);
label inletPatchID = mesh.boundaryMesh().findPatchID(patchName);

suspensionModel mixture(mesh, transportProperties, p, U, c, phi);
mixture.update();

#include "createFvOptions.H"
