Info<< "Reading/calculating face flux field phi\n\n" << endl;

#include "fvsPatchField.H"

wordList phiTypes
(
    U.boundaryField().size(),
    calculatedFvPatchScalarField::typeName
);

wordList UTypes
(
    U.boundaryField().size()
);

forAll(U.boundaryField(), i)
{
    UTypes[i] = U.boundaryField()[i].type();

    if (isA<fixedValueFvPatchVectorField>(U.boundaryField()[i]))
    {
        phiTypes[i] = fixedValueFvPatchScalarField::typeName;
    }
}

surfaceScalarField phi
(
    IOobject
    (
        "phi",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    fvc::flux(U),
    phiTypes
);
