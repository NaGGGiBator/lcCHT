/*---------------------------------------------------------------------------*\

\*---------------------------------------------------------------------------*/

#include "lcCHT.H"
#include "fluidThermo.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(lcCHT, 0);
    addToRunTimeSelectionTable(functionObject, lcCHT, dictionary);
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

double Foam::functionObjects::lcCHT::findPatchID(word patchName)
{
    Info << endl << endl << endl;
    Info << patchName << endl;
    Info << endl << endl << endl;

    const polyBoundaryMesh& patches = mesh_.boundaryMesh();

    auto patchID = -1;

    forAll(patches, i)
    {
        if (patches[i].name() == patchName)
        {
            patchID = i;
            break;
        }
    }

    if (patchID == -1)
    {
        FatalErrorInFunction << "Patch " << patchName << " not found!" << exit(FatalError);
    }

    return 0;
}

bool Foam::functionObjects::lcCHT::calc()
{
    Info << endl;
    Info << "************************************************************" << endl;
    Info << "******************************  " << runTime_.timeName() << "  ******************************" << endl;
    Info << "************************************************************" << endl;
    Info << endl;

    const volScalarField& T = lookupObject<volScalarField>(fieldName_);

    findPatchID(patchName);


//    return store
//    (
//        resultName_,
//        T + dimensionedScalar("Tinc", T.dimensions(), 1000)
//    );

    return true; //false;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::lcCHT::lcCHT
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fieldExpression(name, runTime, dict, "T"),
    runTime_(runTime)
{
    time_step_fluid = readScalar(dict.lookup("time_step_fluid"));
    time_step_solid = readScalar(dict.lookup("time_step_solid"));
    patchName = word(dict.lookup("patchName"));

    //setResultName("T2", "T");
}


// ************************************************************************* //

