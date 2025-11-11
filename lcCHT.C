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

bool Foam::functionObjects::lcCHT::calc()
{
    Info << endl;
    Info << "************************************************************" << endl;
    Info << "******************************  " << runTime_.timeName() << "  ******************************" << endl;
    Info << "************************************************************" << endl;
    Info << endl;

    const volScalarField& T = 
        lookupObject<volScalarField>(fieldName_);

    return store
    (
        resultName_,
        T + dimensionedScalar("Tinc", T.dimensions(), 1000)
    );
    //}

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
    setResultName("T2", "T");
}


// ************************************************************************* //
