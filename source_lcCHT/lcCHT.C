/*---------------------------------------------------------------------------*\

\*---------------------------------------------------------------------------*/

#include "lcCHT.H"
#include "fluidThermo.H"
#include "addToRunTimeSelectionTable.H"
#include "rhoReactionThermo.H"
#include "CombustionModel.H"
#include "basicChemistryModel.H"
#include "IOmanip.H" 
#include "Pstream.H" 

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

    return patchID;
}

double Foam::functionObjects::lcCHT::changeTime()
{
    const volScalarField& T = lookupObject<volScalarField>(fieldName_);

    auto patchID = findPatchID(patchName);

    const scalarField& areaMag      = mesh_.boundary()[patchID].magSf(); 
    const scalarField& Tp           = T.boundaryField()[patchID];


    if(firstIter)
    {
        T_old = Tp;
        firstIter = false;
    }
    else
    {
        scalarField deltaT(Tp.size());
        scalar T_old_sum = 0;
        
        forAll(Tp, i)
        {
            deltaT[i] = Tp[i] - T_old[i];
            T_old_sum += T_old[i];
        }

        scalar upperSqrIntegral     = gSum(deltaT * deltaT * areaMag);
        scalar downSqrIntegral      = gSum(Tp * Tp * areaMag);
        scalar areaSum              = gSum(areaMag);
        // scalar rn                   = sqrt(upperSqrIntegral / areaSum / downSqrIntegral /
        //                             areaSum);

        scalar rn                   = sqrt(upperSqrIntegral / areaSum);

        // Info << "----------------------------Your new rn = " << rn << nl;

        if(rn != 0)
        {
            newTimeStep_ = pow((1 / rn), 0.175) * newTimeStep_;
        }
        
        Info << endl;
        Info << "************************************************************" << endl;
        Info << "*********************** " << fixed << setprecision(10) << " rn = " << rn << " STEP = " << newTimeStep_ << " ***********************" << endl;
        Info << "************************************************************" << endl;
        Info << endl;

        T_old = Tp;
    }

    return newTimeStep_;
}


bool Foam::functionObjects::lcCHT::calc()
{
    Foam::Time& runTime = const_cast<Foam::Time&>(obr_.time()); // only to change the time step (because by default runTime_ is const)
    Info << "*********************** " << fixed << setprecision(10) << " rn = "  << " STEP = " << print_timeStep << " ***********************" << endl;  
    if(calc_fluid)
    {
        runTime.setDeltaT(time_step_fluid);

        IOdictionary fvSolFile
        (
            IOobject
            (
                "fvSolution",
                mesh_.time().system(),
                mesh_,
                IOobject::MUST_READ,
                IOobject::AUTO_WRITE
            )
        );

        dictionary& solvers = fvSolFile.subDict("solvers");
        solvers.subDict("p_rgh").set<label>("maxIter", 1000);
        solvers.subDict("p_rghFinal").set<label>("maxIter", 1000);
        solvers.subDict("U").set<label>("maxIter", 1000);
        solvers.subDict("h").set<label>("maxIter", 1000);
        solvers.subDict("R").set<label>("maxIter", 1000);
        solvers.subDict("Yi").set<label>("maxIter", 1000);
        solvers.subDict("UFinal").set<label>("maxIter", 1000);
        solvers.subDict("hFinal").set<label>("maxIter", 1000);
        solvers.subDict("RFinal").set<label>("maxIter", 1000);
        solvers.subDict("YiFinal").set<label>("maxIter", 1000);

        dictionary& pimpleDict = fvSolFile.subDict("PIMPLE");
        pimpleDict.set<label>("nCorrectors", 3);
        pimpleDict.set<bool>("frozenFlow", false);

        fvSolFile.regIOobject::write();
        fvSolution& fvSol = const_cast<fvSolution&>
        (
            mesh_.thisDb().lookupObject<fvSolution>("fvSolution")
        );

        fvSol.read();

        ++iterator_for_fluid;

        if(iterator_for_fluid >= iteration_fluid)
        {
            calc_fluid = false;
            iterator_for_fluid = 0;
            runTime.setDeltaT(time_step_solid);
            iteration_solid = std::ceil(changeTime() / time_step_solid);
            print_timeStep = changeTime();

            IOdictionary fvSolFile
            (
                IOobject
                (
                "fvSolution",
                mesh_.time().system(),
                mesh_,
                IOobject::MUST_READ,
                IOobject::AUTO_WRITE
                )
            );

            dictionary& solvers = fvSolFile.subDict("solvers");
            solvers.subDict("p_rgh").set<label>("maxIter", 0);
            solvers.subDict("p_rghFinal").set<label>("maxIter", 0);
            solvers.subDict("U").set<label>("maxIter", 0);
            solvers.subDict("h").set<label>("maxIter", 0);
            solvers.subDict("R").set<label>("maxIter", 0);
            solvers.subDict("Yi").set<label>("maxIter", 0);
            solvers.subDict("UFinal").set<label>("maxIter", 0);
            solvers.subDict("hFinal").set<label>("maxIter", 0);
            solvers.subDict("RFinal").set<label>("maxIter", 0);
            solvers.subDict("YiFinal").set<label>("maxIter", 0);

            dictionary& pimpleDict = fvSolFile.subDict("PIMPLE");
            pimpleDict.set<label>("nCorrectors", 0);
            pimpleDict.set<bool>("frozenFlow", true);

            fvSolFile.regIOobject::write();
            
            fvSolution& fvSol = const_cast<fvSolution&>
            (
                mesh_.thisDb().lookupObject<fvSolution>("fvSolution")
            );

            fvSol.read();

        }

        // const fileName chemDictPath = mesh_.time().constant()/"chemistryProperties";
        // IOdictionary chemDict
        // (
        //     IOobject
        //     (
        //         "chemistryProperties",
        //         mesh_.time().constant(),
        //         mesh_,
        //         IOobject::MUST_READ_IF_MODIFIED,
        //         IOobject::AUTO_WRITE
        //     )
        // );

        // chemDict.set("chemistry", "on1");
        // chemDict.regIOobject::write();
    } 
    else
    {
        --iteration_solid;

        if(iteration_solid <= 0)
        {
            calc_fluid = true;
        }
    }



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
    runTime_(runTime),
    calc_fluid(true),
    firstIter(true),
    iterator_for_fluid(0),
    T_old(0),
    print_timeStep(0)
{
    time_step_fluid     = readScalar(dict.lookup("time_step_fluid"));
    time_step_solid     = readScalar(dict.lookup("time_step_solid"));
    patchName           = word(dict.lookup("patchName"));
    iteration_fluid     = readScalar(dict.lookup("iteration_fluid"));
    newTimeStep_        = readScalar(dict.lookup("startTime"));

    //setResultName("T2", "T");
}


// ************************************************************************* //

