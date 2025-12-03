/*---------------------------------------------------------------------------*\


\*---------------------------------------------------------------------------*/

#ifndef functionObjects_lcCHT_H
#define functionObjects_lcCHT_H

#include "fieldExpression.H"
#include "PtrList.H"
#include "CombustionModel.H"
#include "rhoReactionThermo.H"
#include "compressibleTurbulenceModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                          Class clCHT Declaration
\*---------------------------------------------------------------------------*/

class lcCHT : public fieldExpression
{
private:

    virtual bool calc();
    double  findPatchID(word patchName);
    double  changeTime();

    const Time& runTime_;
    scalar      time_step_solid;
    scalar      time_step_fluid;
    scalar      iteration_fluid;
    scalar      iterator_for_fluid;
    scalar      iteration_solid;
    scalarField T_old;
    word        patchName;
    bool        calc_fluid;
    bool        firstIter;
    double      newTimeStep_;
    double      print_timeStep;


public:

    TypeName("lcCHT");

    lcCHT
    (
        const word& name,
        const Time& runTime,
        const dictionary& dict
    );

    lcCHT(const lcCHT&) = delete;

    void operator=(const lcCHT&) = delete;

    virtual ~lcCHT() = default;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
