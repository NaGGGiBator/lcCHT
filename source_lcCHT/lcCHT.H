/*---------------------------------------------------------------------------*\


\*---------------------------------------------------------------------------*/

#ifndef functionObjects_lcCHT_H
#define functionObjects_lcCHT_H

#include "fieldExpression.H"

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

    const Time& runTime_;
    scalar  time_step_solid;
    scalar  time_step_fluid;
    scalar  iteration_fluid;
    word    patchName;

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
