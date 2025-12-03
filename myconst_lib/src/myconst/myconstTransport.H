/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2020 OpenCFD Ltd.
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
    Foam::constTransport

Group
    grpSpecieTransport

Description
    Constant properties Transport package.
    Templated into a given thermodynamics package (needed for thermal
    conductivity).

SourceFiles
    constTransportI.H
    constTransport.C

\*---------------------------------------------------------------------------*/

#ifndef myconstTransport_H
#define myconstTransport_H

#include "scalar.H"
#include "dictionary.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations

template<class Thermo> class myconstTransport;

template<class Thermo>
inline myconstTransport<Thermo> operator+
(
    const myconstTransport<Thermo>&,
    const myconstTransport<Thermo>&
);

template<class Thermo>
inline myconstTransport<Thermo> operator*
(
    const scalar,
    const myconstTransport<Thermo>&
);

template<class Thermo>
Ostream& operator<<
(
    Ostream&,
    const myconstTransport<Thermo>&
);


/*---------------------------------------------------------------------------*\
                           Class constTransport Declaration
\*---------------------------------------------------------------------------*/

template<class Thermo>
class myconstTransport
:
    public Thermo
{
    // Private Data

        //- Constant dynamic viscosity [Pa.s]
        //scalar mu_;

        //- Reciprocal Prandtl Number [-]
        //scalar rPr_;

        //- Constant dynamic viscosity [Pa.s]
        scalar mu_;
        scalar alpha_;
        scalar T0_;

        //- Is the Prandtl Number constant?
        bool constPr_;

        //- Reciprocal Prandtl Number [-]
        scalar rPr_;

        //- Constant thermal conductivity [W/m/K]
        scalar kappa_;


    // Private Member Functions


    // Private Member Functions

        //- Construct from components
        inline myconstTransport
        (
            const Thermo& t,
            const scalar mu,
            const scalar alpha,
            const scalar T0,
            const bool constPr,
            const scalar rPr,
            const scalar kappa
        );


public:

    // Generated Methods: copy construct, copy assignment


    // Constructors

        //- Construct as named copy
        inline myconstTransport(const word&, const myconstTransport&);

        //- Construct from dictionary
        explicit myconstTransport(const dictionary& dict);

        //- Construct and return a clone
        inline autoPtr<myconstTransport> clone() const;

        // Selector from dictionary
        inline static autoPtr<myconstTransport> New(const dictionary& dict);


    // Member Functions

        //- Return the instantiated type name
        static word typeName()
        {
            return "myconst<" + Thermo::typeName() + '>';
        }

        //- Dynamic viscosity [kg/ms]
        inline scalar mu(const scalar p, const scalar T) const;

        //- Thermal conductivity [W/mK]
        inline scalar kappa(const scalar p, const scalar T) const;
        
        inline scalar Pr(const scalar p, const scalar T) const;
        inline scalar T0(const scalar p, const scalar T) const;

        //- Thermal diffusivity of enthalpy [kg/ms]
        inline scalar alpha(const scalar p, const scalar T) const;

        // Species diffusivity
        //inline scalar D(const scalar p, const scalar T) const;

        //- Write to Ostream
        void write(Ostream& os) const;


    // Member Operators

        inline void operator+=(const myconstTransport&);

        inline void operator*=(const scalar);


    // Friend Operators

        friend myconstTransport operator+ <Thermo>
        (
            const myconstTransport&,
            const myconstTransport&
        );

        friend myconstTransport operator* <Thermo>
        (
            const scalar,
            const myconstTransport&
        );


    // IOstream Operators

        friend Ostream& operator<< <Thermo>
        (
            Ostream&,
            const myconstTransport&
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "myconstTransportI.H"

#ifdef NoRepository
    #include "myconstTransport.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
