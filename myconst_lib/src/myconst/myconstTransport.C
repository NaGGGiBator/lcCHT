/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
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

\*---------------------------------------------------------------------------*/

#include "myconstTransport.H"
#include "IOstreams.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Thermo>
Foam::myconstTransport<Thermo>::myconstTransport(const word& name, const dictionary& dict)
:
    Thermo(name, dict)
{
    const dictionary& transportDict = dict.subDict("transport");

    mu_ = transportDict.lookup<scalar>("mu0");
    alpha_ = transportDict.lookup<scalar>("alpha");
    T0_ = transportDict.lookup<scalar>("T0");

    const bool foundPr = transportDict.found("Pr");
    const bool foundKappa = transportDict.found("kappa");

    if (foundPr == foundKappa)
    {
        FatalIOErrorInFunction(dict)
            << "Either Pr or kappa must be specified, but not both."
            << exit(FatalIOError);
    }

    constPr_ = foundPr;
    rPr_ = constPr_ ? 1/transportDict.lookup<scalar>("Pr") : NaN;
    kappa_ = constPr_ ? NaN : transportDict.lookup<scalar>("kappa");
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Thermo>
void Foam::myconstTransport<Thermo>::myconstTransport::write(Ostream& os) const
{
   
    os.beginBlock(this->name());

    Thermo::write(os);

    // Entries in dictionary format
    {
        os.beginBlock("transport");
        os.writeEntry("mu0", mu_);
        os.writeEntry("alpha", alpha_);
        os.writeEntry("T0", T0_);
        os.endBlock();

        if (constPr_)
        {
            os.writeEntry("Pr", 1.0/rPr_);
        }
        else
        {
            os.writeEntry("kappa", kappa_);
        }
    }

    os.endBlock();

}


// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

template<class Thermo>
Foam::Ostream& Foam::operator<<(Ostream& os, const myconstTransport<Thermo>& ct)
{
    ct.write(os);
    return os;
}


// ************************************************************************* //
