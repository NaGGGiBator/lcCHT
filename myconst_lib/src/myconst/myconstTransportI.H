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

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Thermo>
inline Foam::myconstTransport<Thermo>::myconstTransport
(
    const Thermo& t,
    const scalar mu,
    const scalar alpha,
    const scalar T0,
    const bool constPr,
    const scalar rPr,
    const scalar kappa
)
:
    Thermo(t),
    mu_(mu),
    alpha_(alpha),
    T0_(T0),
    constPr_(constPr),
    rPr_(rPr),
    kappa_(kappa)
{}


template<class Thermo>
inline Foam::myconstTransport<Thermo>::myconstTransport
(
    const word& name,
    const myconstTransport& ct
)
:
    Thermo(name, ct),
    mu_(ct.mu_),
    alpha_(ct.alpha_),
    T0_(ct.T0_),
    constPr_(ct.constPr_),
    rPr_(ct.rPr_),
    kappa_(ct.kappa_)
{}


template<class Thermo>
inline Foam::autoPtr<Foam::myconstTransport<Thermo>>
Foam::myconstTransport<Thermo>::clone() const
{
    return autoPtr<myconstTransport<Thermo>>::New(*this);
}


template<class Thermo>
inline Foam::autoPtr<Foam::myconstTransport<Thermo>>
Foam::myconstTransport<Thermo>::New
(
    const dictionary& dict
)
{
    return autoPtr<myconstTransport<Thermo>>::New(dict);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Thermo>
inline Foam::scalar Foam::myconstTransport<Thermo>::mu
(
    const scalar p,
    const scalar T
) const
{
    return mu_ * pow((T / T0_), alpha_);
}


template<class Thermo>
inline Foam::scalar Foam::myconstTransport<Thermo>::kappa
(
    const scalar p,
    const scalar T
) const
{
    return constPr_ ? this->Cp(p, T) * mu_ * pow((T / T0_), alpha_) * rPr_ : kappa_;
}



// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class Thermo>
inline void Foam::myconstTransport<Thermo>::operator+=
(
    const myconstTransport<Thermo>& st
)
{
    scalar Y1 = this->Y();

    Thermo::operator+=(st);

    if (mag(this->Y()) > SMALL)
    {
        if
        (
            myconstTransport<Thermo>::debug && (constPr_ != st.constPr_)
        )
        {
            FatalErrorInFunction
                << "Constant " << (constPr_ ? "Pr" : "kappa") << " for "
                << (this->name().size() ? this->name() : "others") << " but "
                << "constant " << (st.constPr_ ? "Pr" : "kappa") << " for "
                << (st.name().size() ? st.name() : "others")
                << exit(FatalError);
        }
        Y1 /= this->Y();
        scalar Y2 = st.Y()/this->Y();

        mu_ = Y1*mu_ + Y2*st.mu_;
        rPr_ = constPr_ ? 1/(Y1/rPr_ + Y2/st.rPr_) : NaN;
        kappa_ = constPr_ ? NaN : Y1*kappa_ + Y2*st.kappa_;
    }
}


template<class Thermo>
inline void Foam::myconstTransport<Thermo>::operator*=
(
    const scalar s
)
{
    Thermo::operator*=(s);
}


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

template<class Thermo>
inline Foam::myconstTransport<Thermo> Foam::operator+
(
    const myconstTransport<Thermo>& ct1,
    const myconstTransport<Thermo>& ct2
)
{
    Thermo t
    (
        static_cast<const Thermo&>(ct1) + static_cast<const Thermo&>(ct2)
    );

    if (mag(t.Y()) < SMALL)
    {
        return myconstTransport<Thermo>
        (
            t,
            0,
            ct1.constPr_,
            ct1.rPr_
        );
    }
    else
    {
        scalar Y1 = ct1.Y()/t.Y();
        scalar Y2 = ct2.Y()/t.Y();

        if
        (
            myconstTransport<Thermo>::debug
         && (ct1.constPr_ != ct2.constPr_)
        )
        {
            FatalErrorInFunction
                << "Constant " << (ct1.constPr_ ? "Pr" : "kappa") << " for "
                << (ct1.name().size() ? ct1.name() : "others") << " but "
                << "constant " << (ct2.constPr_ ? "Pr" : "kappa") << " for "
                << (ct2.name().size() ? ct2.name() : "others")
                << exit(FatalError);
        }

        return myconstTransport<Thermo>
        (
            t,
            Y1*ct1.mu_ + Y2*ct2.mu_,
            ct1.constPr_,
            ct1.constPr_ ? 1/(Y1/ct1.rPr_ + Y2/ct2.rPr_) : NaN,
            ct1.constPr_ ? NaN : Y1*ct1.kappa_ + Y2*ct2.kappa_
        );
    }
}


template<class Thermo>
inline Foam::myconstTransport<Thermo> Foam::operator*
(
    const scalar s,
    const myconstTransport<Thermo>& ct
)
{
    return myconstTransport<Thermo>
    (
        s*static_cast<const Thermo&>(ct),
        ct.mu_,
        ct.alpha_,
        ct.T0_,
        ct.constPr_,
        ct.rPr_,
        ct.kappa_
    );
}


// ************************************************************************* //
