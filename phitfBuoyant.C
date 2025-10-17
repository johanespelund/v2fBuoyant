/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2012-2023 OpenFOAM Foundation
     \\/     M anipulation  |
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

#include "phitfBuoyant.H"
#include "fvModels.H"
#include "fvConstraints.H"
#include "bound.H"
#include "fvmSup.H"
#include "fvMatrices.H"
#include "fvcGrad.H"
#include "fvcLaplacian.H"
#include "wallDist.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace RASModels
{

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

template<class BasicMomentumTransportModel>
void phitfBuoyant<BasicMomentumTransportModel>::correctNut()
{
    this->nut_ = Cmu_*phit_*k_*T_;
    this->nut_.correctBoundaryConditions();
    fvConstraints::New(this->mesh_).constrain(this->nut_);
}


template<class BasicMomentumTransportModel>
tmp<volScalarField> phitfBuoyant<BasicMomentumTransportModel>::Ts() const
{
    // (LUU:Eq. 7)
    return
        max
        (
            k_/epsilon_,
            CT_*sqrt
            (
                max
                (
                    this->nu(),
                    dimensionedScalar(this->nu()().dimensions(), Zero)
                )/epsilon_
            )
        );
}

template<class BasicMomentumTransportModel>
tmp<volScalarField> phitfBuoyant<BasicMomentumTransportModel>::Ls() const
{
    // (LUU:Eq. 7)
    return
        CL_*max
        (
            pow(k_, 1.5)/epsilon_,
            Ceta_*pow025
            (
                pow3
                (
                    max
                    (
                        this->nu(),
                        dimensionedScalar(this->nu()().dimensions(), Zero)
                    )
                )/epsilon_
            )
        );
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class BasicMomentumTransportModel>
phitfBuoyant<BasicMomentumTransportModel>::phitfBuoyant
(
    const alphaField& alpha,
    const rhoField& rho,
    const volVectorField& U,
    const surfaceScalarField& alphaRhoPhi,
    const surfaceScalarField& phi,
    const viscosity& viscosity,
    const word& type
)
:
    eddyViscosity<RASModel<BasicMomentumTransportModel>>
    (
        type,
        alpha,
        rho,
        U,
        alphaRhoPhi,
        phi,
        viscosity
    ),
    v2fBase(),
    includeNu_
    (
        Switch::lookupOrAddToDict
        (
            "includeNu",
            this->coeffDict_,
            true
        )
    ),

    Cmu_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Cmu",
            this->coeffDict_,
            0.22
        )
    ),
    CmuKEps_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "CmuKEps",
            this->coeffDict_,
            0.09
        )
    ),
    Cf1_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Cf1",
            this->coeffDict_,
            1.4
        )
    ),
    Cf2_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Cf2",
            this->coeffDict_,
            0.3
        )
    ),
    C1_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "C1",
            this->coeffDict_,
            1.4
        )
    ),
    C2_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "C2",
            this->coeffDict_,
            0.3
        )
    ),
    CT_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "CT",
            this->coeffDict_,
            6.0
        )
    ),
    CL_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "CL",
            this->coeffDict_,
            0.25
        )
    ),
    Ceta_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Ceta",
            this->coeffDict_,
            110.0
        )
    ),
    Ceps1a_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Ceps1a",
            this->coeffDict_,
            1.4
        )
    ),
    Ceps1b_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Ceps1b",
            this->coeffDict_,
            1.0
        )
    ),
    Ceps1c_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Ceps1c",
            this->coeffDict_,
            0.05
        )
    ),
    Ceps2_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Ceps2",
            this->coeffDict_,
            1.9
        )
    ),
    Ceps3_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Ceps3",
            this->coeffDict_,
            -0.33
        )
    ),
    sigmaK_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "sigmaK",
            this->coeffDict_,
            1.0
        )
    ),
    sigmaEps_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "sigmaEps",
            this->coeffDict_,
            1.3
        )
    ),
    sigmaPhit_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "sigmaPhit",
            this->coeffDict_,
            1.0
        )
    ),
    Cg_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Cg",
            this->coeffDict_,
            1/0.9
        )
    ),
    Cphi_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Cphi",
            this->coeffDict_,
            0.3
        )
    ),
    THFM_
    (
            "GGDH"
    ),
    k_
    (
        IOobject
        (
            this->groupName("k"),
            this->runTime_.name(),
            this->mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh_
    ),
    epsilon_
    (
        IOobject
        (
            this->groupName("epsilon"),
            this->runTime_.name(),
            this->mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh_
    ),
    phit_
    (
        IOobject
        (
            this->groupName("phit"),
            this->runTime_.name(),
            this->mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh_
    ),
    f_
    (
        IOobject
        (
            this->groupName("f"),
            this->runTime_.name(),
            this->mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh_
    ),
    T_
    (
        IOobject
        (
            "T",
            this->runTime_.name(),
            this->mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        ),
        this->mesh_,
        dimensionedScalar(dimTime, Zero)
    ),
    g_
    (
      this->mesh_.objectRegistry::template
        lookupObject<uniformDimensionedVectorField>("g")
    ),
    phitMin_(dimensionedScalar(phit_.dimensions(), small)),
    fMin_(dimensionedScalar(f_.dimensions(), 0)),
    epsilonMin_(dimensionedScalar(epsilon_.dimensions(), small)),
    TMin_(dimensionedScalar("TMin", dimTime, small)),
    L2Min_(dimensionedScalar("L2Min", sqr(dimLength), small))
{
    this->coeffDict().readIfPresent("THFM", THFM_);
    Info << "  Turbulence heat flux model: " << THFM_ << nl;
    bound(k_, this->kMin_);
    bound(epsilon_, epsilonMin_);
    bound(phit_, phitMin_);
    bound(f_, fMin_);
    maxMagPb_ = 0;


    if (type == typeName)
    {
        this->printCoeffs(type);
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class BasicMomentumTransportModel>
bool phitfBuoyant<BasicMomentumTransportModel>::read()
{
    if (eddyViscosity<RASModel<BasicMomentumTransportModel>>::read())
    {
        Cmu_.readIfPresent(this->coeffDict());
        CmuKEps_.readIfPresent(this->coeffDict());
        C1_.readIfPresent(this->coeffDict());
        C2_.readIfPresent(this->coeffDict());
        CL_.readIfPresent(this->coeffDict());
        Ceta_.readIfPresent(this->coeffDict());
        Ceps2_.readIfPresent(this->coeffDict());
        Ceps3_.readIfPresent(this->coeffDict());
        sigmaK_.readIfPresent(this->coeffDict());
        sigmaEps_.readIfPresent(this->coeffDict());
        Cg_.readIfPresent(this->coeffDict());
        Cphi_.readIfPresent(this->coeffDict());

        // Read THFM is found in coeffDict

        return true;
    }
    else
    {
        return false;
    }
}


template<class BasicMomentumTransportModel>
void phitfBuoyant<BasicMomentumTransportModel>::correct()
{
    if (!this->turbulence_)
    {
        return;
    }

    // Local references
    const alphaField& alpha = this->alpha_;
    const rhoField& rho = this->rho_;
    const surfaceScalarField& alphaRhoPhi = this->alphaRhoPhi_;
    const volVectorField& U = this->U_;
    volScalarField& nut = this->nut_;
    const Foam::fvModels& fvModels(Foam::fvModels::New(this->mesh_));
    const Foam::fvConstraints& fvConstraints
    (
        Foam::fvConstraints::New(this->mesh_)
    );

    eddyViscosity<RASModel<BasicMomentumTransportModel>>::correct();

    volScalarField::Internal divU(fvc::div(fvc::absolute(this->phi(), U)));

    const dimensionedScalar epsSmall("epsSmall", epsilon_.dimensions(), SMALL);
    const dimensionedScalar k0("k0", k_.dimensions(), SMALL);

    const volVectorField gradRho(fvc::grad(rho)/rho);


    Info << "  Turbulence heat flux model: " << THFM_ << nl;
    volScalarField Pb = THFM_ == "SGDH" ?
      -nut*Cg_*g_ & gradRho :
      /* -(3/2)*Cg_*nut/(k_ + k0)*(this->sigma() & gradRho) & g_; */
      (3/2)*Cg_*Cphi_*(k_/(epsilon_ + epsSmall))*(this->sigma() & gradRho) & g_;

    scalar maxMagPbOld = maxMagPb_;
    volScalarField magPb = mag(Pb);
    maxMagPb_ = gMax(magPb);

    Info << "old max(mag(Pb)): " << maxMagPbOld << endl;
    Info << "new max(mag(Pb)): " << maxMagPb_ << endl;

    tmp<volSymmTensorField> tS(symm(fvc::grad(U)));
    volScalarField G(this->GName(), nut*(2.0*(dev(tS()) && tS())));
    tS.clear();

    T_ = Ts();
    bound(T_, TMin_);

    const volScalarField L2(typedName("L2"), sqr(Ls()) + L2Min_);

    const volScalarField::Internal Ceps1
    (
        typedName("Ceps1"),
        Ceps1a_*(Ceps1b_ + Ceps1c_*sqrt(1.0/phit_()))
    );

    const dimensionedScalar smallnum("smallNum", U.dimensions()*g_.dimensions(), SMALL);
    const volScalarField sinTheta = sqrt(1 - pow((U & g_)/(mag(U) * mag(g_) + smallnum),2)) * sign(U & g_);

    Info << "Min/Max sinTheta: " << gMin(sinTheta) << " / " << gMax(sinTheta) << endl;


    // Write Pk and Pb
    volScalarField Pk_
    (
        IOobject
        (
            "Pk",
            this->mesh_.time().name(),
            this->mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        G
    );
    volScalarField Pb_
    (
        IOobject
        (
            "Pb",
            this->mesh_.time().name(),
            this->mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        Pb
    );

    const volVectorField& wallNormal = wallDist::New(this->mesh_).n();
    volScalarField limiter = 1 *  mag(g_ & wallNormal) / (mag(g_) * mag(wallNormal));

    volScalarField limiter_
    (
        IOobject
        (
            "wall_limiter",
            this->mesh_.time().name(),
            this->mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        limiter
    );
    // volScalarField limiter
    // (
    //     IOobject("limiter", this->mesh_),
    //     2 * mag(g_ & wallNormal) / (mag(g_) * mag(wallNormal))
    // );

    // Pb *= 0;

    if (this->mesh_.time().write())
    {
        Pk_.write();
        Pb_.write();
        limiter_.write();
    }

    // Print min/max of G and Pb
    Info << "  G: min/max = "
         << gMin(G) << " / " << gMax(G) << nl
          << "  Pb: min/max = "
          << gMin(Pb) << " / " << gMax(Pb) << nl;

    // Print min/max of k, epsilon, v2 and f
    Info << "  k: min/max = "
         << gMin(k_) << " / " << gMax(k_) << nl
          << "  epsilon: min/max = "
          << gMin(epsilon_) << " / " << gMax(epsilon_) << nl
          << "  phi: min/max = "
          << gMin(phit_) << " / " << gMax(phit_) << nl
          << "  f: min/max = "
          << gMin(f_) << " / " << gMax(f_) << nl;
    // Update epsilon (and possibly G) at the wall
    epsilon_.boundaryFieldRef().updateCoeffs();


    // const volScalarField posPb = volScalarField("posPb", pos0(Pb), Pb.dimensions());
    // Info << "posPb dimnesion: " << posPb.dimensions() << endl;

    // Dissipation equation
    tmp<fvScalarMatrix> epsEqn
    (
        fvm::ddt(alpha, rho, epsilon_)
      + fvm::div(alphaRhoPhi, epsilon_)
      - fvm::laplacian(alpha*rho*DepsilonEff(), epsilon_)
      ==
        // alpha()*rho()*Ceps1*(G() + Pb())/T_()
       // alpha()*rho()*Ceps1*(max(Pb(),Pb()*0))*sinTheta/T_()
       // alpha()*rho()*Ceps1*(max(Pb(),Pb()*0))/T_()
       alpha()*rho()*Ceps1*Pb()/T_()
        // alpha()*rho()*Ceps1*posPb*sinTheta/T_()
      + alpha()*rho()*Ceps1*(G())/T_()
      - fvm::SuSp
        (
            (2.0/3.0*Ceps1)*(alpha()*rho()*divU),
            epsilon_
        )
      - fvm::Sp(alpha()*rho()*Ceps2_/T_(), epsilon_)
      + fvModels.source(alpha, rho, epsilon_)
    );

    epsEqn.ref().relax();
    fvConstraints.constrain(epsEqn.ref());
    epsEqn.ref().boundaryManipulate(epsilon_.boundaryFieldRef());
    solve(epsEqn);
    fvConstraints.constrain(epsilon_);
    bound(epsilon_, epsilonMin_);


    // Turbulent kinetic energy equation
    tmp<fvScalarMatrix> kEqn
    (
        fvm::ddt(alpha, rho, k_)
      + fvm::div(alphaRhoPhi, k_)
      - fvm::laplacian(alpha*rho*DkEff(), k_)
      ==
        // alpha()*rho()*(G() + Pb())
        alpha()*rho()*(G())
      - alpha()*rho()*fvm::SuSp(-Pb()/k_, k_)
      - fvm::SuSp(2.0/3.0*alpha()*rho()*divU, k_)
      - fvm::Sp(alpha()*rho()*(1.0/T_()), k_)
      + fvModels.source(alpha, rho, k_)
    );

    kEqn.ref().relax();
    fvConstraints.constrain(kEqn.ref());
    solve(kEqn);
    fvConstraints.constrain(k_);
    bound(k_, this->kMin_);


    // Relaxation function equation
    tmp<fvScalarMatrix> fEqn
    (
      - fvm::laplacian(f_)
     ==
      - fvm::Sp(1.0/L2(), f_)
      - (
            (Cf1_ - 1.0)*(phit_() - 2.0/3.0)/T_()
           -(Cf2_*(G()))/k_()
           +(Cf2_*(2.0/3.0)*divU)
           -(2.0*this->nu()*(fvc::grad(phit_) & fvc::grad(k_)))()/k_()
           -(this->nu()*fvc::laplacian(phit_))()
        )/L2()
    );


    fEqn.ref().relax();
    fvConstraints.constrain(fEqn.ref());
    solve(fEqn);
    fvConstraints.constrain(f_);
    /* bound(f_, fMin_); */

    // Turbulence stress normal to streamlines equation
    Info << "Solving phitEqn" << endl;
    tmp<fvScalarMatrix> phitEqn
    (
        fvm::ddt(alpha, rho, phit_)
      + fvm::div(alphaRhoPhi, phit_)
      - fvm::laplacian(alpha*rho*(this->nu() + this->nut_/sigmaPhit_), phit_)
      ==
        alpha()*rho()*f_()
      // + alpha()*rho()*fvm::SuSp((2 - phit_())*Pb()/k_(), phit_)
      // + alpha()*rho()*(limiter - phit_)*Pb()*phit_()/k()
      // + alpha()*rho()*fvm::SuSp(Pb()/k_(), phit_)
      // + alpha()*rho()*(limiter - phit_())*(Pb()/k_())*phit_
      // - alpha()*rho()*Pb()/k_()*phit_
      // - fvm::SuSp(alpha()*rho()*Pb()/k_(),phit_)
      // - fvm::SuSp(alpha()*rho()*(limiter - phit_())*Pb()/k_()/phit_(),phit_)
      // - fvm::SuSp(alpha()*rho()*(1 - phit_())*Pb()/k_(),phit_)
      - fvm::SuSp
        (
            alpha()*rho()*
            (
                // (G() + Pb())/k_()
                (max(Pb(),Pb()*0))/k_()
              + (G())/k_()
              - (2.0/3.0)*divU
              - (2.0*nut*(fvc::grad(phit_) & fvc::grad(k_)))()
                /(k_()*sigmaPhit_*phit_())
            )
          , phit_
        )
      + fvModels.source(alpha, rho, phit_)
    );
    phitEqn.ref().relax();
    fvConstraints.constrain(phitEqn.ref());
    solve(phitEqn);
    fvConstraints.constrain(phit_);
    bound(phit_, phitMin_);

    correctNut();


}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace Foam

// ************************************************************************* //
