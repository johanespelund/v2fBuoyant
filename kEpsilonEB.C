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

#include "kEpsilonEB.H"
#include "fvModels.H"
#include "fvConstraints.H"
#include "bound.H"
#include "fvmSup.H"
#include "fvcDdt.H"
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
void kEpsilonEB<BasicMomentumTransportModel>::correctNut()
{
    volScalarField magS(sqrt(2.0)*mag(symm(fvc::grad(this->U_))));

    this->nut_ = Cmu_*phit_*k_*
        min
        (
            T_,
            scalar(1.0)/
            max
            (
                dimensionedScalar(pow(dimTime,-1),VSMALL),
                Cmu_.value()*sqrt(3.0)*phit_*magS

            )
        );
    this->nut_.correctBoundaryConditions();
    fvConstraints::New(this->mesh_).constrain(this->nut_);
}


// Compute the limit turbulent time scale (TLLP:Eq.14)
template<class BasicMomentumTransportModel>
tmp<volScalarField> kEpsilonEB<BasicMomentumTransportModel>::Ts() const
{
    return
        sqrt
        (
            sqr(k_/epsilon_) + sqr(Ct_)*
            max
            (
                this->nu()/epsilon_,
                dimensionedScalar(sqr(dimTime), Zero)
            )
        );
}

// Compute the turbulent length scale (TLLP:Eq.12)
template<class BasicMomentumTransportModel>
tmp<volScalarField> kEpsilonEB<BasicMomentumTransportModel>::Ls() const
{
    return
        CL_*sqrt
        (
            max
            (
                pow3(k_)/sqr(epsilon_),
                dimensionedScalar(sqr(dimLength), Zero)
            )
            + sqr(Ceta_)*
            sqrt
            (
                max
                (
                    pow3(this->nu())/epsilon_,
                    dimensionedScalar(pow(dimLength,4), Zero)
                )
            )
        );
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class BasicMomentumTransportModel>
kEpsilonEB<BasicMomentumTransportModel>::kEpsilonEB
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
    // v2fBase(),
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
        dimensionedScalar::lookupOrAddToDict
        (
            "Cmu",
            this->coeffDict_,
            0.22
        )
    ),
    C1_
    (
        dimensionedScalar::lookupOrAddToDict
        (
            "C1",
            this->coeffDict_,
            1.7
        )
    ),
    C3_
    (
        dimensionedScalar::lookupOrAddToDict
        (
            "C3",
            this->coeffDict_,
            0.8
        )
    ),
    C4_
    (
        dimensionedScalar::lookupOrAddToDict
        (
            "C4",
            this->coeffDict_,
            0.625
        )
    ),
    C5_
    (
        dimensionedScalar::lookupOrAddToDict
        (
            "C5",
            this->coeffDict_,
            0.2
        )
    ),
    C1s_
    (
        dimensionedScalar::lookupOrAddToDict
        (
            "C1s",
            this->coeffDict_,
            0.9
        )
    ),
    C3s_
    (
        dimensionedScalar::lookupOrAddToDict
        (
            "C3s",
            this->coeffDict_,
            0.65
        )
    ),
    CK_
    (
        dimensionedScalar::lookupOrAddToDict
        (
            "CK",
            this->coeffDict_,
            2.3
        )
    ),
    Ct_
    (
        dimensionedScalar::lookupOrAddToDict
        (
            "Ct",
            this->coeffDict_,
            4.0
        )
    ),
    Ceps1_
    (
        dimensionedScalar::lookupOrAddToDict
        (
            "Ceps1_",
            this->coeffDict_,
            1.44
        )
    ),
    Ceps2_
    (
        dimensionedScalar::lookupOrAddToDict
        (
            "Ceps2",
            this->coeffDict_,
            1.9
        )
    ),
    CL_
    (
        dimensionedScalar::lookupOrAddToDict
        (
            "CL",
            this->coeffDict_,
            0.164
        )
    ),
    Ceta_
    (
        dimensionedScalar::lookupOrAddToDict
        (
            "Ceta",
            this->coeffDict_,
            75.0
        )
    ),
    sigmaK_
    (
        dimensionedScalar::lookupOrAddToDict
        (
            "sigmaK",
            this->coeffDict_,
            1.0
        )
    ),
    sigmaEps_
    (
        dimensionedScalar::lookupOrAddToDict
        (
            "sigmaEps",
            this->coeffDict_,
            1.2
        )
    ),
    sigmaPhit_
    (
        dimensionedScalar::lookupOrAddToDict
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
    ebf_
    (
        IOobject
        (
            this->groupName("ebf"),
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
    ebfMin_(dimensionedScalar("ebfMin", ebf_.dimensions(), Zero)),
    epsilonMin_(dimensionedScalar(epsilon_.dimensions(), small)),
    // TMin_(dimensionedScalar("TMin", dimTime, small)),
    L2Min_(dimensionedScalar("L2Min", sqr(dimLength), small))
{
    this->coeffDict().readIfPresent("THFM", THFM_);
    Info << "  Turbulence heat flux model: " << THFM_ << nl;
    bound(k_, this->kMin_);
    bound(epsilon_, epsilonMin_);
    bound(phit_, phitMin_);
    bound(ebf_, ebfMin_);
    maxMagPb_ = 0;


    if (type == typeName)
    {
        this->printCoeffs(type);
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class BasicMomentumTransportModel>
bool kEpsilonEB<BasicMomentumTransportModel>::read()
{
    if (eddyViscosity<RASModel<BasicMomentumTransportModel>>::read())
    {
        Cmu_.readIfPresent(this->coeffDict());
        C1_.readIfPresent(this->coeffDict());
        C3_.readIfPresent(this->coeffDict());
        C4_.readIfPresent(this->coeffDict());
        C5_.readIfPresent(this->coeffDict());
        C1s_.readIfPresent(this->coeffDict());
        C3s_.readIfPresent(this->coeffDict());
        CK_.readIfPresent(this->coeffDict());
        Ct_.readIfPresent(this->coeffDict());
        Ceps1_.readIfPresent(this->coeffDict());
        Ceps2_.readIfPresent(this->coeffDict());
        CL_.readIfPresent(this->coeffDict());
        Ceta_.readIfPresent(this->coeffDict());
        sigmaK_.readIfPresent(this->coeffDict());
        sigmaEps_.readIfPresent(this->coeffDict());
        sigmaPhit_.readIfPresent(this->coeffDict());
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
void kEpsilonEB<BasicMomentumTransportModel>::correct()
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
    // bound(T_, TMin_);

    const volScalarField L2(typedName("L2"), sqr(Ls()) + L2Min_);

    // const volScalarField::Internal Ceps1
    // (
    //     typedName("Ceps1"),
    //     Ceps1a_*(Ceps1b_ + Ceps1c_*sqrt(1.0/phit_()))
    // );

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

    // Print min/max of k, epsilon, v2 and ebf
    Info << "  k: min/max = "
         << gMin(k_) << " / " << gMax(k_) << nl
          << "  epsilon: min/max = "
          << gMin(epsilon_) << " / " << gMax(epsilon_) << nl
          << "  phi: min/max = "
          << gMin(phit_) << " / " << gMax(phit_) << nl
          << "  ebf: min/max = "
          << gMin(ebf_) << " / " << gMax(ebf_) << nl;
    // Update epsilon (and possibly G) at the wall
    epsilon_.boundaryFieldRef().updateCoeffs();

    // const volScalarField posPb = volScalarField("posPb", pos0(Pb), Pb.dimensions());
    // Info << "posPb dimnesion: " << posPb.dimensions() << endl;
    //
    // Compute strain, vorticity and anisotropy tensors
    tmp<volTensorField> tgradU = fvc::grad(U);

    // Mean strain rate tensor
    const volSymmTensorField S
    (
        symm(tgradU())
    );

    // Wall-normal vectors defined through the elliptic blending factor
    const volVectorField n
    (
        fvc::grad(ebf_)/
        max
        (
            mag(fvc::grad(ebf_)),
            dimensionedScalar(dimless/dimLength, SMALL)
        )
    );

    volVectorField magTermE
    (
        IOobject
        (
            "magTermE",
            this->mesh_.time().name(),
            this->mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mag(2.0*S & n)*n
    );

    const scalar CK_ = 2.3;

    // Additional production term in epsilon eq. (TLLP:Eq.7)
    const volScalarField E
    (
        CK_*pow3(scalar(1.0) - ebf_)*this->nu()*nut*sqr(fvc::div(magTermE))
    );

    // Mean vorticity tensor
    const volTensorField W
    (
        0.5*(tgradU() - tgradU().T())
    );
    tgradU.clear();

    const volSymmTensorField DSDiv_(fvc::ddt(S) + fvc::div(this->phi(), S) - S * fvc::div(U));
                                
    const volTensorField SDS 
    (
        (S & DSDiv_.T())/(2.0*magSqr(S))
    );

    // Spalart-Shur curvature correction for vorticity tensor (TLLP:Eq.20)
    const volTensorField WTilde
    (
        W - 2.0*skew(SDS)
    );
    // Anisotropy tensor (TLLP:Eq.18)
    const dimensionedScalar beta2_ = scalar(2.0)*(scalar(1.0) - C5_)/(C1_ + C1s_ + scalar(1.0));
    volTensorField A
    (
        -scalar(2.0)*nut/k_*(S + scalar(2.0)*beta2_*((S & WTilde) - (WTilde & S))/
        (mag(S + WTilde)))
    );


    // Time scale
    const volScalarField tau
    (
        k_/epsilon_
    );

    // Function for damping phit in region of low strain (TLLP:Eq.17)
    const volScalarField fmu
    (
        (sqrt(2.0)*mag(S)*tau + pow3(ebf_))/
        max
        (
            sqrt(2.0)*mag(S)*tau,
            scalar(1.87)
        )
    );

    // Coefficient using fmu (TLLP:Eq.15)
    const volScalarField Cp3
    (
         fmu/Cmu_*(scalar(2.0/3.0) - C3_/scalar(2.0))
    );

    // Dissipation equation
    tmp<fvScalarMatrix> epsEqn
    (
        fvm::ddt(alpha, rho, epsilon_)
      + fvm::div(alphaRhoPhi, epsilon_)
      - fvm::laplacian(alpha*rho*DepsilonEff(), epsilon_)
      ==
      - fvm::SuSp
        (
            alpha()*rho()*Ceps1_*
            (
                - G()/k_()
                + 2.0/3.0*divU
            ),
            epsilon_
        )
      - fvm::Sp(alpha()*rho()*Ceps2_/tau(), epsilon_)
      + alpha()*rho()*E()
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
        alpha()*rho()*G()
      - fvm::SuSp(2.0/3.0*alpha()*rho()*divU, k_)
      - fvm::Sp(alpha()*rho()/tau(), k_)
      + fvModels.source(alpha, rho, k_)
    );

    kEqn.ref().relax();
    fvConstraints.constrain(kEqn.ref());
    solve(kEqn);
    fvConstraints.constrain(k_);
    bound(k_, this->kMin_);


    // Relaxation function equation
    tmp<fvScalarMatrix> ebfEqn
    (
      -fvm::laplacian(ebf_)
      ==
      -fvm::SuSp(1.0/L2(), ebf_) + 1.0/L2()
    );


    ebfEqn.ref().relax();
    fvConstraints.constrain(ebfEqn.ref());
    solve(ebfEqn);
    fvConstraints.constrain(ebf_);
    bound(ebf_, ebfMin_);

    // Coefficients to be used in the phitEquation
    const dimensionedScalar Cws = Ceps2_ - scalar(1.0) + scalar(5.0) - scalar(1.0)/Cmu_;
    const dimensionedScalar C1Tilde = C1_ + Ceps2_ - scalar(2.0);
    const dimensionedScalar Cp1 = scalar(2.0) - Ceps1_;
    const dimensionedScalar Cp2 = C3s_/sqrt(2.0);
    const dimensionedScalar C4s = scalar(2.0)/Cmu_*(scalar(1.0) - C4_);
    const dimensionedScalar C5s = scalar(2.0)/Cmu_*(scalar(1.0) - C5_);

    // Turbulence stress normal to streamlines equation
    Info << "Solving phitEqn" << endl;
    tmp<fvScalarMatrix> phitEqn
    (
        fvm::ddt(alpha, rho, phit_)
      + fvm::div(alphaRhoPhi, phit_)
      - fvm::laplacian(alpha*rho*DphitEff(), phit_)
      ==
      - fvm::SuSp
        (
            alpha()*rho()*
            (
                (1.0 - pow3(ebf_()))*Cws/tau()
                + pow3(ebf_())*(C1Tilde + C1s_*(G() - (2.0/3.0)*k_()*divU)
                /epsilon_())/tau()
                + Cp1*(G()/k_() - (2.0/3.0)*divU)
                - pow3(ebf_())*Cp2*sqrt(2.0)*mag(S())
            )
          , phit_
        )
        + alpha()*rho()*
        (
            pow3(ebf_())/tau()/(2.0*magSqr(S()))*((C4s*(A() & S())
                - C5s*(A() & WTilde())) && S())
            + pow3(ebf_())*Cp3()/tau()
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
