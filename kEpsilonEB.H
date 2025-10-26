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

Class
    Foam::RASModels::kEpsilonEB

Description
    Lien and Kalitzin's v2-f turbulence model for incompressible and
    compressible flows, with a limit imposed on the turbulent viscosity given
    by Davidson et al.

    The model solves for turbulence kinetic energy k and turbulence dissipation
    rate epsilon, with additional equations for the turbulence stress normal to
    streamlines, v2, and elliptic damping function, f.

    The variant implemented employs N=6, such that f=0 on walls.

    Wall boundary conditions are:

        k       = kLowReWallFunction
        epsilon = epsilonWallFunction
        phit      = v2WallFunction
        ebf       = fixedValue 0;

    These are applicable to both low- and high-Reynolds number flows.

    Inlet values can be approximated by:

        phit      = 2/3 k
        ebf       = zero-gradient

    References:
    \verbatim
        Lien, F. S., & Kalitzin, G. (2001).
        Computations of transonic flow with the v2f turbulence model.
        International Journal of Heat and Fluid Flow, 22(1), 53-61.

        Davidson, L., Nielsen, P., & Sveningsson, A. (2003).
        Modifications of the v2-f model for computing the flow in a
        3D wall jet.
        Turbulence, Heat and Mass Transfer, 4, 577-584
    \endverbatim

    The default model coefficients are
    \verbatim
        kEpsilonEBCoeffs
        {
            Cmu         0.22;
            CmuKEps     0.09;
            C1          1.4;
            C2          0.3;
            CL          0.23;
            Ceta        70;
            Ceps2       1.9;
            Ceps3       -0.33;
            sigmaEps    1.3;
            sigmaK      1;
        }
    \endverbatim

    Note:
        If the kLowReWallFunction is employed, a velocity variant of the
        turbulent viscosity wall function should be used, e.g. nutUWallFunction.
        Turbulence k variants (nutk...) for this case will not behave correctly.

See also
    Foam::RASModels::v2fBase
    Foam::RASModels::kEpsilon
    Foam::kLowReWallFunctionFvPatchScalarField
    Foam::epsilonWallFunctionFvPatchScalarField
    Foam::v2WallFunctionFvPatchScalarField
    Foam::fWallFunctionFvPatchScalarField

SourceFiles
    kEpsilonEB.C

\*---------------------------------------------------------------------------*/

#ifndef kEpsilonEB_H
#define kEpsilonEB_H

// #include "v2fBase.H"
#include "RASModel.H"
#include "eddyViscosity.H"
#include "uniformDimensionedFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace RASModels
{

/*---------------------------------------------------------------------------*\
                             Class kEpsilonEB Declaration
\*---------------------------------------------------------------------------*/

template<class BasicMomentumTransportModel>
class kEpsilonEB
:
    public eddyViscosity<RASModel<BasicMomentumTransportModel>>
    // public v2fBase
{

protected:

    // Protected data

        Switch includeNu_;

        // Model coefficients

            dimensionedScalar Cmu_;
            dimensionedScalar C1_;
            dimensionedScalar C3_;
            dimensionedScalar C4_;
            dimensionedScalar C5_;
            dimensionedScalar C1s_;
            dimensionedScalar C3s_;
            dimensionedScalar CK_;
            dimensionedScalar Ct_;
            dimensionedScalar Ceps1_;
            dimensionedScalar Ceps2_;
            dimensionedScalar CL_;
            dimensionedScalar Ceta_;
            dimensionedScalar sigmaK_;
            dimensionedScalar sigmaEps_;
            dimensionedScalar sigmaPhit_;

          scalar maxMagPb_;

        //- Gravitational acceleration vector [m/s^2]
          const uniformDimensionedVectorField& g_;

        //- Inverse turbulent Prandtl number
          dimensionedScalar Cg_; 

        //- Coefficient for the buoyancy production term in GGDH
          dimensionedScalar Cphi_;

        //- Turbulent heat flux model (THFM) type.
          //  Eiter SGDH (simple gradient diffusion hypothesis, default) or
          //  or GGDH (generalized gradient diffusion hypothesis).
          word THFM_;

        // Fields

            //- Turbulence kinetic energy
            volScalarField k_;

            //- Turbulence dissipation
            volScalarField epsilon_;

            //- Turbulence stress normal to streamlines
            volScalarField phit_;

            //- Damping function
            volScalarField ebf_;

            volScalarField T_;

        // Bounding values

            dimensionedScalar phitMin_;
            dimensionedScalar ebfMin_;
            dimensionedScalar epsilonMin_;
            dimensionedScalar L2Min_;


    // Protected Member Functions
  
                //- Return the effective diffusivity for phit (LUU:Eq. 17)
        tmp<volScalarField> DphitEff() const
        {
            return tmp<volScalarField>
            (
                new volScalarField
                (
                    "DphitEff",
                    this->nut_/sigmaPhit_ + this->nu()/2
                )
            );
        }

        //- Correct the eddy-viscosity nut
        virtual void correctNut();

        //- Return time scale, Ts
        tmp<volScalarField> Ts() const;

        //- Return length scale, Ls
        tmp<volScalarField> Ls() const;


public:

    typedef typename BasicMomentumTransportModel::alphaField alphaField;
    typedef typename BasicMomentumTransportModel::rhoField rhoField;


    //- Runtime type information
    TypeName("kEpsilonEB");


    // Constructors

        //- Construct from components
        kEpsilonEB
        (
            const alphaField& alpha,
            const rhoField& rho,
            const volVectorField& U,
            const surfaceScalarField& alphaRhoPhi,
            const surfaceScalarField& phi,
            const viscosity& viscosity,
            const word& type = typeName
        );


    //- Destructor
    virtual ~kEpsilonEB()
    {}


    // Member Functions

        //- Read RASProperties dictionary
        virtual bool read();

        //- Return the effective diffusivity for k
        tmp<volScalarField> DkEff() const
        {
            return volScalarField::New
            (
                "DkEff",
                this->nut_/sigmaK_ + this->nu()
            );
        }

        //- Return the effective diffusivity for epsilon
        tmp<volScalarField> DepsilonEff() const
        {
            return volScalarField::New
            (
                "DepsilonEff",
                this->nut_/sigmaEps_ + this->nu()
            );
        }

        //- Return the turbulence kinetic energy
        virtual tmp<volScalarField> k() const
        {
            return k_;
        }

        //- Return the turbulence kinetic energy dissipation rate
        virtual tmp<volScalarField> epsilon() const
        {
            return epsilon_;
        }

        //- Return the turbulence specific dissipation rate
        virtual tmp<volScalarField> omega() const
        {
            return volScalarField::New
            (
                "omega",
                epsilon_/(Cmu_*k_)
            );
        }

        //- Return turbulence stress normal to streamlines
        virtual tmp<volScalarField> v2() const
        {
            return phit_;
        }

        //- Return turbulence stress normal to streamlines
        virtual tmp<volScalarField> phit() const
        {
            return phit_;
        }

        //- Return the damping function
        virtual tmp<volScalarField> ebf() const
        {
            return ebf_;
        }

        //- Solve the turbulence equations and correct the turbulence viscosity
        virtual void correct();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "kEpsilonEB.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
