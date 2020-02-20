/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2017 OpenCFD Ltd.
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

#include "phaseChangeReaction.H"
#include "constants.H"
#include "fvcGrad.H"
#include "fvcSnGrad.H"
#include "fvcDiv.H"
#include "surfaceInterpolate.H"
#include "fvcReconstruct.H"
#include "fvm.H"
#include "zeroGradientFvPatchFields.H"

using namespace Foam::constant;

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //




// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::phaseChangeReaction::phaseChangeReaction
(
    const dictionary& dict,
    const fvMesh& mesh,
    const volScalarField& C,
    const volScalarField& alpha
)
:
    Cu_(dict.get<scalar>("Cu")),
    K_("K", dimVelocity, dict),
    Cactivate_("Cactivate", dimMoles/dimVolume, dict),
    Mv_("Mv", dimMass/dimMoles, dict),
    alphaMax_(dict.lookupOrDefault<scalar>("alphaMax", 1.0)),
    alphaMin_(dict.lookupOrDefault<scalar>("alphaMin", 0.5)),
    alphaSolidMin_(dict.lookupOrDefault<scalar>("alphaSolidMin", 0.5)),
    alphaRestMax_(dict.lookupOrDefault<scalar>("alphaRestMax", 0.01)),
    smoothSurface_(dict.lookupOrDefault<bool>("smoothSurface", false)),
    smoothAreaDensity_(dict.lookupOrDefault<scalar>("smoothAreaDensity", 1.0)),
    alpha_(alpha),
    C_(C),
    mesh_(mesh)
{
    // Convert from g/mol to Kg/mol
    Mv_.value() = Mv_.value() * 1e-3;
    Info<< "Constructor of phaseChangeReaction finished" << endl;
}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField> Foam::phaseChangeReaction::Kexp(const volScalarField& field)
{
    const volScalarField& from = alpha_;

    const volScalarField to
    (
        "to",
        1.0-alpha_
    );

    Info<< "alpha from max/min: "<< max(from).value() << ", " << min(from).value() << endl;
    Info<< "alpha to max/min: "<< max(to).value() << ", " << min(to).value() << endl;
    const volVectorField gradFrom(fvc::grad(from));
    const volVectorField gradTo(fvc::grad(to));
    const volScalarField areaDensitySmooth
    (
        IOobject
        (
            "areaDensitySmooth",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("areaDensitySmoothInput", dimless/dimLength, smoothAreaDensity_)
    );

    const volScalarField areaDensityGrad
    (
        "areaDensityGrad",
        2*mag(gradFrom)
    );
    //areaDensityGrad = areaDensityGradTmp;
    Info<< "areaDensityGrad min/max : "<<min(areaDensityGrad).value()<<", "<< max(areaDensityGrad).value()<<endl;
    
    const volScalarField gradAlphaf(gradFrom & gradTo);

    volScalarField Cmask("Cmask", from*0.0);
    Cmask.correctBoundaryConditions();
    Info<< "calculate Cmask..." <<endl;

    // mark cells that is next to the solid phase cell
    // within same processor
    forAll(Cmask, celli)
    {
        if (gradAlphaf[celli] < 0)
        {
            if (from[celli] > alphaMin_ && from[celli] < alphaMax_)
            {
                {
                    scalar alphaRes = 1.0 - from[celli] - to[celli];
                    if (alphaRes < alphaRestMax_)
                    {
                        Cmask[celli] = 1.0;
                    }
                }
            }
        }
        //- check nearby cell within same processor
        bool flag = false;
        forAll(mesh_.cellCells()[celli],cellj)
        {
    
            if (to[mesh_.cellCells()[celli][cellj]] > alphaSolidMin_)
            {
                flag = true;
            }
        }
        if (flag)
        {
            Cmask[celli] = 1.0;
        }
        else
        {
            Cmask[celli] = 0.0;
        }
    }

    //- check mesh boundaries of the processor
    forAll(mesh_.boundaryMesh(), patchI)
    {
        const fvPatchScalarField& pf = to.boundaryField()[patchI];
        const labelList& faceCells = pf.patch().faceCells();

        //- Coupled boundaries (processor, cylic, etc)
        if(pf.coupled())
        {
            scalarField neighbors = pf.patchNeighbourField();

            forAll(faceCells, faceI)
            {
                if (neighbors[faceI] > alphaSolidMin_)
                {
                    Cmask[faceCells[faceI]] = 1.0;
                }
            }
        }
    } 

    tmp<volScalarField> tRhom
    (
        new volScalarField
        (
            IOobject
            (
                "trhom",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh_,
            dimensionedScalar(dimDensity, Zero)
        )
    );
    volScalarField& rhom = tRhom.ref();

    tmp<volScalarField> tTdelta
    (
        new volScalarField
        (
            IOobject
            (
                "tTdelta",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh_,
            dimensionedScalar(dimless, Zero)
        )
    );
    volScalarField& tDelta = tTdelta.ref();
    Info<< "calculate rhom&tDelta..."<<endl;

    if (sign(K_.value()) > 0)
    {
        //- convert mol/m3 to kg/m3 for OpenFOAM
        rhom =
            (C_-Cactivate_)*Mv_;

        tDelta = max
        (
            pos(C_*Cmask - Cactivate_),//dimensionedScalar("C1", dimMoles/dimVolume, 1.0),
            dimensionedScalar("C0", dimless, Zero)
        );
//            Info<<"tDelta: " << tDelta << endl;
    }
    else
    {
        Info<< "Unit K_ is negative or zero!!!" << endl;
    }
    Info<< "calculate massFluxPrec..."<<endl;


    volScalarField massFluxPrec
    (
        "massFluxPrec",
        K_
        * rhom
        * tDelta
        * pos(from-0.01)
    );

    if (mesh_.time().outputTime())
    {
        areaDensityGrad.write();
        Cmask.write();
        to.write();
//            volScalarField mKGasDot
//            (
//                "mKGasDot",
//                massFluxPrec*areaDensity*Nl*from
//            );
//            mKGasDot.write();
    }


    //- Constrain the reaction rate to prevent negative concentration
    if(smoothSurface_)
    {
        forAll(mesh_.C(), cellI)
        {
            if(massFluxPrec[cellI]*areaDensitySmooth[cellI]*mesh_.time().deltaTValue()
                > mag(C_[cellI]-Cactivate_.value())*Mv_.value())
                {
                    massFluxPrec[cellI] = mag(C_[cellI]-Cactivate_.value())*Mv_.value()
                                        / (mesh_.time().deltaTValue()*areaDensitySmooth[cellI]+VSMALL);
                }
        }   
    }
    else
    {
        forAll(mesh_.C(), cellI)
        {
            if(massFluxPrec[cellI]*areaDensityGrad[cellI]*mesh_.time().deltaTValue()
                > mag(C_[cellI]-Cactivate_.value())*Mv_.value())
                {
                    massFluxPrec[cellI] = mag(C_[cellI]-Cactivate_.value())*Mv_.value()
                                        / (mesh_.time().deltaTValue()*areaDensityGrad[cellI]+VSMALL);
                }
        }
    }

    Info<< "precipitate return: " << min(massFluxPrec*areaDensitySmooth).value() << ", "<< max(massFluxPrec*areaDensitySmooth).value() << endl;

    if(smoothSurface_)
    {
        Info<< "Smooth surface areaDensity method used!" << endl;
        return massFluxPrec * areaDensitySmooth;
    }
    else
    {
        Info<< "Gradient surface areaDensity method used!" << endl;
        return massFluxPrec * areaDensityGrad;            
    }
}

Foam::tmp<Foam::fvScalarMatrix> 
Foam::phaseChangeReaction::preRect
(
    const volScalarField& C,
    volScalarField& SuOutPut
)
{
    Info<< "Start preRect calculation" << endl;
    
    tmp<fvScalarMatrix> cEqnPtr
    (
        new fvScalarMatrix(C, dimMoles/dimTime)
    );

    fvScalarMatrix& eqn = cEqnPtr.ref();
    // Net mass transfer from k to i phase
    tmp<volScalarField> cdmdtNetki
    (
        new volScalarField
        (
            IOobject
            (
                "cdmdtYki",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh_,
            dimensionedScalar(dimDensity/dimTime, Zero)
        )
    );
    volScalarField& dmdtNetki = cdmdtNetki.ref();

    // Explicit mass transfer rate
    Info<< "Calculate KExp" << endl;
    tmp<volScalarField> KExp = Kexp(C);
    Info<< "End of calculation of KExp" << endl;

    if (KExp.valid())
    {
        Info<< "Fluid to solid calculated!" << endl;
        dmdtNetki += KExp.ref();
        Info<< "Kexp output min/max: " << min(KExp.ref()) << ", " << max(KExp.ref()) << endl;
    }

    SuOutPut = dmdtNetki;
    eqn -= (dmdtNetki/Mv_);

    return cEqnPtr;
}

void Foam::phaseChangeReaction::addInterfacePorosity(fvVectorMatrix& UEqn)
{
    const scalarField& Vc = mesh_.V();
    scalarField& Udiag = UEqn.diag();

    //- Voller Prakash interfacial porosity model

    const volScalarField& liquidAlpha = alpha_;

    tmp<volScalarField> STerm 
    (
        new volScalarField
        (
            IOobject
            (
                "STerm",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh_,
            dimensionedScalar(dimless, Zero)
        )
    );
    volScalarField& STermRef = STerm.ref();

    STermRef = Cu_*sqr(1.0-liquidAlpha)/(pow3(liquidAlpha) + 1e-3);

    Udiag += Vc*STermRef;
}

Foam::phaseChangeReaction::~phaseChangeReaction(){};

// ************************************************************************* //
