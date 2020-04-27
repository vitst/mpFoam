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
#include <random>

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
    pKsp_(dict.get<scalar>("pKsp")),
    lnA_(dict.get<scalar>("lnA")),
    gamma_(dict.get<scalar>("gamma")),
    Vmol_(dict.get<scalar>("Vmol")),
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
    debug_(dict.lookupOrDefault<bool>("debug", false)),
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

    dimensionedScalar unitConv
    (
        "unitConv",
        dimensionSet(0,3,0,0,-1,0,0),
        scalar(1.0)
    );

    if (sign(K_.value()) > 0)
    {
        //- convert mol/m3 to kg/m3 for OpenFOAM
        rhom =
            (C_-Cactivate_)*(C_-Cactivate_)*Mv_*unitConv;

        tDelta = max
        (
            pos(C_*Cmask - Cactivate_),//dimensionedScalar("C1", dimMoles/dimVolume, 1.0),
            dimensionedScalar("C0", dimless, Zero)
        );
            Info<<"rhom min/max: " << min(rhom).value() << ", "<< max(rhom).value() << endl;
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
//    if(smoothSurface_)
//    {
//        forAll(mesh_.C(), cellI)
//        {
//            if(massFluxPrec[cellI]*areaDensitySmooth[cellI]*mesh_.time().deltaTValue()
//                > mag(C_[cellI]-Cactivate_.value())*Mv_.value())
//                {
//                    massFluxPrec[cellI] = mag(C_[cellI]-Cactivate_.value())*Mv_.value()
//                                        / (mesh_.time().deltaTValue()*areaDensitySmooth[cellI]+VSMALL);
//                }
//        }   
//    }
//    else
//    {
//        forAll(mesh_.C(), cellI)
//        {
//            if(massFluxPrec[cellI]*areaDensityGrad[cellI]*mesh_.time().deltaTValue()
//                > mag(C_[cellI]-Cactivate_.value())*Mv_.value())
//                {
//                    massFluxPrec[cellI] = mag(C_[cellI]-Cactivate_.value())*Mv_.value()
//                                        / (mesh_.time().deltaTValue()*areaDensityGrad[cellI]+VSMALL);
//                }
//        }
//    }

    Info<< "precipitate return: " << min(massFluxPrec).value() << ", "<< max(massFluxPrec).value() << endl;

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

void Foam::phaseChangeReaction::nuSiteCal
(
    volScalarField& nuSite,
    volScalarField& cryDomain,
    vectorList& nuSiteList
)
{
    //Calculate the possibilities for nucleation
    tmp<volScalarField> nuRate = nuRateCal();
    Info<< "Max/Min nucleation rate(n/m2/s): " << max(nuRate.ref()).value() << ", " << min(nuRate.ref()).value()<< endl;

    //Loop through wall boundary patches to estimate nucleation
    const fvPatchList& patches = mesh_.boundary();
    labelList wallList;
    wallList.clear();

    forAll(patches, patchI)
    {
        const fvPatch& p = patches[patchI];
        if(p.type()=="wall")
        {
            forAll(p, pFaceI)
            {
                label faceCelli = p.faceCells()[pFaceI];
                wallList.resize(wallList.size()+1);
                wallList[wallList.size()-1]=faceCelli;
            }
        }
    }

    std::random_device rd;  //Will be used to obtain a seed for the random number engine
    std::mt19937 gen(rd());

    forAll(wallList, i)
    {
        label poss = rand() % 1000000000000 + 1;
        const cell& faces = mesh_.cells()[wallList[i]];
        scalar faceArea = 0.0;
        forAll(faces, faceI)
        {
            if(mesh_.magSf()[faceI]>faceArea)
            {
                faceArea = mesh_.magSf()[faceI];
            }
        }

        std::uniform_real_distribution<> dis(0.0, 1.0);

        if(debug_)
        {
            Info<< "Max faceArea: " << faceArea << endl;
            Info<< "Random integer generated: " << poss << endl;
            Info<< "Nucleation rate: " << nuRate.ref()[wallList[i]] << endl;
            Info<< "Computational time: " << mesh_.time().value() << endl;
        }

        if((nuRate.ref()[wallList[i]]/1e12*1e2)>(dis(gen)*faceArea*1e12))
        {
            nuSite[wallList[i]] = 1.0;
        }
        // mark cells for check, remove in running code
        // nuSite[wallList[i]] = 1.0;
    }

    cryCons(nuSite, cryDomain, nuSiteList);
}

Foam::tmp<Foam::volScalarField> 
Foam::phaseChangeReaction::CtoSI
()
{
    tmp<volScalarField> SItmp 
    (
        new volScalarField
        (
            IOobject
            (
                "SItmp",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh_,
            dimensionedScalar(dimless, Zero)
        )
    );
    volScalarField& SItmpRef = SItmp.ref();

    scalar Ksp = pow(10,-pKsp_);
    Info<< "KSP value: " << Ksp << endl;

    // mol/m3 to SI
    forAll(mesh_.C(), cellI)
    {
        SItmpRef[cellI] = max(0.001, log10((pow(C_[cellI],2.0)+VSMALL)/Ksp));
    };

    return SItmp;
}

Foam::tmp<Foam::volScalarField> 
Foam::phaseChangeReaction::nuRateCal
()
{
    tmp<volScalarField> nuRate 
    (
        new volScalarField
        (
            IOobject
            (
                "nuRate",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh_,
            dimensionedScalar(dimArea/dimTime, Zero)
        )
    );
    volScalarField& nuRateRef = nuRate.ref();

    // obtain SI field from concentration
    tmp<volScalarField> SI = CtoSI();    
    Info<< "SI output min/max: " << min(SI.ref()) << ", " << max(SI.ref()) << endl;

    // define constants
    scalar kb = 1.38064852e-23;
    scalar NA = 6.02214076e23;

    forAll(mesh_.C(),cellI)
    {
        nuRateRef[cellI] = exp(lnA_)*exp((-16.0*3.14159*pow(gamma_,3.0)*pow(Vmol_/NA,2.0))/(3.0*2.303*2.303*pow(SI.ref()[cellI],2.0)*pow(kb*293,3.0)));
    }

    return nuRate;
}

Foam::vectorList  
Foam::phaseChangeReaction::extractNuSite
(
    const volScalarField& nuSite,
    vectorList& nuSiteList
)
{

    // store the nucleation site coordinate
    forAll(mesh_.C(), cellI)
    {
        if(nuSite[cellI]==1.0)
        {
            nuSiteList.resize(nuSiteList.size()+1);
            nuSiteList[nuSiteList.size()-1]=mesh_.C()[cellI];
            //Info<< "Nucleation site coordinates: " << nuSiteList[nuSiteList.size()-1] << endl;
        }
    }

    return nuSiteList;
}

void Foam::phaseChangeReaction::cryCons
(
    const volScalarField& nuSite,
    volScalarField& cryDomain,
    vectorList& nuSiteList
)
{
    tmp<volScalarField> tmpCryDomain 
    (
        new volScalarField
        (
            IOobject
            (
                "tmpCryDomain",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh_,
            dimensionedScalar(dimless, Zero)
        )
    );
    volScalarField& tmpCryDomainRef = tmpCryDomain.ref();

    vectorList nuSiteListTmp = extractNuSite(nuSite, nuSiteList);
    Info<< "Total amount of nucleation sites: " << nuSiteListTmp.size() << endl;

    // calculate the volume occupied by the constrained crystal domain from cell center
    forAll(mesh_.C(), cellI)
    {
        scalar xCell = mesh_.C()[cellI].component(0);
        scalar yCell = mesh_.C()[cellI].component(1);
        scalar zCell = mesh_.C()[cellI].component(2);
        for(int nuI=0; nuI<nuSiteListTmp.size(); nuI++)
        {
            if(cellI==0)
            {
                Info<< "operation #: " << nuI << endl;
            }
            scalar xCoorMin = nuSiteListTmp[nuI][0]-4.1e-6;
            scalar xCoorMax = nuSiteListTmp[nuI][0]+4.1e-6;
            scalar yCoorMin = nuSiteListTmp[nuI][1]-4.1e-6;
            scalar yCoorMax = nuSiteListTmp[nuI][1]+4.1e-6;
            scalar zCoorMin = nuSiteListTmp[nuI][2]-2.1e-6;
            scalar zCoorMax = nuSiteListTmp[nuI][2]+2.1e-6;
            if((((((xCoorMin<=xCell) && (xCell<=xCoorMax)) && (yCoorMin<=yCell)) && (yCell<=yCoorMax)) && (zCoorMin<=zCell)) && (zCell<=zCoorMax))
            {
                tmpCryDomainRef[cellI] = 1.0;
            }
        }
    }

    forAll(mesh_.C(), cellI)
    {
        cryDomain[cellI] = tmpCryDomainRef[cellI];
    }
}

Foam::phaseChangeReaction::~phaseChangeReaction(){};

// ************************************************************************* //
