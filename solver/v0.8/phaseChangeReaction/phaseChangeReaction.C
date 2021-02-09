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
    const volScalarField& alpha,
    volScalarField& Cmask
)
:
    dirGrowthFlag_(dict.lookupOrDefault<bool>("dirGrowth", false)),
    dirGrowthNum_(dict.lookupOrDefault<scalar>("dirGrowthNum", 1.0)),
    dirGrowthNorm_(dict.lookupOrDefault<List<vector>>("dirGrowthNorm", {vector::zero})),
    dirFactor_(dict.lookupOrDefault<List<scalar>>("dirFactor", {1.0})),
    Cu_(dict.lookupOrDefault<scalar>("Cu", 1e15)),
    pKsp_(dict.lookupOrDefault<scalar>("pKsp", 1.0)),
    lnA_(dict.lookupOrDefault<scalar>("lnA", 1.0)),
    gamma_(dict.lookupOrDefault<scalar>("gamma", 1.0)),
    Vmol_(dict.lookupOrDefault<scalar>("Vmol", 1.0)),
    reacModel_(dict.lookupOrDefault<word>("reactionModel", "linear")),
    K_("K", dimMoles/dimArea/dimTime, dict),
    KList_(dict.lookupOrDefault<List<scalar>>("KList", {1.0})),
    Cactivate_("Cactivate", dimMoles/dimVolume, dict),
    Mv_("Mv", dimMass/dimMoles, dict),
    alphaMax_(dict.lookupOrDefault<scalar>("alphaMax", 1.1)),
    alphaMin_(dict.lookupOrDefault<scalar>("alphaMin", 0.01)),
    alphaSolidMin_(dict.lookupOrDefault<scalar>("alphaSolidMin", 0.5)),
    alphaRestMax_(dict.lookupOrDefault<scalar>("alphaRestMax", 0.01)),
    smoothSurface_(dict.lookupOrDefault<bool>("smoothSurface", false)),
    smoothAreaDensity_(dict.lookupOrDefault<scalar>("smoothAreaDensity", 1.0)),
    cornerCell_(dict.lookupOrDefault<bool>("cornerCell", false)),
    gradLim_(dict.lookupOrDefault<scalar>("gradientLimit", 1e20)),
    alpha_(alpha),
    Cmask_(Cmask),
    debug_(dict.lookupOrDefault<bool>("debug", false)),
    nucleationFlag_(dict.lookupOrDefault<bool>("nucleation", false)),
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
    volScalarField from
    (
        "from",
        alpha_
    );
    from.correctBoundaryConditions();

    volScalarField to
    (
        "to",
        1.0-alpha_
    );
    to.correctBoundaryConditions();

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
    Info<< "areaDensityGrad min/max : "<<min(areaDensityGrad).value()<<", "<< max(areaDensityGrad).value()<<endl;
    
    const volScalarField gradAlphaf(gradFrom & gradTo);

    dimensionedScalar VSUnit
    (
        "VSUnit",
        dimensionSet(0,-1,0,0,0,0,0),
        VSMALL
    );

    //- construct directional gradient to track interface norm
    if (debug_)
    {
        const vector xDir(1, 0, 0);
        const vector yDir(0, 1, 0);
        const vector zDir(0, 0, 1);
        const vector dir210(0.447, 0.894, 0);
        const vector dir210Neg(0.447, -0.894, 0);

        const volScalarField xDirGradSolid
        (
            "xDirGradSolid",
            (gradTo & xDir) 
        );

        const volScalarField yDirGradSolid
        (
            "yDirGradSolid",
            (gradTo & yDir) 
        );

        const volScalarField zDirGradSolid
        (
            "zDirGradSolid",
            (gradTo & zDir) 
        );

        const volScalarField dir210GradSolid
        (
            "dir210GradSolid",
            (gradTo & dir210) 
        );

        const volScalarField dir210NegGradSolid
        (
            "dir210NegGradSolid",
            (gradTo & dir210Neg) 
        );

        if (mesh_.time().outputTime())
        {
            xDirGradSolid.write();
            yDirGradSolid.write();
            zDirGradSolid.write();
            dir210GradSolid.write();
            dir210NegGradSolid.write();
        }
    }
    
    //- Update the reaction constant field basing on surface norm
    tmp<volScalarField> tkConst
    (
        new volScalarField
        (
            IOobject
            (
                "tkConst",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh_,
            dimensionedScalar(dimMoles/dimArea/dimTime, Zero)
        )
    );
    volScalarField& kConst = tkConst.ref();

    if(dirGrowthFlag_)
    {
        forAll(kConst, cellI)
        {
            //- Calculate similarity of surface norm and directional vectors
            scalar tmpDir = 0.0;
            for (size_t i = 0; i < dirGrowthNum_; i++)
            {
                scalar tCosTmp = ((gradTo[cellI] & dirGrowthNorm_[i])/(mag(gradTo[cellI])+SMALL))*dirFactor_[i];
                if (tmpDir<mag(tCosTmp))
                {
                    kConst[cellI] = KList_[i]*Cmask_[cellI];
                    tmpDir = mag(tCosTmp);
                }
            }
        }
    }
    else
    {
        forAll(kConst, cellI)
        {
            kConst[cellI] = K_.value();
        }
    }

    //- Define units for kinetic equation to suit input format
    FixedList<scalar,7> dims;
    if(reacModel_ == "linear")
    {
        dims = {1, -3, 0, 0, 0, 0, 0};
    }
    else if(reacModel_ == "SI")
    {
        dims = {1, 0, 0, 0, -1, 0, 0};
    }

    dimensionSet unit(dims);

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
            dimensionedScalar(unit, Zero)
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
        if(reacModel_=="linear")
        {
            //- convert mol/m3 to kg/m3 for OpenFOAM
            rhom =
                (C_-Cactivate_)*Mv_;

            tDelta = max
            (
                pos(C_*Cmask_ - Cactivate_),//dimensionedScalar("C1", dimMoles/dimVolume, 1.0),
                dimensionedScalar("C0", dimless, Zero)
            );
            Info<<"rhom min/max: " << min(rhom).value() << ", "<< max(rhom).value() << endl; 
        }
        else if(reacModel_=="SI")
        {
            //- convert mol/m3 to kg/m3 for OpenFOAM
            rhom =
                (pow(C_,2.0)/pow(Cactivate_,2.0)-1.0)*Mv_;

            tDelta = max
            (
                pos(C_*Cmask_ - Cactivate_),//dimensionedScalar("C1", dimMoles/dimVolume, 1.0),
                dimensionedScalar("C0", dimless, Zero)
            );
            Info<<"rhom min/max: " << min(rhom).value() << ", "<< max(rhom).value() << endl; 
        }
    }
    else
    {
        Info<< "Unit K_ is negative or zero!!!" << endl;
    }
    Info<< "calculate massFluxPrec..."<<endl;

    volScalarField massFluxPrec
    (
        "massFluxPrec",
        kConst//K_
        * rhom
        * tDelta
        * pos(from-0.01)
    );


    if (mesh_.time().outputTime())
    {
        areaDensityGrad.write();
        //Cmask.write();
        to.write();
        kConst.write();
//            volScalarField mKGasDot
//            (
//                "mKGasDot",
//                massFluxPrec*areaDensity*Nl*from
//            );
//            mKGasDot.write();
    }

    dimensionedScalar totReactionRate_(0.0);

    if(reacModel_=="linear")
    {
        totReactionRate_ = gSum((mesh_.V()*massFluxPrec*areaDensityGrad*unitConv)());
    }
    else if(reacModel_=="SI")
    {
        totReactionRate_ = gSum((mesh_.V()*massFluxPrec*areaDensityGrad)());
    }


    Info<< "precipitate return: " << min(massFluxPrec).value() << ", "<< max(massFluxPrec).value() << endl;
    Info<< "Total reaction rate: " << totReactionRate_.value() << endl;

    if(smoothSurface_)
    {
        Info<< "Smooth surface areaDensity method used!" << endl;
        return massFluxPrec * areaDensitySmooth;
    }
    else
    {
        Info<< "Gradient surface areaDensity method used!" << endl;
        if(reacModel_=="linear")
        {
            return massFluxPrec * areaDensityGrad * unitConv;
        }
        else if(reacModel_=="SI")
        {
            return massFluxPrec * areaDensityGrad;
        }
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

    STermRef = Cu_*sqr(1.0-liquidAlpha)/(pow3(liquidAlpha) + 1e-3)*Cmask_;

    Udiag += Vc*STermRef;
}

void Foam::phaseChangeReaction::updateCmask()
{
    //- Readjust in next version for better coding
    const volScalarField& from = alpha_;

    volScalarField to
    (
        "to",
        1.0-alpha_
    );
    to.correctBoundaryConditions();

    const volVectorField gradFrom(fvc::grad(from));
    const volVectorField gradTo(fvc::grad(to));
    const volScalarField gradAlphaf(gradFrom & gradTo);

    Info<< "calculate Cmask..." <<endl;

    // mark cells that is next to the solid phase cell
    // within same processor
    forAll(Cmask_, celli)
    {
        if (gradAlphaf[celli] < 0)
        {
            if (from[celli] > alphaMin_ && from[celli] < alphaMax_)
            {
                if(mag(gradFrom[celli])>gradLim_)
                {
                    Cmask_[celli] = 1.0;
                }
            }
        }
        //- check nearby cell within same processor
        bool flag = false;
        labelListList cellCellsHolder;
        cellCellsHolder.clear();

        forAll(mesh_.cellCells()[celli],cellj)
        {
            if (to[celli] > alphaSolidMin_)
            {
                Cmask_[mesh_.cellCells()[celli][cellj]] = 1.0;

                //- store cellCells lists of adjacent cells
                label dirNearbyCellMarker = mesh_.cellCells()[celli][cellj];
                cellCellsHolder.resize(cellCellsHolder.size()+1);
                labelList dirNearbyCellCells = mesh_.cellCells()[dirNearbyCellMarker];
                cellCellsHolder.append(dirNearbyCellCells);
                //Info<< dirNearbyCellCells << endl;

                flag = true;
            }
        }

        //- enable corner cell growth for cartisian grid
        if(flag == true && cornerCell_ == true)
        {
            labelList cornerCellList;
            cornerCellList.clear();
            forAll(cellCellsHolder, iterI)
            {
                forAll(cellCellsHolder, iterJ)
                {
                    if(iterI != iterJ)
                    {
                        forAll(cellCellsHolder[iterI], labelI)
                        {
                            if(cellCellsHolder[iterJ].found(cellCellsHolder[iterI][labelI]))
                            {
                                if(!cornerCellList.found(cellCellsHolder[iterI][labelI]))
                                {
                                    cornerCellList.append(cellCellsHolder[iterI][labelI]);
                                }
                            }                            
                        }
                    }
                }
            }
            //Info<< "corner cell list: " << cornerCellList << endl;


            forAll(cornerCellList, cellI)
            {
                Cmask_[cornerCellList[cellI]] = 1.0;
            }
        }

        //- nucleation site enable growth
        // if(nucleationFlag_)
        // {
        //     if(nuSite[celli]==1.0)
        //     {
        //         Cmask_[celli] = 1.0;
        //     }
        // }
        //if (flag)
        //{
        //    Cmask_[celli] = 1.0;
        //}
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
                //if (neighbors[faceI]>alphaSolidMin_)
                //{
                //    Info<< faceI << endl;
                //    Info<< neighbors[faceI] << endl;
                //}
                
                if (neighbors[faceI] >= alphaSolidMin_)
                {
                    Cmask_[faceCells[faceI]] = 1.0;
                }
            }
        }
    } 

    Cmask_.correctBoundaryConditions();
}

void Foam::phaseChangeReaction::nuSiteCal
(
    volScalarField& nuSite,
    volScalarField& nuRateOut,
    vectorList& nuSiteList,
    dimensionedScalar& nuTotal_
)
{
    //Create tmp nucleation field
    tmp<volScalarField> nuSitePerStep
    (
        new volScalarField
        (
            IOobject
            (
                "nuSitePerStep",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh_,
            dimensionedScalar(dimless, Zero)
        )
    );
    volScalarField& nuSitePerStepRef = nuSitePerStep.ref();

    //tmp<volScalarField> wallMarker
    //(
    //    new volScalarField
    //    (
    //        IOobject
    //        (
    //            "wallMarker",
    //            mesh_.time().timeName(),
    //            mesh_,
    //            IOobject::NO_READ,
    //            IOobject::AUTO_WRITE
    //        ),
    //        mesh_,
    //        dimensionedScalar(dimless, Zero)
    //    )
    //);
    //volScalarField& wallMarkerRef = wallMarker.ref();

    tmp<volScalarField> faceAreaTmp
    (
        new volScalarField
        (
            IOobject
            (
                "faceAreaTmp",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            mesh_,
            dimensionedScalar(dimless, Zero)
        )
    );
    volScalarField& faceAreaTmpRef = faceAreaTmp.ref();

    //Calculate the possibilities for nucleation
    tmp<volScalarField> nuRate = nuRateCal();
    nuRateOut = nuRate.ref();
    Info<< "Max/Min nucleation rate(n/m2/s): " << max(nuRate.ref()).value() << ", " << min(nuRate.ref()).value()<< endl;

    //Loop through wall boundary patches to estimate nucleation
    const fvPatchList& patches = mesh_.boundary();
    labelList wallList;
    wallList.clear();

    scalarList faceAreaList;
    faceAreaList.clear();

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
                scalar faceCellAreai = mesh_.magSf().boundaryField()[patchI][pFaceI];
                faceAreaList.resize(faceAreaList.size()+1);
                faceAreaList[faceAreaList.size()-1]=faceCellAreai;
            }
        }
    }
    //Pout << "Wall list size: "<< wallList.size() << endl;

    std::random_device rd;  //Will be used to obtain a seed for the random number engine
    std::default_random_engine generator(rd());
    
    dimensionedScalar totSufArea(0.0);
    dimensionedScalar nuRateAccum(0.0);
    dimensionedScalar averNuRate(0.0);

    forAll(wallList, i)
    {
        //label poss = rand() % 1000000000000 + 1;
        //const labelListList& cellFaces = mesh_.cellFaces()[wallList[i]];

        //Mark out cells next to wall BC
        //wallMarkerRef[wallList[i]] = 1.0;

        const cell& faces = mesh_.cells()[wallList[i]];

        //scalar faceArea = mesh_.magSf()[wallList[i]];
        //Info << "Face ID: " << faces << endl;
        //forAll(faces, faceI)
        //{
        //    //Info << faceI << ": " << mesh_.magSf()[faces[faceI]] << endl;
        //    if(mesh_.magSf()[faceI]>faceArea)
        //    {
        //        faceArea = mesh_.magSf()[faceI];
        //    }
        //}

        //Info << faces << endl;

        //Info << mesh_.magSf()[faces[0]] << endl;

        faceAreaTmpRef[wallList[i]] = faceAreaList[i];
        //forAll(faces, faceI)
        //{
        //    //Info << faceI << ": " << mesh_.magSf()[faces[faceI]] << endl;
        //    if((mesh_.magSf()[faces[faceI]]>faceAreaTmpRef[wallList[i]]) && (mesh_.magSf()[faces[faceI]] < 2e-12))
        //    {
        //        faceAreaTmpRef[wallList[i]] = mesh_.magSf()[faces[faceI]];
        //    }
        //}
        //Pout << faceAreaTmpRef[wallList[i]] << endl;

        std::random_device rd;  //Will be used to obtain a seed for the random number engine
        std::default_random_engine generator(rd());
        scalar upperLim = 1.0/(faceAreaTmpRef[wallList[i]]+VSMALL);
        std::uniform_int_distribution<long long unsigned> dis(1, static_cast<long long unsigned>(upperLim));
        scalar randNum = dis(generator);

        if(debug_)
        {
            //Info<< "Max faceArea: " << faceArea << endl;
            //cout<< "Total cell counts per unit area: " << upper << endl;
            Info<< "Random integer generated: " << randNum << endl;
            Info<< "Nucleation rate: " << nuRate.ref()[wallList[i]] << endl;
            Info<< "Computational time: " << mesh_.time().value() << endl;
        }

        if((nuRate.ref()[wallList[i]])>randNum)
        {
            nuSitePerStepRef[wallList[i]] = 1.0;
            nuSite[wallList[i]] = 1.0;
            Cmask_[wallList[i]] = 1.0;
            if(alpha_[wallList[i]] > 0.5)
            {
                //nuTotal_.value() += 1.0;
                Info<< "New nucleation site found: " << nuRate.ref()[wallList[i]] << " > " << randNum << endl;
            }
        }
        // mark cells for check, remove in running code
        // nuSite[wallList[i]] = 1.0;

        // calculate the available surface area
        //totSufArea.value() += faceArea*alpha_[wallList[i]];
        //nuRateAccum.value() += nuRate.ref()[wallList[i]]*faceArea*alpha_[wallList[i]];
    }

    Info<< "Max/Min cell surface Area(m2): " << max(faceAreaTmpRef).value() << ", " << min(faceAreaTmpRef).value()<< endl;

    nuTotal_.value() += gSum((nuSitePerStepRef)());

    totSufArea.value() = gSum((faceAreaTmpRef*alpha_*pos(alpha_-0.1))())+VSMALL;
    nuRateAccum.value() = gSum((faceAreaTmpRef*alpha_*nuRate.ref()*pos(alpha_-0.1))());

    Info<< "Total nucleation sites: " << nuTotal_.value() << endl;
    averNuRate.value() = nuRateAccum.value()/totSufArea.value();

    // output available surface area
    if(Pstream::master() && mesh_.time().outputTime())
    {
        string fileName = "nucleationSurfaceArea.csv";
        std::ifstream file(fileName);
        std::fstream fout;
        fout.open(fileName, std::ios::out | std::ios::app);
        if (file.peek() == std::ifstream::traits_type::eof())
        {
            fout << "Time[s],totSufArea,averNuRate" << "\n";
        }
        fout << mesh_.time().timeName() << ", "
            << totSufArea.value() << ", "
            << averNuRate.value() << "\n";
    }

    if(debug_)
    {
        if (mesh_.time().outputTime())
        {
            faceAreaTmpRef.write();
            //wallMarkerRef.write();
        }
    }
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
        SItmpRef[cellI] = max(0.001, log10((pow(C_[cellI]/1000,2.0)+VSMALL)/Ksp));
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
            dimensionedScalar(dimless/dimArea/dimTime, Zero)
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

Foam::phaseChangeReaction::~phaseChangeReaction(){};

// ************************************************************************* //
