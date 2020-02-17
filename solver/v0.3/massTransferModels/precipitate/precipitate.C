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

#include "precipitate.H"
#include "constants.H"
#include "fvcGrad.H"
#include "fvcSnGrad.H"
#include "fvcDiv.H"
#include "surfaceInterpolate.H"
#include "fvcReconstruct.H"
#include "fvm.H"
#include "zeroGradientFvPatchFields.H"

using namespace Foam::constant;

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Thermo, class OtherThermo>
Foam::meltingEvaporationModels::precipitate<Thermo, OtherThermo>
::precipitate
(
    const dictionary& dict,
    const phasePair& pair
)
:
    InterfaceCompositionModel<Thermo, OtherThermo>(dict, pair),
    C_("C", dimVelocity, dict),
    Tactivate_("Tactivate", dimTemperature, dict),
    Cactivate_("Cactivate", dimMoles/dimVolume, dict),
    Mv_("Mv", dimMass/dimMoles, dict),
    alphaMax_(dict.lookupOrDefault<scalar>("alphaMax", 1.0)),
    alphaMin_(dict.lookupOrDefault<scalar>("alphaMin", 0.5)),
    alphaSolidMin_(dict.lookupOrDefault<scalar>("alphaSolidMin", 0.5)),
    alphaRestMax_(dict.lookupOrDefault<scalar>("alphaRestMax", 0.01)),
    smoothSurface_(dict.lookupOrDefault<bool>("smoothSurface", false)),
    smoothAreaDensity_(dict.lookupOrDefault<scalar>("smoothAreaDensity", 1.0))
{
    if (this->transferSpecie() != "none")
    {
        word fullSpeciesName = this->transferSpecie();
        auto tempOpen = fullSpeciesName.find('.');
        const word speciesName(fullSpeciesName.substr(0, tempOpen));

        // Get the "to" thermo
        const typename OtherThermo::thermoType& toThermo =
            this->getLocalThermo
            (
                speciesName,
                this->toThermo_
            );

         // Convert from g/mol to Kg/mol
        Mv_.value() = toThermo.W()*1e-3;
    }


    if (Mv_.value() == -1)
    {
        FatalErrorInFunction
            << " Please provide the molar weight (Mv) of vapour [g/mol] "
            << abort(FatalError);
    }
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

template<class Thermo, class OtherThermo>
Foam::tmp<Foam::volScalarField>
Foam::meltingEvaporationModels::precipitate<Thermo, OtherThermo>
::Kexp(label variable, const volScalarField& field)
{
    if (true)
    {
        const volScalarField& to = this->pair().to();

        const volScalarField& from = this->pair().from();

        const fvMesh& mesh = this->mesh_;

        const volScalarField& C =
            mesh.lookupObject<volScalarField>("C").oldTime();

//        const volScalarField& areaDensityGrad = 
//            mesh.lookupObject<volScalarField>("areaDensityGrad");

        const dimensionedScalar convConst
        (
            "convConst",
            dimLength/dimTime,
            1.0
        );

        const dimensionedScalar convConstVol
        (
            "convConstVol",
            dimArea,
            1.0
        );

        const dimensionedScalar convConstLen
        (
            "convConstLen",
            dimless/dimLength,
            1.0
        );

        const dimensionedScalar convConstUnitLess
        (
            "convConstUnitLess",
            dimLength/dimLength,
            1.0
        );

        word fullSpeciesName = this->transferSpecie();
        auto tempOpen = fullSpeciesName.find('.');
        const word speciesName(fullSpeciesName.substr(0, tempOpen));


        const volVectorField gradFrom(fvc::grad(from));
        const volVectorField gradTo(fvc::grad(to));

        const volScalarField areaDensitySmooth
        (
            IOobject
            (
                "areaDensitySmooth",
                mesh.time().timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            mesh,
            dimensionedScalar("areaDensitySmoothInput", dimless/dimLength, smoothAreaDensity_)
        );

        Info<< "areaDensitySmooth min/max : "
            <<min(areaDensitySmooth).value()
            <<", "<< max(areaDensitySmooth).value()
            <<endl;

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
	        forAll(mesh.cellCells()[celli],cellj)
	        {
		
	    	    if (to[mesh.cellCells()[celli][cellj]] > alphaSolidMin_)
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
        forAll(mesh.boundaryMesh(), patchI)
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
                    mesh.time().timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                mesh,
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
                    mesh.time().timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                mesh,
                dimensionedScalar(dimless, Zero)
            )
        );
        volScalarField& tDelta = tTdelta.ref();
	    Info<< "calculate rhom&tDelta..."<<endl;

        if (sign(C_.value()) > 0)
        {
            //- convert mol/m3 to kg/m3 for OpenFOAM
            rhom =
                (C-Cactivate_)*Mv_*1e-3;

            tDelta = max
            (
                pos(C*Cmask - Cactivate_),//dimensionedScalar("C1", dimMoles/dimVolume, 1.0),
                dimensionedScalar("C0", dimless, Zero)
            );
//            Info<<"tDelta: " << tDelta << endl;
        }
        else
        {
            rhom =
                this->pair().to().rho()*this->pair().from().rho()
              / (this->pair().to().rho() - this->pair().from().rho());

            tDelta = max
            (
                Cmask*(Cactivate_ - C)/dimensionedScalar("C1", dimMoles/dimVolume, 1.0),
                dimensionedScalar("C0", dimless, Zero)
            );
        }
	    Info<< "calculate massFluxPrec..."<<endl;

        volScalarField massFluxPrec
        (
            "massFluxPrec",
            C_
          * rhom
          * tDelta
          * pos(from-0.01)
        );

        Info<< "massFluxPrec min/max: " << min(massFluxPrec).value() << " " << max(massFluxPrec).value() << endl;

        // 'from' phase normalization
        // WIP: Normalization could be convinient for cases where the area were
        // the source term is calculated is uniform
//        const dimensionedScalar Nl
//        (
//            gSum((areaDensity*mesh.V())())
//           /(
//               gSum
//               (
//                   ((areaDensity*from)*mesh.V())()
//               )
//             + dimensionedScalar("SMALL", dimless, VSMALL)
//            )
//        );


        if (mesh.time().outputTime())
        {
            areaDensityGrad.write();
            Cmask.write();
//            volScalarField mKGasDot
//            (
//                "mKGasDot",
//                massFluxPrec*areaDensity*Nl*from
//            );
//            mKGasDot.write();
        }

//    	Info<< "massFluxPrec dimensions: "<<massFluxPrec.dimensions()<<endl;

        
        //- Constrain the reaction rate to prevent negative concentration
        if(smoothSurface_)
        {
            forAll(mesh.C(), cellI)
            {
                if(massFluxPrec[cellI]*areaDensitySmooth[cellI]*mesh.time().deltaTValue()
                    > mag(C[cellI]-Cactivate_.value())*Mv_.value()*1e-3)
                    {
                        massFluxPrec[cellI] = mag(C[cellI]-Cactivate_.value())*Mv_.value()*1e-3
                                            / (mesh.time().deltaTValue()*areaDensitySmooth[cellI]+VSMALL);
                    }
            }   
        }
        else
        {
            forAll(mesh.C(), cellI)
            {
                if(massFluxPrec[cellI]*areaDensityGrad[cellI]*mesh.time().deltaTValue()
                    > mag(C[cellI]-Cactivate_.value())*Mv_.value()*1e-3)
                    {
                        massFluxPrec[cellI] = mag(C[cellI]-Cactivate_.value())*Mv_.value()*1e-3
                                            / (mesh.time().deltaTValue()*areaDensityGrad[cellI]+VSMALL);
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
    else
    {
        return tmp<volScalarField> ();
    }
}


template<class Thermo, class OtherThermo>
const Foam::dimensionedScalar&
Foam::meltingEvaporationModels::precipitate<Thermo, OtherThermo>
::Tactivate() const
{
    return Tactivate_;
}


// ************************************************************************* //
