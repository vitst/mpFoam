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
#include "addToRunTimeSelectionTable.H"


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Thermo, class OtherThermo>
Foam::meltingEvaporationModels::precipitate<Thermo, OtherThermo>::precipitate
(
    const dictionary& dict,
    const phasePair& pair
)
:
    InterfaceCompositionModel<Thermo, OtherThermo>(dict, pair),
    C_("C", inv(dimTime), dict),
    Tactivate_("Tactivate", dimTemperature, dict),
    Cactivate_("Cactivate", dimMoles/dimVolume, dict),
    alphaMin_(dict.lookupOrDefault<scalar>("alphaMin", 0))
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

template<class Thermo, class OtherThermo>
Foam::tmp<Foam::volScalarField>
Foam::meltingEvaporationModels::precipitate<Thermo, OtherThermo>::Kexp
(
    label variable,
    const volScalarField& refValue
)
{
    const fvMesh& mesh = this->mesh_;
//    const volScalarField& Ctest =
//          mesh.lookupObject<volScalarField>("C").oldTime();


    Info << endl << "C field min/max: " << min(refValue).value()<< " / " << max(refValue).value() << endl << endl;

    tmp<volScalarField> rate;

        volScalarField from
        (
            min(max(this->pair().from(), scalar(0)), scalar(1))
        );

        volScalarField to
        (
            min(max(this->pair().to(), scalar(0)), scalar(1))
        );


        if (sign(C_.value()) > 0)
        {
            rate = 
            (
                C_
              * from
              * this->pair().from().rho() 
              * (refValue.oldTime() - Cactivate_)
              * pos(from - alphaMin_)
              * pos(refValue.oldTime() - Cactivate_)/Cactivate_
            );
        }
        else
        {
            rate =
            (
               -C_
              * from
              * this->pair().from().rho() 
              * pos(from - alphaMin_)
              * (Cactivate_ - refValue.oldTime())
              * pos(Cactivate_ - refValue.oldTime())/Cactivate_
            );

        }
    
    label cCount = 0;

    // Only retain reaction rate next to solid phase cell
    forAll(mesh.C(),i)
    {
    	bool flag = true;
    	forAll(mesh.cellCells()[i], j)
    	{ 
    		if(to[mesh.cellCells()[i][j]]>0.1)
    		{
			cCount += 1;
			flag = false;
			break;
    		}
    	}
    	if(flag)
    	{
    		rate.ref()[i] = 0.0;
	}
    }

    Pout<< "Find total solid surface cell: " << cCount <<endl;
    Info<< "Reaction rate calculation: " << min(rate.ref()).value() << " / " << max(rate.ref()).value()<<endl;

    return rate;
}


template<class Thermo, class OtherThermo>
const Foam::dimensionedScalar&
Foam::meltingEvaporationModels::precipitate<Thermo, OtherThermo>::Tactivate() const
{
    return Tactivate_;
}


// ************************************************************************* //
