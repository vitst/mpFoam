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

Class
    Foam::meltingEvaporationModels::precipitate

Description

Usage
    Example usage:
    \verbatim
        massTransferModel
        (
            (solid to liquid)
            {
                type            precipitate;
                C               40;
                Cactivate       302.78;
            }
        );
    \endverbatim

    Where:

    \table
        Property    | Description            | Required    | Default value
        Cactivate   | Activation temperature | yes
        C           | Model constant         | yes
    \endtable

SourceFiles
    precipitate.C

\*---------------------------------------------------------------------------*/

#ifndef meltingEvaporationModels_precipitate_H
#define meltingEvaporationModels_precipitate_H

#include "InterfaceCompositionModel.H"


// * * * * * * * * * * * * * * * * *  * * * * * * * * * * * * * * * * * * * *//

namespace Foam
{
namespace meltingEvaporationModels
{

/*---------------------------------------------------------------------------*\
                                  Class precipitate
\*---------------------------------------------------------------------------*/

template<class Thermo, class OtherThermo>
class precipitate
:
    public InterfaceCompositionModel<Thermo, OtherThermo>
{
    // Private data

        //- Condensation coefficient [1/s]
        dimensionedScalar C_;

        //- Phase transition concentration 
        const dimensionedScalar Tactivate_;

        //- Phase minimum value for activation
        scalar alphaMin_;


public:

    //- Runtime type information
    TypeName("precipitate");


    // Constructors

        //- Construct from components
        precipitate
        (
            const dictionary& dict,
            const phasePair& pair
        );


    //- Destructor
    virtual ~precipitate() = default;


    // Member Functions

        //- Explicit mass transfer coefficient
        virtual tmp<volScalarField> Kexp
        (
            label variable,
            const volScalarField& field
        );

        //- Return T transition between phases
        virtual const dimensionedScalar& Tactivate() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace meltingEvaporationModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "precipitate.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
