/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | OpenQBMM - www.openqbmm.org
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Code created 2015 by Matteo Icardi and 2017 by Alberto Passalacqua
    Contributed 2018-07-31 to the OpenFOAM Foundation
    Copyright (C) 2018 OpenFOAM Foundation
    Copyright (C) 2019 Alberto Passalacqua
-------------------------------------------------------------------------------
2017-03-28 Alberto Passalacqua: Adapted to single scalar calculation.
2020-03-12 Alberto Passalacqua: Added access to fvMesh.
-------------------------------------------------------------------------------
License
    This file is derivative work of OpenFOAM.

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
    Foam::populationBalanceSubModels::growthModels::constantGrowth

Description
    Constant growth model. The value of the constant is read from dictionary
    or assumed equal to 1. The minimum and maximum value of the abscissa to
    which the model is applied can be defined by the user or it is assumed to
    be 1.0e-10 and 1.0e10.

SourceFiles
    constantGrowth.C

\*---------------------------------------------------------------------------*/

#ifndef classicGrowth_H
#define classicGrowth_H

#include "growthModel.H"
#include "nucleationModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace populationBalanceSubModels
{
namespace growthModels
{

/*---------------------------------------------------------------------------*\
                    Class constantGrowth Declaration
\*---------------------------------------------------------------------------*/

class classicGrowth
:
    public growthModel
{
private:

    // Private data

        //- Minimum value of the abscissa to apply the model
        scalar minAbscissa_;

        //- Maximum value of the abscissa to apply the model
        scalar maxAbscissa_;

        //- 分子摩尔质量
        dimensionedScalar moleDensity_;

        //- 密度
        dimensionedScalar particleRho_;

        //- Concentration change in solution
        //- need to be defined previously
        const volScalarField& dcdt_;

        //- total number of particles in cell 
        //- need to be defined previously
        //- celli中颗粒总数量
        const volScalarField& totalN_;

        //- Nucleation model
        autoPtr
        <
            Foam::populationBalanceSubModels::nucleationModel
        > nucleationModel_;

public:

        //- Runtime type information
        TypeName("classicGrowth");


    // Constructors

        //- Construct from components
        classicGrowth
        (
            const dictionary& dict,
            const fvMesh& mesh
        );


    //- Destructor
    virtual ~classicGrowth();


    // Member Functions

        //- Growth model, not used
        virtual scalar Kg
        (
            const scalar& abscissa,
            const bool lengthBased,
            const label environment = 0
        ) const;

        //- private growth model actually used
        scalar Kg_new
        (
            const scalar& abscissa,
            const label celli,
            const bool lengthBased,
            const label environment = 0
        );

        //- Return phase space convection
        virtual scalar phaseSpaceConvection
        (
            const labelList& momentOrder,
            const label celli,
            const scalarQuadratureApproximation& quadrature
        );

        //- Return phase space convection
        virtual scalar phaseSpaceConvection
        (
            const labelList& momentOrder,
            const label celli,
            const velocityQuadratureApproximation& quadrature
        );

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace growthModels
} // End namespace populationBalanceSubModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //


