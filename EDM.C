/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2015 OpenFOAM Foundation
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

#include "EDM.H"
#include "turbulenceModel.H"
#include "turbulentFluidThermoModel.H"
#include "volFields.H"

namespace Foam
{
namespace combustionModels
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CombThermoType, class ThermoType>
EDM<CombThermoType, ThermoType>::EDM
(
    const word& modelType,
    const fvMesh& mesh,
    const word& combustionProperties,
    const word& phaseName
)
:
    singleStepCombustion<CombThermoType, ThermoType>
    (
        modelType,
        mesh,
        combustionProperties,
        phaseName
    ),

    C_(readScalar(this->coeffs().lookup("C_EDM"))),

    delta_(this->coeffs().lookupOrDefault("delta", 1.0)),

    CDelta(readScalar(this->coeffs().lookup("C_Delta"))),

    chist_ref
    (
        dimensionedScalar("chist_ref",
        dimensionSet(0,0,-1,0,0,0,0),
        readScalar(this->coeffs().lookup("chist_ref")))
    ),

    D0
    (
        dimensionedScalar("D0", 
        dimensionSet(0,2,-1,0,0,0,0),
        readScalar(this->coeffs().lookup("D0")))
    ),

    CDiff
    (
        IOobject
        (
            "CDIff",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("zero", dimless, 0.0)
    ),    

    rtTurb
    (
        IOobject
        (
            "rtTurb",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("zero", dimensionSet(0, 0, -1, 0 , 0, 0, 0), 0.0)
    ),  

    rtDiff
    (
        IOobject
        (
            "rtDiff",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("zero", dimensionSet(0, 0, -1, 0 , 0, 0, 0), 0.0)
    ),

    deltaCell
    (
        IOobject
        (
            "deltaCell",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("deltaCell", dimensionSet(0, 1, 0, 0 , 0, 0, 0), delta_)
    )       

{

}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class CombThermoType, class ThermoType>
EDM<CombThermoType, ThermoType>::~EDM()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

template<class CombThermoType, class ThermoType>
void EDM<CombThermoType, ThermoType>::correct()
{
    this->wFuel_ ==
        dimensionedScalar("zero", dimMass/pow3(dimLength)/dimTime, 0.0);

    if (this->active())
    {
        this->singleMixturePtr_->fresCorrect();

        const label fuelI = this->singleMixturePtr_->fuelIndex();

        const volScalarField& YFuel = this->thermoPtr_->composition().Y()[fuelI];
        const volScalarField& YO2 = this->thermoPtr_->composition().Y("O2");

        const dimensionedScalar s = this->singleMixturePtr_->s();

        scalar YO2Inf = readScalar(this->thermoPtr_->lookup("YO2Inf"));
        scalar YFInf = readScalar(this->thermoPtr_->lookup("YFInf"));


        if (YFuel.db().foundObject<compressible::LESModel>(turbulenceModel::propertiesName))
        {
            const compressible::LESModel& lesModel =
                YFuel.db().lookupObject<compressible::LESModel>
                (
                 turbulenceModel::propertiesName
                );

            deltaCell == lesModel.delta();   
        }

        //Inverse of turbulent time-scale
        rtTurb == C_*this->turbulence().epsilon()/
              max(this->turbulence().k(),
              dimensionedScalar("SMALL",dimVelocity*dimVelocity,SMALL));

        //Inverse of diffusion time-scale
        CDiff == 2.0/(CDelta*sqrt(2.0*D0)) * sqr(s*YFInf+YO2Inf)/(s*YFInf*YO2Inf) 
            * sqrt(chist_ref) * deltaCell; 

        rtDiff == CDiff * this->thermoPtr_->alpha()/this->rho() /sqr(deltaCell);    

        this->wFuel_ ==
                  this->rho()
                * min(YFuel, YO2/s.value())
                * max(rtTurb,rtDiff);
            
    }
}


template<class CombThermoType, class ThermoType>
bool EDM<CombThermoType, ThermoType>::read()
{
    if (singleStepCombustion<CombThermoType, ThermoType>::read())
    {
        this->coeffs().lookup("C") >> C_ ;
        return true;
    }
    else
    {
        return false;
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace combustionModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
