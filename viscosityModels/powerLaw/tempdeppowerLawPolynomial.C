/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2020 OpenFOAM Foundation
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

#include "tempdeppowerLawPolynomial.H"
#include "addToRunTimeSelectionTable.H"
#include "surfaceFields.H"
#include "fvCFD.H"
//#include <csignal> // for gdb
// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace viscosityModels
{
    defineTypeNameAndDebug(tempdeppowerLawPolynomial, 0);

    addToRunTimeSelectionTable
    (
        viscosityModel,
        tempdeppowerLawPolynomial,
        dictionary
    );
}
}


// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::viscosityModels::tempdeppowerLawPolynomial::calcNu() const
{
//Info<< "it's now calculating viscosity "<< nl << endl;
//added for temperature dependent
const volScalarField& T = U_.mesh().lookupObject<volScalarField>("T");
//added for temperature dependent

/*
Info<< "(k_-kslope_*(T-Tbase_))" << (k_-kslope_*(T-Tbase_)) << nl << endl;
Info<< "dimensionedScalar(dimTime, 1.0)*strainRate()" << dimensionedScalar(dimTime, 1.0)*strainRate() << nl << endl;
Info<< "dimless "<< dimless << nl << endl;
Info<< "n_.value() - scalar(1)"<< (n_.value() - scalar(1)) << nl << endl;
*/
/*
	volTensorField gradUx_yz = fvc::grad(U_);
	forAll(gradUx_yz, i)
	{
//	Info<< "i" << i<< endl;	
		forAll(gradUx_yz[i],j)
		{
//			Info<< "j" << gradUx_yz[i][j]<< endl;
			if ((j == 3) or (j == 6))
			{continue;}
			else
			{gradUx_yz[i][j] = 0;}
		}
	}
	volScalarField strainRate_new = dimensionedScalar(dimTime, 1.0)*Foam::sqrt(2.0)*mag(symm(gradUx_yz)) ;
	
	Info<< "gradUx_yz"<< gradUx_yz << endl;
	Info<< "shear rate"<< Foam::sqrt(2.0)*mag(symm(gradUx_yz)) << endl;
*/	
	
    return max
    (
        nuMin_,
        min
        (
            nuMax_,
            //changed for temperature dependent
//            (k_-kslope_*(T-Tbase_))*pow
            (k0+k1*(T-dimensionedScalar(dimTemperature, 273.15))+k2*pow((T-dimensionedScalar(dimTemperature, 273.15)),2))*k_rho*pow
            //changed for temperature dependent
            (
                max
                (
//                    dimensionedScalar(dimTime, 1.0)*strainRate(),
//			strainRate_new,
			dimensionedScalar(dimTime, 1.0)*Foam::sqrt(2.0)*mag(symm(fvc::grad(U_))),			

                    dimensionedScalar(dimless, small)
                ),
                (n0+n1*T+n2*pow(T,2) - 1)
            )
        )
    );
}



// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::viscosityModels::tempdeppowerLawPolynomial::tempdeppowerLawPolynomial
(
    const word& name,
    const dictionary& viscosityProperties,
    const volVectorField& U,
    const surfaceScalarField& phi
)
:
    viscosityModel(name, viscosityProperties, U, phi),
    tempdeppowerLawPolynomialCoeffs_(viscosityProperties.optionalSubDict(typeName + "Coeffs")),
    k0("k0", dimViscosity, tempdeppowerLawPolynomialCoeffs_),
    k1("k1", dimViscosity/dimTemperature, tempdeppowerLawPolynomialCoeffs_),    
    k2("k2", dimViscosity/dimTemperature/dimTemperature, tempdeppowerLawPolynomialCoeffs_),     
    k_rho("k_rho", dimless, tempdeppowerLawPolynomialCoeffs_),   
       
    n0("n0", dimless, tempdeppowerLawPolynomialCoeffs_),
    n1("n1", dimless/dimTemperature, tempdeppowerLawPolynomialCoeffs_),
    n2("n2", dimless/dimTemperature/dimTemperature, tempdeppowerLawPolynomialCoeffs_),
            
    //added for temperature dependent
//    kslope_("kslope", dimViscosity/dimTemperature, tempdeppowerLawPolynomialCoeffs_),
//    Tbase_("Tbase", dimTemperature, tempdeppowerLawPolynomialCoeffs_),
    //added for temperature dependent
    
    nuMin_("nuMin", dimViscosity, tempdeppowerLawPolynomialCoeffs_),
    nuMax_("nuMax", dimViscosity, tempdeppowerLawPolynomialCoeffs_),
    nu_
    (
        IOobject
        (
            name,
            U_.time().timeName(),
            U_.db(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        calcNu()
    )
{}//print nu, Info<< "calcNu()"<< calcNu() << nl << endl;


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

bool Foam::viscosityModels::tempdeppowerLawPolynomial::read
(
    const dictionary& viscosityProperties
)
{
//    std::raise(SIGINT); // for gdb
//    Info<< "read read read read"<< nl << endl;
//    std::raise(SIGINT); // for gdb
    
    viscosityModel::read(viscosityProperties);

    tempdeppowerLawPolynomialCoeffs_ = viscosityProperties.optionalSubDict(typeName + "Coeffs");

    tempdeppowerLawPolynomialCoeffs_.lookup("k0") >> k0;
    tempdeppowerLawPolynomialCoeffs_.lookup("k1") >> k1;    
    tempdeppowerLawPolynomialCoeffs_.lookup("k2") >> k2; 
    tempdeppowerLawPolynomialCoeffs_.lookup("k_rho") >> k2;        
    tempdeppowerLawPolynomialCoeffs_.lookup("n0") >> n0;
    tempdeppowerLawPolynomialCoeffs_.lookup("n1") >> n1;
    tempdeppowerLawPolynomialCoeffs_.lookup("n2") >> n2;     
       
    //added for temperature dependent
//    tempdeppowerLawPolynomialCoeffs_.lookup("kslope") >> kslope_;
//    tempdeppowerLawPolynomialCoeffs_.lookup("Tbase") >> Tbase_;        
    //added for temperature dependent
    
    tempdeppowerLawPolynomialCoeffs_.lookup("nuMin") >> nuMin_;
    tempdeppowerLawPolynomialCoeffs_.lookup("nuMax") >> nuMax_;

    return true;
}


// ************************************************************************* //
