newPackage(
     "BuchsbaumEisenbudMultipliers",
     Version => "0.3",
     Date => "August 24, 2020",
     AuxiliaryFiles => true,
     Authors => {{Name => "Federico Galetto",
     	       Email => "galetto.federico@gmail.com",
	       HomePage => "https://math.galetto.org"},
	   {Name => "Keller VandeBogert",
     	       Email => "kellerlv@math.sc.edu",
	       HomePage => "https://sites.google.com/view/kellervandebogert/home?authuser=0"}},
     HomePage => "https://github.com/galettof/BEMultipliers",
     Headline => "Buchsbaum-Eisenbud multipliers of free resolutions",
     PackageImports => {"SimpleDoc"}
     )


export {
    -- from multipliers.m2
    "aMultiplier",--method
    "cMultiplier",--method
    "ComputeRanks",--option
    "exteriorDuality",--method
    -- from complexes.m2
    "Kcomplex",
    "Lcomplex"
    }

load "./BuchsbaumEisenbudMultipliers/multipliers.m2"

load "./BuchsbaumEisenbudMultipliers/complexes.m2"

end



uninstallPackage "BuchsbaumEisenbudMultipliers"
restart
installPackage "BuchsbaumEisenbudMultipliers"
installPackage("BuchsbaumEisenbudMultipliers",RemakeAllDocumentation=>true)
check BuchsbaumEisenbudMultipliers
