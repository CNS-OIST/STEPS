with import <BBPpkgs> { };


let 
    steps_build = steps.overrideDerivation (oldAttr: rec {
        name = "steps-DEV_ENV";
        src = ./.;

        makeFlags = [ "VERBOSE=1" ];

        doCheck = true;

        checkPhase = [ ''
		export LD_LIBRARY_PATH=$PWD/src:$LD_LIBRARY_PATH
		ctest -V	
        ''];
    });
in
{ steps = steps_build; }
