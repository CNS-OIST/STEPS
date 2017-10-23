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

    steps-py3_build = steps_build.override {
        python = python3;
        pythonPackages = python3Packages;
        cython = python3Packages.cython;
        numpy = python3Packages.numpy;
    };
in
{
    steps = steps_build;
    steps-py3 = steps-py3_build;
    steps-src-validation = steps-validation.override {
        steps = steps_build;
    };
    steps-src-py3-validation = steps-validation.override {
        steps = steps-py3_build;
    };
}
