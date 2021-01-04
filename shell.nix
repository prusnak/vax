with import <nixpkgs> { };

let MyPython = python3.withPackages(ps: with ps; [ black dnachisel ipython ]);

in

stdenv.mkDerivation ({
  name = "vax";
  buildInputs = [
    MyPython
  ];
})
