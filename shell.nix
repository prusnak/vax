with import <nixpkgs> { };

let MyPython = python3.withPackages(ps: with ps; [ dnachisel ]);

in

stdenv.mkDerivation ({
  name = "vax";
  buildInputs = [
    MyPython
  ];
})
