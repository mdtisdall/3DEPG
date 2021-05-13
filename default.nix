# default.nix

    with import <nixpkgs> {};
    stdenv.mkDerivation {
      name = "dev-environment"; 
      buildInputs = [
				gcc
				gdb
				mkl
      ];
    }
