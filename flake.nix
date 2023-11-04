{
  description = "Flake for development of the ddml R package";
  inputs.nixpkgs.url = "github:NixOS/nixpkgs/nixpkgs-unstable";
  inputs.flake-utils.url = "github:numtide/flake-utils";

  outputs = { self, nixpkgs, flake-utils }:
    flake-utils.lib.eachDefaultSystem (system: let
      
      pkgs = nixpkgs.legacyPackages.${system};

      # R packages
      my-R-packages = with pkgs.rPackages; [
        devtools
        pkgdown
        AER
        MASS
        Matrix
        nnls
        quadprog
        glmnet
        ranger
        xgboost
        sandwich
        covr
        testthat
        knitr
        markdown
        rmarkdown
        readstata13
      ];
      my-R = [pkgs.R my-R-packages];

    in {
      devShells.default = pkgs.mkShell {
        nativeBuildInputs = [ pkgs.bashInteractive ];
        buildInputs = [ 
          my-R
          pkgs.rstudio
          pkgs.quarto # needed for rstudio
        ];
       };
    });
}