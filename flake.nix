{
  description = "Flake for development of the ddml R package";
  inputs.nixpkgs.url = "github:NixOS/nixpkgs/nixpkgs-unstable";
  inputs.flake-utils.url = "github:numtide/flake-utils";

  outputs = { self, nixpkgs, flake-utils }:
    flake-utils.lib.eachDefaultSystem (system: let

      pkgs = nixpkgs.legacyPackages.${system};

      # Install local version of ddml
      ddml = pkgs.rPackages.buildRPackage {
        name = "ddml";
        src = ./.;
        # ddml dependencies
        propagatedBuildInputs = with pkgs.rPackages; [AER MASS Matrix nnls quadprog glmnet ranger xgboost sandwich];
      };

      # R packages
      my-R-packages = with pkgs.rPackages; [
        # this package
        ddml
        # dependencies
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
        # packages used in the vignettes
        did
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
