{
  description = "arizona searches by imputation";
  nixConfig.bash-prompt = "\[az-search\]$ ";

  inputs = {
    nixpkgs.url = "github:NixOS/nixpkgs/release-23.05";
    flake-utils.url = "github:numtide/flake-utils";
    impute.url = "github:EgorLappo/arizona-searches-beagle-impute";
    matches.url = "github:EgorLappo/arizona-searches-compute-matches";
  };

  outputs = { self, nixpkgs, flake-utils, impute, matches }:
    flake-utils.lib.eachDefaultSystem (system:
      let
        pkgs = (import nixpkgs) {
          inherit system;
        };

        dontTestPackage = drv: drv.overridePythonAttrs (old: { doCheck = false; });
        python-env = pkgs.python311.withPackages (ps: with ps; [
          tqdm
          numpy
          scipy
          pandas
          matplotlib
          (dontTestPackage seaborn) # tests fail on darwin due to different numerical results on intel vs ARM
        ]);
        pypath = "${python-env}/bin/python";

        R-env = pkgs.rWrapper.override {
          packages = with pkgs.rPackages; [
            readr
            dplyr
            data_table
            purrr
            doParallel
          ];
        };
        rscriptpath = "${R-env}/bin/Rscript";

        prepare = pkgs.writeShellScriptBin "prepare" ''
          ${pkgs.cowsay}/bin/cowsay Thanks for trying to replicate this analysis!

          cd raw_data/vcf

          echo "writing STR loci info (name and chromosome)"
          ls -l *.vcf | awk '{print $NF}' | cut -d "_" -f 1 > ../STR_loci.csv

          touch ../STR_loci_chr.csv
          for f in *.vcf; do
              cat $f | grep -v "^#" | head -n1 | cut -f 1 >> ../STR_loci_chr.csv
          done

          echo "indexing vcf files"
          for f in *.vcf; do
              bcftools view $f -Oz -o "$f.gz" 
              bcftools index "$f.gz" 
          done

          echo "writing sample ids"
          bcftools query -l CSF1PO_halfwindow500000WithCODIS.vcf > ../sample_ids.csv

          cd ../..

          mkdir -p processed_data
          mkdir -p beagle && cd beagle
          echo "downloading beagle"
          wget -q https://faculty.washington.edu/browning/beagle/beagle.22Jul22.46e.jar
          wget -q http://faculty.washington.edu/browning/beagle/bref3.22Jul22.46e.jar
          mkdir -p plink && cd plink
          wget -q https://bochet.gcc.biostat.washington.edu/beagle/genetic_maps/plink.GRCh37.map.zip 
          ouch decompress plink.GRCh37.map.zip
          mv plink.GRCh37/* .
          rm -rf plink.GRCh37

          cd ../..
          mkdir -p imputation_results

          echo "done with the setup stage!"
        '';

        impute-bin = impute.packages."${system}".default;

        matches-bin = matches.defaultPackage."${system}";

        run-analysis = pkgs.writeShellScriptBin "run-analysis" ''
          time prepare
          time impute
          time matches
          
          cd make_figures
          time ${rscriptpath} --vanilla run_tests.R
          time ${pypath} make_figures.py
        '';

      in
      rec {
        # make a default shell with default packages from impute and compute-matches flakes
        devShells.default = with pkgs; mkShell {
          name = "arizona_searches";

          LC_ALL = "C";

          buildInputs = [
            prepare
            run-analysis
            matches-bin
            impute-bin
            jdk
            bcftools
            python-env
            R-env
            wget
            ouch
          ];

          shellHook = ''
            alias Rs="${rscriptpath} --vanilla"
            ${pkgs.cowsay}/bin/cowsay -f small hi! welcome to the az-search shell!
          '';
        };
      }
    );
}
