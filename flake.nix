{
  description = "Flake for genmol evaluation (FPG) on Wille system";

  inputs.nixpkgs.url = "github:NixOS/nixpkgs/nixpkgs-unstable";

  outputs = { self, nixpkgs }:
    let
      pkgs = import nixpkgs
        {
          system = "x86_64-linux";
        };
    in
    {
      devShells.x86_64-linux.default = pkgs.mkShell {
        packages = [
          pkgs.vim
          pkgs.clingo
          pkgs.python311
          pkgs.python311.pkgs.matplotlib
          pkgs.python311Packages.dataclasses-json
          pkgs.libstdcxx5
          pkgs.wineWowPackages.minimal
        ];

        shellHook = ''
          alias piclasp="python2.7 $PWD/piclasp-1.2.1-x86_64-linux/piclasp.py";
          alias breakID="$PWD/BreakID-2.5";
          alias python="python3";
          alias molgen="wine $PWD/molgen50/mgen.exe";
          alias > $PWD/.bash_aliases
        '';
      };
    };
}

