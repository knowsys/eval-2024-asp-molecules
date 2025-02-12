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
          pkgs.python311Packages.networkx
        ];
        shellHook = ''
          export LD_LIBRARY_PATH="${pkgs.stdenv.cc.cc.lib}/lib"
        '';
      };
    };
}

