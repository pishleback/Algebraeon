{ pkgs, ... }:

{
  languages.rust = {
    enable = true;
    channel = "stable";
  };

  packages = with pkgs; [
    clang
    llvmPackages.bintools
  ];
}