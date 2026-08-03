{ pkgs, ... }:

{
  languages.rust = {
    enable = true;
    channel = "stable";
  };

  packages = with pkgs; [
    clang
    llvmPackages.bintools
    mdbook
  ];

  processes.guide = {
    exec = "mdbook serve --port 3000";
    cwd = "./guide";
    ready = {
      http.get = {
        host = "localhost";
        port = 3000;
        path = "/";
      };
      period = 5;
    };
  };
}