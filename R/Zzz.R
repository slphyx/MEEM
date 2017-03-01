.onAttach <- function(libname, pkgname) {
  packageStartupMessage(
    paste("\nThis is MEEM version", packageVersion("MEEM"), "\n"))
}
