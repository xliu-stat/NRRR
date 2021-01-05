utils::globalVariables(c("Var1", "Var2", "value"))


.onUnload <- function (libpath) {
  library.dynam.unload("NRRR", libpath)
}
