# first and last
.First.lib <- function(lib, pkg)
{
#require("epitools","TeachingDemos","tcltk",quietly=TRUE,warn.conflicts=FALSE)
# if (.Platform$OS.type=="windows")
see <- packageDescription(pkg,fields="Version")
cat("'DiagnosisMed' library",see," loaded\n",sep=" ")
}

.Last.lib <- function(libpath)
{
# nothing so far
}
