
Use `R CMD check` to  check this package, not **RStudio**! 

Indeed, in the check, **RStudio** deletes the `inst/doc` directory
without warning. However this directory used here is conforming to
section 1.4 in the [Writing R Extensions
manual](https://cran.r-project.org/doc/manuals/R-exts.html)

"*In addition to the help files in Rd format, R packages allow the
inclusion of documents in arbitrary other formats. The standard
location for these is subdirectory `inst/doc` of a source package, ...*"

In the present case `inst/doc` contains `RenextGuide.pdf` as well as
its source files.
