\name{string functions}

\alias{string_functions}
\alias{string functions}
\alias{collapse.vec}
\alias{upper.first}

\title{
Functions for Creating or Manipulating Strings
}

\description{
See Details for a description of the individual functions.
}

\usage{
collapse.vec(x, space=TRUE, quote=FALSE)
upper.first(s)
}

\arguments{
\item{x}{a vector.}
\item{space}{logical; if \code{TRUE}, items will be separated by spaces after commas.}
\item{quote}{logical; if \code{TRUE}, items will be quoted.}
\item{s}{a string.}
}

\details{
The function \code{collapse.vec} returns a string representing a vector, intended for printing or for entry as a value.  The elements of the vector (``items'') are separated by commas and enclosed by parentheses.

The function \code{upper.first} returns a version of the string \code{s} with the first letter (if any) in the string in upppercase.
}

\value{A string.}
\author{
Daniel Dvorkin
}

\examples{
collapse.vec(1:3) # "(1, 2, 3)"

upper.first(collapse.vec(c("spam", "eggs"))) # "(Spam, eggs)"
}

\keyword{manip}
\keyword{character}
