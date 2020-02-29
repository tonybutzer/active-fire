\name{l8fire_detect.R}
\alias{l8fire_detect.R}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
l8fire_detect.R(b7, b6, b5, b4, b3, b2, fname = "out.tiff")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{b7}{
%%     ~~Describe \code{b7} here~~
}
  \item{b6}{
%%     ~~Describe \code{b6} here~~
}
  \item{b5}{
%%     ~~Describe \code{b5} here~~
}
  \item{b4}{
%%     ~~Describe \code{b4} here~~
}
  \item{b3}{
%%     ~~Describe \code{b3} here~~
}
  \item{b2}{
%%     ~~Describe \code{b2} here~~
}
  \item{fname}{
%%     ~~Describe \code{fname} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (b7, b6, b5, b4, b3, b2, fname = "out.tiff") 
{
    require(raster)
    require(rgdal)
    dyn.load("l8fire_detection.dll")
    out <- band7
    out[, ] <- 0
    jcol <- dim(band7)[2]
    irow <- dim(band7)[1]
    bs <- .C("l8fire_detection", jcol = as.integer(jcol), irow = as.integer(irow), 
        data7 = as.double(b7[, ]), data6 = as.double(b6[, ]), 
        data5 = as.double(b5[, ]), data4 = as.double(b4[, ]), 
        data3 = as.double(b3[, ]), data2 = as.double(b2[, ]), 
        out = as.integer(out[, ]))
    dyn.unload("l8fire_detection.dll")
    out[, ] <- bs$out
    print(unique(bs$out))
    print(fname)
    writeRaster(out, fname, format = "GTiff", overwrite = T, 
        datatype = "INT1U")
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
