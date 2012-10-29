
".First.lib" <- function(lib, pkg) library.dynam("CHsharp", pkg, lib)

"sharp3d" <-
function (x, y, z, hspace=1, htime=1, v=1) 
{
n <- length(x)
xsharp <- x
ysharp <- y
zsharp <- z
hsharp <- hspace
z <- .Fortran("sharp3d", as.integer(n), as.double(hsharp), as.double(htime),
as.double(x), as.double(y), as.double(z), as.double(xsharp),
as.double(ysharp), as.double(zsharp), as.integer(v), PACKAGE="CHsharp")

names(z) <- c("n", "hsharp", "htime", "x", "y", "z", "xsharp", "ysharp", "zsharp", "v")
data.frame(x.sharp=z$xsharp, y.sharp=z$ysharp, z.sharp=z$zsharp)
}

"sharp3dB" <-
function (x, y, z, hspace=1, htime=1, v=1) 
{
n <- length(x)
xsharp <- x
ysharp <- y
zsharp <- z
hsharp <- hspace
z <- .Fortran("sharp3dB", as.integer(n), as.double(hsharp), 
as.double(htime),
as.double(x), as.double(y), as.double(z), as.double(xsharp),
as.double(ysharp), as.double(zsharp), as.integer(v), PACKAGE="CHsharp")

names(z) <- c("n", "hsharp", "htime", "x", "y", "z", "xsharp", "ysharp", "zsharp", "v")
data.frame(x.sharp=z$xsharp, y.sharp=z$ysharp, z.sharp=z$zsharp)
}

