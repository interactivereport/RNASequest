suppressMessages(require(Hmisc))
Hmisc.rcorr <- function (x, y, type = "pearson")
{
    type <- match.arg(type)#c("pearson", "spearman")
    if (!missing(y))
        x <- cbind(x, y)
    x[is.na(x)] <- 1e+50
    storage.mode(x) <- "double"
    p <- as.integer(ncol(x))
    if (p < 1)
        stop("must have >1 column")
    n <- as.integer(nrow(x))
    if (n < 5)
        stop("must have >4 observations")
    h <- .Fortran(Hmisc:::F_rcorr, x, n, p, itype = as.integer(1 + (type =="spearman")),
                  hmatrix = double(p * p), npair = integer(p * p), double(n),
                  double(n), double(n), double(n), double(n),
                  integer(n))
    npair <- matrix(h$npair, ncol = p)
    h <- matrix(h$hmatrix, ncol = p)
    h[h > 1e+49] <- NA
    nam <- dimnames(x)[[2]]
    dimnames(h) <- list(nam, nam)
    dimnames(npair) <- list(nam, nam)
    #https://stackoverflow.com/questions/63994852/r-in-sqrt1-h-h-nans-produced-from-within-rcorr-full-sample-data-avai
    #P <- matrix(2 * (1 - pt(abs(h) * sqrt(npair - 2)/sqrt(1 - h * h), npair - 2)), ncol = p)
    P <- matrix(2 * (1 - pt(abs(h) * sqrt(npair - 2)/max(0, 1-h^2), npair - 2)), ncol = p)
    P[abs(h) == 1] <- 0
    diag(P) <- NA
    dimnames(P) <- list(nam, nam)
    structure(list(r = h, n = npair, P = P), class = "rcorr")
}
