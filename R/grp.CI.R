grp.CI <-
function (maUp, maLow, position, barLen, col=NULL,...) 
{
    vLen = length(position)
    lapply(1:vLen, FUN = function(i, up, low, xPos, barLen, col,...) {
        segments((xPos - barLen/2)[i], low[i], (xPos + barLen/2)[i], 
            low[i], col,...)
        segments(xPos[i], low[i], xPos[i], up[i],col, ...)
        segments((xPos - barLen/2)[i], up[i], (xPos + barLen/2)[i], 
            up[i],col, ...)
    }, up = maUp, low = maLow, xPos = position, barLen = barLen, col)
    return(NULL)
}

