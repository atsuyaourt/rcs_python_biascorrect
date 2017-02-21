require(dplyr)
require(qmap)

fit_qmap <- function(df) {
    if(sum(!is.na(df$obs)) > 10) {
        return(qmap::fitQmapQUANT(df$obs, df$mod))
    }
}

do_qmap <- function(df, fit) {
    if(missing(fit)) {
        df$bc <- df$mod
    } else {
        df$bc <- qmap::doQmap(df$mod, fit)
    }
    return(df)
}
