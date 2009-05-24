###
### From spida-test
###

# Problem: tab not showing NaNs

tab <- function(x,...) UseMethod('tab')


tab.default <-
function (..., total.margins = TRUE, useNA = "ifany")
{
    aa <- list(...)
    if (length(aa) == 1 && is.list(aa[[1]])) {
        aa <- c(aa[[1]], total.margins = total.margins, useNA = useNA)
        return(do.call("tab", aa[[1]]))
    }
    if (is.null(names(aa))) {
        nns = names(match.call())
        names(aa) = nns[2:(1 + length(aa))]
    }
    for (ii in 1:length(aa)) aa[[ii]] <- factor(aa[[ii]], exclude = NULL)
    #for (ii in 1:length(aa)) aa[[ii]] <- as.character(aa[[ii]])
    aa[['useNA']] <- useNA
    ret <- do.call("table", aa)
    if (total.margins)
        ret = atotal(ret)
    ret
}

tab.data.frame <-
function (dd, fmla, total.margins = TRUE, useNA= "ifany")
{
    if (missing(fmla))
        return(do.call("tab", as.list(dd, total.margins = total.margins,useNA=useNA)))
    xx = model.frame(fmla, dd, na.action = na.include)
    xx = c(xx, total.margins = total.margins, useNA = useNA)
    do.call("tab", xx)
}

tab.formula <-
function (fmla, dd, total.margins = TRUE, useNA = 'ifany', ...)
{
    tab(dd, fmla, total.margins = total.margins, useNA = useNA, ...)
}


