#' Generate mboost prediction
#'
#' @param fit fitted mboost model
#' @param data training data
#' @param newdata test data
#' @param bl base learner
#' @param link link function
#' @param cols columns to keep in the output
#'
#' @returns A data.table with predictions and other relevant information
#' @export

predict_result <- function(fit, data, newdata, bl, link = "logit",
							cols = c("ID", "Status", "FamilyID", "Clinical_Findings")){
  # prediction of each base leaner
  prd <- predict(fit, type = "link", newdata = newdata, which = bl)

  rsp <- predict(fit, type = "response", newdata = newdata)

  # logit link function
  if(is.character(link)){
  	link_fun <- make.link(link)$linkinv
  } else{
  	link_fun <- link
  }

  return(data.table(data[, ..cols],
                    offset = attr(prd, "offset"),
                    prd,
                    eta = attr(prd, "offset") + rowSums(prd),
                    response = as.vector(rsp)))
}

#' Title
#'
#' @param data
#' @param measurevar
#' @param groupvars
#' @param na.rm
#' @param conf.interval
#' @param .drop
#'
#' @returns
#' @export
#'
#' @examples
summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE,
                      conf.interval=.95, .drop=TRUE) {

    # New version of length which can handle NA's: if na.rm==T, don't count them
    length2 <- function (x, na.rm=FALSE) {
        if (na.rm) sum(!is.na(x))
        else       length(x)
    }

    # This does the summary. For each group's data frame, return a vector with
    # N, mean, and sd
    datac <- plyr::ddply(data, groupvars, .drop=.drop,
                         .fun = function(xx, col) {
                             c(N    = length2(xx[[col]], na.rm=na.rm),
                               mean = mean   (xx[[col]], na.rm=na.rm),
                               sd   = sd     (xx[[col]], na.rm=na.rm)
                               )
                             },
                         measurevar
                         )

    # Rename the "mean" column
    datac <- rename(datac, c("mean" = measurevar))

    datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean

    # Confidence interval multiplier for standard error
    # Calculate t-statistic for confidence interval:
    # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
    ciMult <- qt(conf.interval/2 + .5, datac$N-1)
    datac$ci <- datac$se * ciMult

    return(datac)
}
