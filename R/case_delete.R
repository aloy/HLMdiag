#' Implementing Case Deletion.
#'
#' This function is used to iteratively delete groups corresponding to the levels of a 
#' two-level hierarchical linear model. It uses \code{lmer()} to fit the models for each
#' deleted case (i.e. uses brute force). To investigate numerous levels of
#' the model, the function will need to be called multiple times, specifying
#' the group (level) of interest each time.
#'
#' @param model the original hierarchical model fit using \code{lmer()}
#' @param group a variable used to define the group for which cases will be deleted.
#'   If this is left \code{FALSE}, then the function will delete individual observations.
#' @param type the part of the model for which you are obtaining deletion diagnostics:
#'   the fixed effects (\code{fixef}), variance components (\code{varcomp}), or \code{both}
#' @return a list with the following compontents: 
#'   \item{fixef.original}{the original fixed effects}
#'   \item{ranef.original}{the origingal random effects}
#'   \item{vcov.original}{the original variance-covariance parameters}
#'   \item{varcomp.original}{the original variance components}
#'   \item{fixef.delete}{a list of the fixed effects obtained through case deletion}
#'   \item{ranef.delete}{a list of the random effects obtained through case deletion}
#'   \item{vcov.delete}{a list of the variance-covariance parameters obtained 
#'   through case deletion}
#'   \item{fitted.delete}{a list of the fitted values obtained through case deletion}
#'   \item{varcomp.delete}{a list of the variance components obtained through case
#'   deletion}
#' @author Adam Loy \email{aloy@@istate.edu}
#' @references Christensen, R., Pearson, L.M., and Johnson, W. (1992),
#' ``Case-Deletion Diagnostics for Mixed Models,'' \emph{Technometrics},
#' 34, 38 -- 45.
#'
#' Schabenberger, O. (2004),``Mixed Model Influence Diagnostics,''
#' in \emph{Proceedings of the Twenty-Ninth SAS Users Group International Conference},
#' SAS Users Group International.
#' @examples
#' data(Oxboys, package = 'mlmRev')
#' fm <- lmer(formula = height ~ age + I(age^2) + (age + I(age^2)| Subject), data = Oxboys)
#' fmDel <- case_delete(model = fm, group = TRUE, type = "both")
#'
#' \dontrun{library(mlmRev)
#' exm1 <- lmer(normexam ~ standLRT + sex + schgend + (1 | school), data = Exam)
#' exm1DEL <- case_delete(model = exm1, group = TRUE, type = "both")}
#' @export
case_delete <- function(model, group = FALSE, type = c("both", "fixef", "varcomp")){
  fixef.delete <- NULL
  vcov.delete <- NULL
  varcomp.delete <- NULL
  ranef.delete <- NULL
  fitted.delete <- NULL
  
  type <- match.arg(type) #default is "both"
  if(group == FALSE){ # SINGLE CASE DELETION DIAGNOSTICS
    model@frame$INDEX <- row(model@frame)[,1]
#
#    if(type %in% c("both", "fixef")){
#      fixef.delete <- NULL
#      vcov.delete <- NULL
#    }
#
#    if(type %in% c("both", "varcomp")){
#      varcomp.delete <- NULL
#    }

#    ranef.delete <- NULL
#    fitted.delete <- NULL

    for(i in model@frame$INDEX){
      model.delete <- lmer(formula = formula(model), data = subset(model@frame, INDEX != i))

      if(type %in% c("both", "varcomp")){
        ranef.delete[[i]] <- data.frame(deleted = i, id = rownames(ranef(model.delete)[[1]]), ranef(model.delete)[[1]])
        sig <- attr(VarCorr(model.delete), "sc")
        varcomp.delete[[i]] <- c(sigma2 = sig^2, diag(VarCorr(model.delete)[[names(model@flist)]]))
      }

      if(type %in% c("both", "fixef")){
        fixef.delete[[i]] <- fixef(model.delete)
        vcov.delete[[i]] <- as.matrix(vcov(model.delete))
      }

      fitted.delete[[i]] <- data.frame(deleted = i, model.delete@frame, fitted(model.delete))
    }
  }

  else{ # MULTIPLE CASE DELETION DIAGNOSTICS
    column.groups <- which(names(model@frame) == names(model@flist))
    data.delete <- split(model@frame, model@frame[,names(model@flist)])
    data.delete <- lapply(data.delete, function(df){
      data.delete[[unique(df[,names(model@flist)])]] <- NULL
      do.call('rbind', data.delete)
    })

    model.delete <- lapply(data.delete, lmer, formula = formula(model))

    if(type %in% c("both", "varcomp")){
      ranef.delete <- lapply(model.delete, function(x){
        data.frame(deleted = setdiff(model@frame[, names(model@flist)], x@frame[, names(model@flist)]),
                   id = rownames(ranef(x)[[1]]), ranef(x)[[1]])
      })
      varcomp.delete <- lapply(model.delete, function(x){
        sig <- attr(VarCorr(x), "sc")
        varcomp <- c(sigma2 = sig^2, diag(VarCorr(x)[[names(model@flist)]]))
        return(varcomp)
      })
    }

    if(type %in% c("both", "fixef")){
      fixef.delete <- lapply(model.delete, fixef)
    
      vcov.delete <- lapply(model.delete, vcov)
      vcov.delete <- lapply(vcov.delete, as.matrix)
    }


    fitted.delete <- lapply(model.delete, function(x){
      data.frame(deleted = setdiff(model@frame[, column.groups], x@frame[, column.groups]),
                 x@frame, fitted(x))
    })
  }

  # Organizing results
  if(type %in% c("both", "fixef")){
    fitted.delete <- do.call('rbind', fitted.delete)
    fixef.delete <- do.call('rbind', fixef.delete)
  }

  if(type %in% c("both", "varcomp")){
    ranef.delete <- do.call('rbind', ranef.delete)
  }

  fixef.original <- model@fixef
  ranef.original <- ranef(model)[[names(model@flist)]]
  vcov.original <- as.matrix(vcov(model))

  sigma.original <- attr(VarCorr(model), "sc")
  varcomp.original <- c(sigma2 = sigma.original^2, diag(VarCorr(model)[[names(model@flist)]]))

  val <- list(fixef.original = fixef.original, ranef.original = ranef.original,
              vcov.original = vcov.original, varcomp.original = varcomp.original,
              fixef.delete = fixef.delete, ranef.delete = ranef.delete,
              vcov.delete = vcov.delete, fitted.delete = fitted.delete,
              varcomp.delete = varcomp.delete)

  attr(val, "type") <- type

  return(val)
}
