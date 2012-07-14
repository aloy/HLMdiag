#' Dot plots for influence diagnostics
#'
#' This is a function that can be used to create modified dotplots for the
#' diagnostic measures.  The plot allows the user to understand the distribution
#' of the diagnostic measure and visually identify unusual cases.
#'
#' @param data an object containing the output from \code{diagnostics()}.
#' @param type the part of the model the diagnostic corresponds to, either \code{"fixef"} or \code{"varcomp"}.
#' @param name specification of which diagnostic to plot (either COOKSD, MDFFITS, COVTRACE, COVRATIO, or rvc).
#' @param cutoff value specifying unusual values of the diagnostic
#' @param modify if \code{TRUE} will produce a space-saving modification
#' @param ... other arguments to be passed to \code{qplot()}
#' @author Adam Loy \email{aloy@@istate.edu}
#' @examples 
#' data(Oxboys, package = 'mlmRev')
#' fm <- lmer(formula = height ~ age + I(age^2) + (age + I(age^2)| Subject), data = Oxboys)
#' fmDel <- case_delete(model = fm, group = "Subject")
#' fmDiag <- diagnostics(model = fm, delete = fmDel, type = "fixef")
#' dotplot_diag(diag.out = fmDiag, type = "fixef", name = "COOKSD", cutoff = "internal", modify = TRUE, xlab = "Subject", ylab = "Cook's Distance")
#' dotplot_diag(diag.out = fmDiag, type = "fixef", name = "COOKSD", cutoff = "internal", modify = FALSE, xlab = "Subject", ylab = "Cook's Distance")
#' @export
dotplot_diag <- function(data, type = c("fixef", "varcomp"), name, cutoff = NULL, modify = FALSE, ... ){
  type <- match.arg(type)
  if(class(data) == "list"){ data <- data[[paste(type,"_diag", sep = "")]] }

  if(modify == FALSE) p <- qplot(x = reorder(IDS, get(name)), y = get(name), data = data, geom = "blank", ... )
  
  if(!is.null(cutoff)){
    if(!is.numeric(cutoff)){cutoff <- internal_cutoff(data=data, type=type, name=name)}

    if(is.numeric(cutoff)){
      if(type != "varcomp" & name != "COVRATIO"){
        data$extreme <- with(data, get(name) > cutoff)

        if(modify == TRUE){
          levels(data$IDS)[levels(data$IDS) %in% data$IDS[which(data$extreme == FALSE)]] <- "within cutoff"
          p <- qplot(x = reorder(IDS, get(name)), y = get(name), data = data, geom = "blank", ... )
        }

        if(sum(data$extreme) > 0){
          p <- p + geom_point(data = subset(data, extreme == TRUE), colour = I("red"), shape = 17) +
            geom_text(data = subset(data, extreme == TRUE), aes(label = IDS, hjust=.5, vjust=1.5, size=3))
        }

        p + geom_point(data = subset(data, extreme == FALSE), colour = I("blue")) + 
          geom_hline(aes(yintercept = cutoff), colour=I("red")) +
           # geom_text(data = subset(data, extreme == TRUE), aes(label = IDS, hjust=.5, vjust=1.5, size=3)) + 
              opts(legend.position = "none") +
                coord_flip()
      }

      else{
        data$extreme <- with(data, get(name) < cutoff[1] | get(name) > cutoff[2])

        if(modify == TRUE){
          levels(data$IDS)[levels(data$IDS) %in% data$IDS[which(data$extreme == FALSE)]] <- "within cutoff"
          p <- qplot(x = reorder(IDS, get(name)), y = get(name), data = data, geom = "blank", ... )
        }
        
        if(sum(data$extreme) > 0){
          p <- p + geom_point(data = subset(data, extreme == TRUE), colour = I("red"), shape = 17) +
            geom_text(data = subset(data, extreme == TRUE), aes(label = IDS, hjust=.5, vjust=1.5, size=3))
        }

        p + geom_point(data = subset(data, extreme == FALSE), colour = I("blue")) + 
          geom_hline(aes(yintercept = cutoff[1]), colour=I("red")) + 
            geom_hline(aes(yintercept = cutoff[2]), colour=I("red")) + 
            #  geom_text(data = subset(data, extreme == TRUE), aes(label = IDS, hjust=.5, vjust=1.5, size=3)) + 
                opts(legend.position = "none") +
                  coord_flip()	
      }
    }
  }	

  else{
    p + geom_point() + coord_flip()
  }
}


#-------------------------------------------------------------------------------------------------------

#' Calculating a cutoff value for diagnostic measures
#'
#' This function provides cutoff values using internal scaling. 
#' In other words, a measure of relative standing is used (3*IQR)
#' to specify unusual values of the diagnostic of interest relative
#' to the vector of diagnostics.
#'
#' @param data an object containing the output from \code{diagnostics()}
#' @param type the part of the model the diagnostic corresponds to, either \code{"fixef"} or \code{"varcomp"}.
#' @param name specification of which diagnostic to plot (either COOKSD, MDFFITS, COVTRACE, COVRATIO, or rvc).
#' @author Adam Loy \email{aloy@@istate.edu}
internal_cutoff <- function(data, type, name){
  series <- data[, name]
  q3 <- quantile(series, p=0.75)
  series.iqr <- IQR(series)
	
  if(name == "COVRATIO" | type == "varcomp"){
    q1 <- quantile(series, p=0.25)
    cutoff <- c(lower = q1 - 3*series.iqr, upper = q3 + 3*series.iqr)	
  }
	
  else{cutoff <- q3 + 3*series.iqr}
	
  return(cutoff)	
}
