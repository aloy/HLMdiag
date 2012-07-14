#' Decomposing formula for an LMM or HLM
.parse_formula <- function(form, type="fixed", incl=TRUE) {
# returns a data frame of the form
#  term   		character of variable name
#  type   		either random, fixed, or group
#  inclusion 	binary: include term?
	if (length(form) == 1) return(c(as.character(form), type, incl))
	
	if((form[[1]] == 'I') | (form[[1]] == '*') | (form[[1]] == ':')) 
		return(c(deparse(form), type, incl))
	
	if(form[[1]] == '~') return(.parse_formula(form[[3]]))
	
	#if(form[[1]] %in% getGroupMembers(Math)) return(c(deparse(form), type, incl))
	
	if ((form[[1]] == '+') | (form[[1]] == '-')) { 
		incl <- form[[1]] == '+'

		resright <- .parse_formula(form[[3]], type=type, incl=incl)
		resleft <- .parse_formula(form[[2]], type=type, incl=incl)
		
		return(rbind(resright, resleft))
	}
	
	if (form[[1]] == "(") {
		return(.parse_formula(form[[2]], type="random"))
	}

	if (form[[1]] == "|") {
		resright <- .parse_formula(form[[3]], type="group", incl=TRUE)
		resleft <- .parse_formula(form[[2]], type="random", incl=TRUE)

		return(rbind(resright, resleft))		
	}
}

#' Formatting a decomposed formula for an LMM or HLM
.parse_f <- function(formula) {
	res <- data.frame(.parse_formula(formula), row.names=NULL)
	names(res) <- c("name","type","included")
	res
}
