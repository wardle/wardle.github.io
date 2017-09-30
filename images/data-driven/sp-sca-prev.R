    sp.sca.prev <- function(sp=data.sca.prevalence(), subset, by=unique, FUN=sp.raw, long=F, ...) {
    	# apply any required subsetting
    	if(!missing(subset)) {
    		e 		<- substitute(subset)
    		subset 	<- eval(e, sp, parent.frame())
    		sp  	<- subset(sp, subset)
    	}
    	sp$unique		<- factor(1:nrow(sp))		# generate a unique ID so we can aggregate on it if no aggregation term is specified
    	# and now parse the aggregation factor
    	by		<- substitute(by)
    	# capitalise names where appropriate
    	names(sp)[2:13] <- toupper(names(sp)[2:13])

    	# special aggregate function that adds up a column and divides by total to give a proportion
    	# performing aggregation and division like this in one step means we can fix the denominator to only include studies that
    	# looked at the mutation in question
    	# you provide a callback function that gets called with the appropriate numerator and denominator
    	sp.agg.num <- function(x, total) {
    		if(all(is.na(x))==TRUE)
    			 sum.x <- NA		# if all are NAs then we return NA
    		else sum.x <- sum(x, na.rm=T)				# calculate numerator
    		if(all(is.na(total))==TRUE)
    			 sum.total <- NA	
    		else sum.total <- sum(total[!is.na(x)], na.rm=T)	# calculate the denominator from the totals that have a valid numerator ONLY
    		# if x==total, then we are looking the the TOTAL column, so we just return the RAW sum and don't call the callback function
    		if((all(x==total, na.rm=T))==TRUE & !is.na(sum.x)) 
    				return(sum(x, na.rm=T))
    		return(FUN(sum.x, sum.total, ...))
    	}

    	sp.agg.str <- function(x) {
    		paste(unique(x), collapse=", ")
    	}
    	# and now aggregate
    	sp.num <- by(sp, list(eval(by, sp, parent.frame())), FUN=function(x) { apply(x[,2:16],2, sp.agg.num, total=x$total)})
    	sp.num <- do.call("rbind", sp.num)		# and convert list to data frame
    	sp.num <- as.data.frame(sp.num)
    	sp.str <- by(sp, list(eval(by, sp, parent.frame())), FUN=function(x) { apply(x[,c(1,17:20)],2,sp.agg.str)})
    	sp.str <- do.call("rbind", sp.str)
    	sp.str <- as.data.frame(sp.str)
    	#sp = as.data.frame(cbind(sp.num, sp.str))
    	sp <- cbind(cite=sp.str[,1],sp.num,sp.str[,-1])		# note: we move the citation column back to column one so we return same data as if aggregation did not occur
	
    	# do we need to turn into long format?
    	if(long==T) {
    		sp <- reshape(sp, direction="long", timevar='mutation',
    			varying = list(names(sp[,2:13])),
    			times   = toupper(names(sp[,2:13])),
    			v.names = 'positive')
    		sp$mutation = as.factor(sp$mutation)		# convert mutation to a factor
    		sp$proportion = sp$positive / sp$total
    	}
    	return(sp)
    }