#' Arguments (included in plsDA function within DiscriMiner package)
#' @aliases isoPLSDA
#' @usage isoPLSDA(obj, var, validation = NULL, learn = NULL, test = NULL,
#'  tol = 0.001, nperm = 400, refinment = FALSE, vip = 1.2)
#' @param obj IsomirDataSeq
#' @param var column name in design data.frame
#' @param validation type of validation, either NULL or "learntest". Default NULL
#' @param learn	optional vector of indices for a learn-set. Only used when validation="learntest". Default NULL
#' @param test	optional vector of indices for a test-set. Only used when validation="learntest". Default NULL
#' @param tol tolerance value based on maximum change of cumulative R-squared coefficient for each additional PLS component. Default tol=0.001
#' @param nperm	number of permutations to compute the PLD-DA p-value based on R2 magnitude. Default nperm=400
#' @param refinment logical indicating whether a refined model, based on filtering out variables with low VIP values
#' @param vip	Variance Importance in Projection threshole value when a refinement precess is considered. Default vip=1.2
#' @return PLS model
#' 
#' @export
isoPLSDA <- function(obj,var ,validation = NULL, learn = NULL, test = NULL, 
                     tol = 0.001, nperm = 400, refinment = FALSE, vip = 1.2)
{
    tryCatch ({
        class(obj@normcounts)
    },error = function(e){
        return("please, run first normIso.")
    })
    variables <- t(obj@normcounts)
    group <- obj@design[,var]
    if (length(group)<6){
        return("this analysis only runs with group larger than 6.")
    }
    # Auxiliar data containing variable names and numeric ID per variable
    dataVariables <- data.frame( variable = colnames(variables), 
                                 id = c(1:dim(variables)[2]))
    
    # PLS-DA model
    model.plsDA <- plsDA(variables, group, validation, learn, test)
    
    # R-squared coefficents (number of components are selected based on 
    # cumulative R-squared coefficient changes)
    if(model.plsDA$R2[dim(model.plsDA$R2)[1],3] < tol){
        a <- which(model.plsDA$R2[,3] < tol)
        if(a[1] > 1){ R2.mat <- model.plsDA$R2[1:(a[1]-1),] }
        if(a[1] == 1){ R2.mat <- model.plsDA$R2[1,] }
    }
    if(model.plsDA$R2[dim(model.plsDA$R2)[1],3] >= tol){
        R2.mat <- model.plsDA$R2
    }
    
    # model-based p-value (permutation analysis)
    R2.perm <- R2PermutationVector(variables, group, validation, learn, 
                                   test, tol, nperm)
    p.val <- sum(R2.mat[dim(R2.mat)[1],4] <= R2.perm, na.rm=TRUE)/nperm
    
    if(refinment == TRUE){
        
        # refine model based on VIPs
        vip.max <- apply(model.plsDA$VIP[,1:dim(R2.mat)[1]], 1, max)	
        dataVIP <- data.frame(variable = row.names(model.plsDA$VIP), 
                              VIP = vip.max)
        dataVIPref <- dataVIP[dataVIP$VIP >= vip,]
        sel.vars <- merge(dataVariables , dataVIPref, 
                          by.x="variable", by.y="variable")
        variables.ref <- variables[,sel.vars[,2]]
        
        # PLS-DA model with contributing variables
        model.plsDA.ref <- plsDA(variables.ref, group, validation, learn, test)
        if(model.plsDA.ref$R2[dim(model.plsDA.ref$R2)[1],3] < tol){
            a.ref <- which(model.plsDA.ref$R2[,3] < tol)
            if(a.ref[1] > 1){ R2.mat.ref <- model.plsDA.ref$R2[1:(a.ref[1]-1),]}
            if(a.ref[1] == 1){ R2.mat.ref <- model.plsDA.ref$R2[1,] }
        }
        if(model.plsDA.ref$R2[dim(model.plsDA.ref$R2)[1],3] >= tol){
            R2.mat.ref <- model.plsDA.ref$R2
        }
        
        # p-value using fixed vip variables while shuffling individuals
        R2.perm.ref1 <- R2PermutationVector(variables.ref, group, 
                                            validation, learn, test, tol, nperm)
        p.val.ref1 <- sum(R2.mat.ref[dim(R2.mat.ref)[1],4] <= R2.perm.ref1,
                          na.rm=TRUE)/nperm
        
        # p-value based on performing the whole refinement process 
        # in  each permutation iteration
        R2.perm.ref2 <- R2RefinedPermutationVector(variables, group, 
            validation, learn, test, tol, nperm, vip)
        p.val.ref2 <- sum(R2.mat.ref[dim(R2.mat.ref)[1],4] <= R2.perm.ref2,
                          na.rm=TRUE)/nperm
        
        res <- list(R2.mat, model.plsDA$components[,1:dim(R2.mat)[1]], 
                    p.val, R2.perm, R2.mat.ref, 
                    model.plsDA.ref$components[,1:dim(R2.mat.ref)[2]], 
                    p.val.ref2, R2.perm.ref2, p.val.ref1, R2.perm.ref1)
        names(res) <- c("R2Matrix", "components", 
                        "p.val", "R2PermutationVector", "R2RefinedMatrix", 
                        "componentsRefinedModel", 
                        "p.valRefined", "R2RefinedPermutationVector",
                        "p.valRefinedFixed", "R2RefinedFixedPermutationVector")
        print(paste0("pval:",res$p.val))
        return(res)
        
        
    }
    
    if(refinment == FALSE){
        res <- list(R2.mat, model.plsDA$components[,1:dim(R2.mat)[1]],
                    p.val, R2.perm)
        names(res) <- c("R2Matrix", "components", "p.val",
                        "R2PermutationVector")
        print(paste0("pval:",res$p.val))
        return(res)
    }
}


#' Function that computes p-values when only individuals 
#' are considered to be permuted
#' @param variables   matrix or data frame with explanatory variables
#' @param group   	vector or factor with group memberships
#' @param validation	type of validation, either NULL or "learntest". 
#' Default NULL
#' @param learn		optional vector of indices for a learn-set. 
#' Only used when validation="learntest". Default NULL
#' @param test		optional vector of indices for a test-set. 
#' Only used when validation="learntest". Default NULL
#' @param tol	 tolerance value based on maximum change of cumulative 
#' R-squared coefficient for each additional PLS component. Default tol=0.001
#' @param nperm	 number of permutations to compute the PLD-DA p-value 
#' based on R2 magnitude. Default nperm=400
#' @return PLS model
#' 

R2PermutationVector <- function(variables, group, validation,
                                learn, test, tol, nperm)
{
	# intitialize R2 vector
	R2perm <- rep(NA,nperm)

	# number of individuals
	n <- length(group)

	# permutation loop
	for(p in 1:nperm){
	
		# shuffle individuals
		mysample <- group[sample(1:length(group), n, replace=FALSE)]
		
		# pls-da with groups shuffled
		model.perm <- plsDA(variables, mysample, validation, learn, test)

		# select components based on R2 contribution 
		if(model.perm$R2[dim(model.perm$R2)[1],3] <= tol){
			a.perm <- which(model.perm$R2[,3] < tol)
			if(a.perm[1] > 1){ R2perm[p] <- model.perm$R2[a.perm[1]-1,4] }
			if(a.perm[1] == 1){ R2perm[p] <- model.perm$R2[1,4] }
		}
		if(model.perm$R2[dim(model.perm$R2)[1],3] > tol){
			R2perm[p] <- model.perm$R2[dim(model.perm$R2)[1],4]
		}

	}

	return(R2perm)

}


#' Function that computes p-values when a refinment process is considered 
#' per each permutation iteration
#' Arguments (included in plsDA function within DiscriMiner package)
#' @param variables   matrix or data frame with explanatory variables
#' @param group   	vector or factor with group memberships
#' @param validation	type of validation, either NULL or "learntest". 
#' Default NULL
#' @param learn		optional vector of indices for a learn-set. 
#' Only used when validation="learntest". Default NULL
#' @param test		optional vector of indices for a test-set. 
#' Only used when validation="learntest". Default NULL
#' @param tol	 tolerance value based on maximum change of cumulative 
#' R-squared coefficient for each additional PLS component. Default tol=0.001
#' @param nperm	 number of permutations to compute the PLD-DA p-value 
#' based on R2 magnitude. Default nperm=400
#' @param vip	 Variance Importance in Projection threshole value 
#' when a refinement precess is considered. Default vip=1.2
#' @return PLS model
#' 
R2RefinedPermutationVector <- function(variables, group, validation, learn,
                                       test, tol, nperm, vip)
{
    
    # Auxiliar data containing variable names and numeric ID per variable
    dataVariables <- data.frame( variable = colnames(variables),
                                 id = c(1:dim(variables)[2]))
    
    # number of individuals
    n <- length(group)
    
    # intitialize R2 vector
    R2Refined.perm <- rep(NA,nperm)
    
    # permutation loop
    for(p in 1:nperm){
        
        # shuffle individuals
        mysample <- group[sample(1:length(group), n, replace=FALSE)]
        # pls-da with groups shuffled
        model.perm <- plsDA(variables, mysample, validation, learn, test)
        
        if(model.perm$R2[dim(model.perm$R2)[1],3] < tol){
            a.perm <- which(model.perm$R2[,3] < tol)
            if(a.perm[1] > 1){ v <- (a.perm[1]-1) }
            if(a.perm[1] == 1){ v <- 1 }
        }
        if(model.perm$R2[dim(model.perm$R2)[1],3] >= tol){
            v <- dim(model.perm$R2)[1]
        }
        
        # refine model based on VIPs
        vip.perm.max <- apply(model.perm$VIP[,1:v], 1, max)	
        dataVIP.perm <- data.frame(variable = row.names(model.perm$VIP),
                                   VIP = vip.perm.max)
        dataVIP.perm.ref <- dataVIP.perm[dataVIP.perm$VIP >= vip,]
        sel.vars.perm <- merge(dataVariables , dataVIP.perm.ref,
                               by.x="variable", by.y="variable")
        variables.perm.ref <- variables[,sel.vars.perm[,2]]
        
        if(!is.null(dim(variables.perm.ref)[2])){
            model.perm.ref <- plsDA(variables.perm.ref, group, 
                autosel = TRUE, comps = 2, validation, learn, test)
            
            # select components based on R2 contribution 
            if(model.perm.ref$R2[dim(model.perm.ref$R2)[1],3] <= tol){
                a.perm.ref <- which(model.perm.ref$R2[,3] < tol)
                if(a.perm.ref[1] > 1){ R2Refined.perm[p] <- model.perm.ref$R2[a.perm.ref[1]-1,4] }
                if(a.perm.ref[1] == 1){ R2Refined.perm[p] <- model.perm.ref$R2[1,4] }
            }
            if(model.perm.ref$R2[dim(model.perm.ref$R2)[1],3] > tol){
                R2Refined.perm[p] <- model.perm.ref$R2[dim(model.perm.ref$R2)[1],4]
            }
        }
    }
    
    return(R2Refined.perm)  
}

#' Components plot (LATTICE PLOT)
#' @aliases isoPLSDAplot
#' @usage isoPLSDAplot(components, groups)
#' @param components PLS-DA components as it comes from isoPLSDA main function
#' @param groups	vector or factor with group memberships
#' @return plot
#' 
#' @export
isoPLSDAplot <- function(components, groups)
{  
    datacomponents <- data.frame(condition = groups, components)
    t <- dim(datacomponents)[2]-1
    n <- length(levels(factor(groups)))
    super.sym <- trellis.par.get("superpose.symbol")
    splom(~datacomponents[,2:dim(datacomponents)[2]], groups = condition, data = datacomponents,
          panel = panel.superpose,
          key = list(title = "PLS-DA components",
                     columns = t, 
                     points = list(pch = super.sym$pch[1:n],
                                   col = super.sym$col[1:n]),
                     text = list(levels(factor(groups)))))
}

