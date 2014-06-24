test_singleFunctions <-
    function()
    {
        data(isomiRs)
        checkTrue(class(degMean(DEGreportSet$deg[,4],
                DEGreportSet$counts))[[2]]=="ggplot")

    }