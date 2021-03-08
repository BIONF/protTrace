args <- commandArgs(trailingOnly = TRUE)

        path = args[1]

        y <- scan(path,quiet = TRUE)
        x = seq(0.1, length(y)*0.1, 0.1)
        y=1-y
        
        myfunc = function(t, params){ (1*params[1]*exp(params[2]*t))/(1+params[1]*(exp(params[2]*t)-1)) }
        
        fit <- function(givenData,givenStartingN)
        {
          return(tryCatch(
          {
              nls(formula=givenData~(1*N*exp(r*x))/(1+N*(exp(r*x)-1)),start = c(N=givenStartingN,r=0.6),control = list(maxiter=500))
          }, error=function(e) NULL))
        }
        
        loopingFitFunc <- function(givenDataset)
        {
          finalCoef = c(N=0.01,r=0.6)
          nseq = 10^(-2:-8)
          counter = 0
          while (counter <= length(nseq)) 
          {
            intermediaryFitModel = fit(givenDataset,nseq[counter])
            if (!is.null(intermediaryFitModel)) {
              intermediaryCoef = coef(intermediaryFitModel)
              if (intermediaryCoef[1] < 0 || intermediaryCoef[2] < 0) {
                intermediaryFitModel = NULL
              } else
              {
                finalCoef['N']=intermediaryCoef[1]
                finalCoef['r']=intermediaryCoef[2]
                break
              }
            }
            counter = counter + 1
          }
          return(finalCoef)
        }
        
        par = loopingFitFunc(y)
        y2 = myfunc(x, par)
        
        #png(paste(path, ".png", sep=""), width=600, height=600)
        #par(mar=c(5,5,2,2))
        #plot(x, 1-y2, type="l", ylab="Fraction of detection", xlab="Simulation distance (in substitutions per site)", cex.lab=2, cex.axis=2)
        #points(x, 1-y)
        #dev.off()

        #pdf(paste(path, ".pdf", sep=""), width=20, height=20)
        #par(mar=c(5,5,2,2))
        #plot(x, 1-y2, type="l", ylab="Fraction of detection", xlab="Simulation distance (in substitutions per site)", cex.lab=2, cex.axis=2)
        #points(x, 1-y)
        #invisible(dev.off())
        
        write(par[1], file=paste(path, "_parameter", sep=""), append = FALSE)
        write(par[2], file=paste(path, "_parameter", sep=""), append = TRUE)

        #residual_sumofsquares = sum(resid(fitModel)^2);
        #write(residual_sumofsquares, file=paste(path, "_residualsumofsquares", sep=""), append = FALSE)

q()
