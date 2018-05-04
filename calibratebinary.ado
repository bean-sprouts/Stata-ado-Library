/*********************************************/
/*THIS PROGRAM PRODUCES THE CALIBRATION PLOTS*/
/*********************************************/
/*
INPUTS:
REQUIRED: TWO VARIABLE NAMES, THE FIRST IS THE BINARY AND THE SECOND IS THE PREDICTED PROBABILITY
OPTIONS: THE NUMBER OF KERNELS DESIRD IN PLOT, DEFAULT IS 10
*/

capture program drop calibratebinary
program calibratebinary
	version 11.0
	syntax varlist(min=2 max=2 numeric) [in] [if] [pweight] [,KERnel(integer 10) GENerate(name) ///
			lowessopts(string) lineopts(string) rcapopts(string) scatteropts(string) functionopts(string) scale100 * ] 
	
	marksample touse
	tokenize `varlist'
	local outcome `1'
	tempvar phat
	qui g `phat'=`2'
	if "`scale100'"!="" replace `phat'=`phat'/100
	
	/*CHECKING THAT THE OUTCOME IS NUMERIC 0/1 BINARY VARIABLE*/
	quietly levelsof `outcome'
	if "`r(levels)'"!="0 1" {
		di as error "Outcome must be numeric zero-one variable"
		exit 198 
	}
	
	
	sort `phat'
	tempvar group
	qui egen `group' = cut(`phat') if `touse' & !mi(`phat') & !mi(`outcome'), group(`kernel') 
	quietly replace `group' = `group' + 1
	tempvar phatmean observemean observelb observeub
	
	qui g `phatmean'=.
	qui g `observemean'=.
	qui g `observelb'=.
	qui g `observeub'=.
	
	qui forvalues i=1(1)`kernel' {
		*calculating mean rate of the observed outcome within each group
		tempname  meanmatrix
		mean `outcome' [`weight' `exp'] if `group'==`i'

			matrix `meanmatrix'=r(table)
			replace `observemean'=`meanmatrix'[1,1] if `group'==`i'
			
			if "`weight'"=="" { /*calculating confidence interval for non-weighted analysis using the binomial distribution*/
				ci `outcome' if `group'==`i', exact binomial
				replace `observelb'=`r(lb)' if `group'==`i'
				replace `observeub'=`r(ub)' if `group'==`i'
			}
			else { /*use Xbar mean and SD for confidence interval if pweights present*/
				replace `observelb'=`meanmatrix'[5,1] if `group'==`i'
				replace `observeub'=`meanmatrix'[6,1] if `group'==`i'
			}
			
		*calculating mean rate of the predicting varaible (phat) within each group
		tempname  meanmatrix
		mean `phat' [`weight' `exp'] if `group'==`i'
			matrix `meanmatrix'=r(table)
			replace `phatmean'=`meanmatrix'[1,1] if `group'==`i'
	}
			
	* get a non-parametric smoothed estimate of the actual results (only for non-weighted analyses)
	tempvar lwsobserve
	if "`weight'"=="" {
	qui lowess `outcome' `phat' if `touse', gen(`lwsobserve') nograph logit `lowessopts'
		qui replace `lwsobserve'=invlogit(`lwsobserve')
	}
	else g `lwsobserve'=.
	
	* calibration plot
	twoway (line `lwsobserve' `phat', sort clcolor(gs0) clpat(dash) `lineopts') || ///
		   (rcap `observelb' `observeub' `phatmean', lcolor(gs8) `rcapopts') || ///
		   (scatter `observemean' `phatmean', mcolor(gs0) msize(1.5) `scatteropts') || ///
		   (function y=x, range(0 1) clcolor(gs0) clpat(solid) `functionopts'), ///	
			title("Calibration Plot") ytitle("Actual", margin(medium)) xtitle("Predicted", margin(medium))  ///
			legend(off) `options'

	if "`generate'"!="" {
		foreach v in group phatmean observemean observelb observeub lwsobserve {
			qui g `generate'`v'=``v''
		}
	}
end


