
// October 12, 2014: Program updated by Emily Vertosick. "nobs(numlist)" option added.
	// The "nobs" option allows you to specify how many estimates are calculated.
	// The total number of observations will be divided by the number specified in "nobs" to get "n"
	// and the estimates will be calculated and saved for every nth observation.

/*
THIS PROGRAM IS IN PROGRESS
NONPARAMETRIC SURVIVAL SMOOTHING FUNCTION
*/

capture program drop stlowess
program stlowess

	version 12
	syntax varlist(numeric max=1) [if] [in] , time(numlist) [nobs(numlist) bwidth(real .80) GENerate(name) Failure nograph *]
	
	*marking observations to use 
	marksample touse
	
	*not using observations that are not st set OR are missing the covariate of interest
	qui replace `touse'=0 if _st==0 | mi(`varlist')
		
	* bandwidth must be between 0 and 1
	if !inrange(`bwidth',0,1) {
		disp as error "Bandwidth must be between 0 and 1."
		exit
	}
	
	*sorting data with utilized variables on top, then by the varlist
	gsort -`touse' `varlist'
	
	*calculating the number of observations to skip when estiamting risk (default is use all variables)
	qui count if `touse'==1
		local touseN=`r(N)'
	if "`nobs'"=="" local byN=1
	else local byN=max(floor(`touseN'/`nobs'),1)
	
	*defining which obs to include based on bandwidth
	local k=floor( (`touseN'*`bwidth' - 0.5)/2 )

	*saving variables to reset stset for each weighted analysis.
	foreach v of varlist _st _d _t _t0 {
		tempvar `v'
		qui g ``v''=`v'
	}

	*initializing lowess estimates variable
	tempvar lwsest
	qui g `lwsest'=.
	
	// for each tagged observation (every observation tagged by default)
	local i = 1 - `byN'
	qui while `i' < `touseN' - `byN' {
		*increasing by specifed increment of analysis
		local i=`i' + `byN'
						
		* calculating the covariate value for the i-th observation
		local `varlist'`i'=`varlist'[`i']
		* calculating the lower limit of covariate values to include in analysis
		local i_minus=max(1,`i' - `k')
		* calculating the upper limit of covariate values to include in analysis
		local i_plus=min(`touseN',`i' + `k')
		* calculating delta for weighted KM analysis (see LOWESS documentation)
		local delta=1.0001*max(`varlist'[`i_plus'] - `varlist'[`i'], `varlist'[`i'] - `varlist'[`i_minus'])

		* generating tricube weighted values
		tempvar weight
		qui g `weight'=( 1 - (abs(`varlist' - `varlist'[`i'])/`delta')^3 )^3 in `i_minus'/`i_plus'
		
		* setting survival data with weights
		qui stset `_t' [pw=`weight'], f(`_d') origin(`_t0')

		* calculating survival estimates
		kmest `time', local(km) `failure'	
		replace `lwsest'=`km' in `i'
	}	

	*resetting previous  stset
	qui stset `_t', f(`_d') origin(`_t0')
	
	* printing graph
	if trim("`graph'")=="" {
		twoway line `lwsest' `varlist', ylabel(0(0.2)1, format("%9.1f")) ///
				ytitle("Locally-Weighted Kaplan-Meier Estimate" "(Analysis Time=`time')") ///
				note("Bandwidth=`bwidth'") `options'
	}
	

	*saving out estimate if requested
	qui if trim("`generate'")!="" rename `lwsest' `generate'
end



capture program drop kmest
program kmest, rclass
	syntax anything, local(string) [Failure]
	preserve
	tempfile survlist
	qui sts list, at(0 `anything') saving(`survlist') `failure'
	
	use `survlist', clear
	drop if time==0
	list
	
	if trim("`failure'")=="" c_local `local'=survivor
	else  c_local `local'=failure
	
	restore
end

