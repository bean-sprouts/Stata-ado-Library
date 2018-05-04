/*
PROGRAM: skkmdiff
PROGRAMMER: Daniel
DESCRIPTION:  This ado file estimates confidence intervals about differences in survival estimates.
				See "Survival Analysis: Techniques for Censored and Truncated Data" Second Edition By John P Klein and Melvin Moeschberger
				Section 7.8 "Test Based on Differences in Outcome at a Fixed Point in Time" page 234 for details about the methodology.
*/

capture program drop stkmdiff
program stkmdiff, rclass
	version 12
	syntax [if] [in], time(numlist) by(varname numeric) [DECimal(integer 0) DIFFDECimal(integer 1) Failure table]
	preserve
	*keeping selected observations
	capture keep `if' `in'
	
	*checking the coding of the by variable
	qui levelsof `by'
	capture assert "`r(levels)'"=="0 1"
	if _rc>0 {
		disp as err "By variable must by coded as 0 and 1, and have observed values for both groups."
		exit 198
	}
	
	*estimating survival prob and se
	tempfile kmest
	qui sts list, at(0 `time') saving(`kmest') by(`by') `failure'
	use `kmest', clear
	qui drop if time==0
	
	*local with name of prob estiamte varaible name
	if trim("`failure'")=="" local probestimate="survivor"
	else local probestimate="failure"
	
	*reshaping to be one line per timepoint
	keep `probestimate' std_err time `by'
	qui reshape wide `probestimate' std_err, i(time) j(`by')

	*calculating CI using Greenwoods se estimates (assumes normality on the 0-1 interval) 
	g probdiff=`probestimate'0 - `probestimate'1
	g probdiffse=sqrt(std_err0^2+std_err1^2)
	g probdiff_lb=probdiff - invnormal(0.975)*probdiffse
	g probdiff_ub=probdiff + invnormal(0.975)*probdiffse

	*formattting results
	g `probestimate'0_fmt=string(`probestimate'0*100,"%9.`decimal'f")+"%"
	g `probestimate'1_fmt=string(`probestimate'1*100,"%9.`decimal'f")+"%"
	g probdiff_fmt=string(probdiff*100,"%9.`diffdecimal'f")+"%"
	g probci_fmt=string(probdiff_lb*100,"%9.`diffdecimal'f")+"%, "+string(probdiff_ub*100,"%9.`diffdecimal'f")+"%"
	
	*return estimates as locals
	return local km_0=`probestimate'0
	return local km_1=`probestimate'1
	return local diff=probdiff
	return local diff_lb=probdiff_lb
	return local diff_ub=probdiff_ub
	return local diffse=probdiffse
	
	
	rename `probestimate'0_fmt `by'_0
	rename `probestimate'1_fmt `by'_1
	rename probdiff_fmt `probestimate'_difference 
	rename probci_fmt `probestimate'_difference_ci
	
	*displaying results
	if trim("`table'")=="" list time `by'_0 `by'_1 `probestimate'_difference  `probestimate'_difference_ci, clean noobs ab(21)
	else listtex time `by'_0 `by'_1 `probestimate'_difference `probestimate'_difference_ci, type headlines(`=subinstr("time `by'_0 `by'_1 `probestimate'_difference `probestimate'_difference_ci"," ","&",.)')
	
	restore
end


