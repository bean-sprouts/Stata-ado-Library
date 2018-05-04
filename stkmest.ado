/*
PROGRAM: stkmest.ado
PROGRAMMER: Daniel
DATE: 4/30/2013
DESCRIPTION: Calculates and returns Kaplan Meier estiamtes
*/


capture program drop stkmest
program stkmest, rclass
	version 12
	syntax [if] [in], time(numlist>0) [Failure DECimal(integer 0) by(varlist) table saving(string asis)]
	preserve
	
	*keeping needed observations
	qui capture keep `if' `in'
	
	
	tempfile kmest
	
	*calculaing KM estimates
	if trim("`by'")=="" qui sts list, at(0 `time') saving(`kmest') `failure'
	else qui sts list, by(`by') at(0 `time') saving(`kmest') `failure'
	*loading the estiamtes into memory
	use `kmest', clear
	/*zero included by default due to bug with sts list if only one timepoint is selected*/
	qui drop if time==0
	
	*defining central estiamte varname
	if trim("`failure'")=="" local kmest="survivor"
	else local kmest="failure"
	
	*formatting results
	qui g kmest=trim(string(`kmest'*100,"%9.`decimal'f"))+"%"
	qui g ci=trim(string(lb*100,"%9.`decimal'f"))+"%, "+trim(string(ub*100,"%9.`decimal'f"))+"%" if !mi(lb) & !mi(ub)
	qui replace ci="NA, NA" if mi(ci)
	
	qui g Kaplan_Meier_estimate=kmest+" (95% CI "+ci+")"
	
	*listing results into command window
	if trim("`table'")=="" list `by' time Kaplan_Meier_estimate, clean noobs
	*outputting results in table format
	else listtex `by' time kmest ci, headlines(`=subinstr("`=trim("`by' time kmest ci")'"," ","&",.)') type
	
	*saving estimates into dataset if requested
	if trim("`saving'")!="" save `saving'
	
	*returning estimates as locals
	qui count
	return local N=`r(N)'
	foreach i of numlist `r(N)'/1 {
		return local full_`i'=Kaplan_Meier_estimate[`i']
		return local cifmt_`i'=ci[`i']
		return local kmestfmt_`i'=kmest[`i']
		return local ub_`i'=ub[`i']
		return local lb_`i'=lb[`i']
		return local kmest_`i'=`kmest'[`i']
		return local time_`i'=time[`i']
		if trim("`by'")!="" {
			foreach v of varlist `by' {
					return local `v'_`i'=`v'[`i']
			}
		}
	}
	
	restore
end
