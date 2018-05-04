/*
	PROGRAM: auctab.ado
	PROGRAMMER: Daniel
	DESCROPTION:  Formats and outputs AUC along with confidence interval
*/

capture program drop auctab
program auctab, rclass
	version 12
	syntax varlist [if] [in] [, text DECimal(integer 2) label(string) local(string) *]
	
	*getting variable label for var
	if trim("`label'")=="" {
		local label: variable label `=word("`varlist'",2)'
		if trim("`label'")=="" local label=trim("`=word("`varlist'",2)'")
	}
	
	*calcluating AUC
	qui roctab `varlist' `if' `in', `options'
	
	*setting format for AUC and CI
	local format="%9.`decimal'f"
	
	*returnning stats adn formatted stats
	foreach stat in area lb ub {
		return local `stat'=`r(`stat')'
		return local `stat'_fmt=trim(string(`r(`stat')',"`format'"))
	}

	*saving out local of formatted AUC	
        if trim("`local'")!="" c_local `local'=trim(string(`r(area)',"`format'"))

	*outputting formatted results
	if trim("`text'")=="" {
		local auc=trim(string(`r(area)',"`format'"))+" ("+trim(string(`r(lb)',"`format'"))+", "+trim(string(`r(ub)',"`format'"))+")"
		return local auc="`auc'"
		disp"`label'&`auc'"
	}
	
	*outputting formatted results to be copied into text
	else {
		local auc="AUC "+trim(string(`r(area)',"`format'"))+" (95% CI "+trim(string(`r(lb)',"`format'"))+", "+trim(string(`r(ub)',"`format'"))+")"
		return local auc="`auc'"
		disp"`label' `auc'"	
	}	
	*returning N and SE
	return local N=`r(N)'
	return local se=`r(se)'




end
