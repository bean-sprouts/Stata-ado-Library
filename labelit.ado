/*
THE FORMATS PROGRAM CREATES A FILE THAT CAN BE RUN IN 
STATA THAT WILL CRAETE THE SPECIFIED FORMATS

INPUT 1: Source file with labels and formats defined
*/

capture program drop labelit
program labelit, rclass
	version 12
	syntax using [, sheet(string)]
	preserve

	if trim(`"`sheet'"')=="" local sheet="Derived Variables"
	import excel `using', firstrow allstring clear sheet(`"`sheet'"')
	quietly keep if trim(varname)!= "" & upper(codeadded)=="X"
		
	*CREATING A RETURN LIST VARIABLE WITH ALL DERIVED VARIABLE NAMES
	quietly vallist varname
	return local drvvarlist `r(list)'

	*CREATING VAR LBL CODE TO APPLY VARIABLE LABELS
	qui g labvarcmd=`"label variable "'+varname+`" ""'+label+`"""'

	*CREATING VARIABLE FOR APPLYING VALUE LABELS
	qui g labvalcmd=`"label value "'+varname+`" "'+formatname if strpos(formatname,"%")==0 & trim(formatname)!=""
	qui replace labvalcmd=`"format "'+varname+`" "'+formatname if strpos(formatname,"%")>0 & trim(formatname)!=""
		
	*number of variables
	local varn=_N
	
	foreach v of numlist 1/`varn' {
		local varlbl`v'=labvarcmd[`v']
		local varval`v'=labvalcmd[`v']
	}
	
	*importing formats
	import excel `using', firstrow allstring clear sheet("Formats")
	qui drop if missing(formatname)
	qui bysort formatname: keep if _n==_N
	qui g labdefcmd=`"label define "'+formatname+`" "'+format+`" , replace"'

	*number of formats
	local fmtn=_N
	
	*cycling through each value label and applying code
	*set trace on
	foreach v of numlist 1/`fmtn' {
		local fmt`v'=labdefcmd[`v']
	}	
	set trace off
	
	*restoring datatset
	restore
	
	foreach v of numlist 1/`fmtn' {	
		disp `"`fmt`v''"'
		*Applying the value label code!
		`fmt`v''
	}	

	*cycling through each variable and applying label and value label	
	foreach v of numlist 1/`varn' {
		disp `"`varlbl`v''"'
		disp `"`varval`v''"'
		`varlbl`v''
		`varval`v''
	}
end	

