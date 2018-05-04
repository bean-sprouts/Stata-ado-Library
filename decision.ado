/*
PROGRAM: decision.ado
DESCRIPTION: this program inputs a decicsion and tabulates the results based on that decision
*/

capture program drop decision
program decision, rclass
	syntax anything [if] [in] [pweight], outcome(varname)  [header order(string) label(string) nper(numlist) format(string)]
	preserve
	*subsetting data
	capture qui keep `if' `in'
	
	**  setting default format
	if "`format'"=="" local format="%9.2f"
	
	**  asserting that each variable listed in order option is a acceptable variable
	if length("`order'")>0 {
		foreach v in `order' {
			capture assert inlist("`v'","testpos","testneg","pos","neg","tp","fn","tn","fp","N") | inlist("`v'","sens","spec","ppv","npv")
			if _rc>0 {
				disp as error "`v' not allowed.  Must be N, Nper, testpos, testneg, pos, neg, tp, fn, tn, fp, sens, spec, ppv, and/or npv."
				exit
			}
		}
	}
	*assigning order vars appear in if none specified
	if "`order'"=="" local order="Nper testpos testneg pos neg tp fn tn fp sens spec ppv npv"
	
	tempvar decision 
	qui g `decision'=`anything'
	
	*the decision variable must be binary and not missing
	capture assert inlist(`decision',0,1) & inlist(`outcome',0,1)
	if _rc>0 {
		disp as error "There can be no missing values for the outcome or the decision variable."
		exit
	}
	
	tempvar testpos testneg pos neg tp fn tn fp 
	qui g `testpos'=`decision'==1
	qui g `testneg'=`decision'==0
	qui g `pos'=`outcome'==1
	qui g `neg'=`outcome'==0
	qui g `tp'=`outcome'==1 & `decision'==1
	qui g `tn'=`outcome'==0 & `decision'==0
	qui g `fn'=`outcome'==1 & `decision'==0
	qui g `fp'=`outcome'==0 & `decision'==1

	
	tempname memhold
	tempfile results
	postfile `memhold' str10(statname) stat using `results'
	post `memhold' ("N") (`=_N')
	if "`nper'"=="" post `memhold' ("Nper") (`=_N')
	else post `memhold' ("Nper") (`nper') 
	
	qui foreach stat in testpos testneg pos neg tp fn tn fp {
		mean ``stat'' [`weight' `exp']
			tempname mean
			matrix `mean'=r(table)
			post `memhold' ("`stat'") (`mean'[1,1])
	}
	postclose `memhold'
	use `results', clear
	g id=1
	
	qui reshape wide stat, i(id) j(statname) string
	
	rename (stat*) (*)
	drop id
	qui foreach stat in testpos testneg pos neg tp fn tn fp {
		replace `stat'=`stat'*Nper
	}
	


	qui g sens=tp/pos
	qui g spec=tn/neg
	qui g ppv=tp/testpos
	qui g npv=tn/testneg
	
	*label and formatting
		label var pos "Positive"
		label var neg "Negative"
		label var sens "Sensitivity"
		label var spec "Specificity"
		label var ppv "PPV"
		label var npv "NPV"
		label var testpos "Test Positive"
		label var testneg "Test Negative"
		label var tp "True Positive"
		label var tn "True Negative"
		label var fp "False Positive"
		label var fn "False Negative"
		label var N "N"
		label var Nper "Standardized N per"

		format `format' sens spec ppv npv 
		format %9.0f testpos testneg pos neg tp fn tn fp
	
	*returning results as locals
	foreach v of varlist * {
		return local `v'=`v'
	}

	*assigning label if none indicated
	if "`label'"=="" local label="`anything'"
	qui g label=`"`label'"'
	label var label "Label"
	order label
	
	*creating header if requested
	foreach v of varlist label `order' {
		local headerdisp="`headerdisp'`:var label `v''&"
	}
	local headerdisp=substr("`headerdisp'",1,length("`headerdisp'")-1)

	*displaying results
	if "`header'"=="" listtex label `order', type
	else listtex label `order', type headlines("`headerdisp'")

	restore 
end

