/*
PROGRAM: pdisplay.ado
PROGRAMMER: Daniel
DATE: 2/8/2012
DESCRIPTION: This is a simple program that displays a formatted p-value and 
				returns it as a local variable.  
				The Bonferroni correction can also be applied.
*/

capture program drop pdisplay
program pdisplay, rclass
	version 11
	syntax anything [, local(string) BONferroni(integer 1) MINimum(real 0.0001) p]
	
	*checking to ensure p-value is between 0 and 1.
	*using the rounded version because Stata can sometimes print out p-values that are slightly larger than one for exact test (eg 1.00000000000000006)
	if !inrange(round(`anything',0.000001),0,1) {
		disp as err "p-value not between 0 and 1"
		exit
	}
	
	*Minumum p-value cannot be less than 0.000001
	if !inrange(`minimum',0.000001,1) {
		disp as err "Minimum p-value not be less than 0.000001 or greater than 1."
		exit
	}
	
	/*performing bonferroni correction to pvalue*/
	local pexact=`anything'*`bonferroni'

	*now format p depending on how large or small it is
	if `pexact'>=.95 local pformat="1"
	else if round(`pexact',0.1)>=0.2 local pformat: disp %9.1f `pexact'
	else if round(`pexact',0.01)>=0.1 local pformat: disp %9.2f `pexact'
	else if round(`pexact',0.001)>=0.001 local pformat: disp %9.3f `pexact'
	else if round(`pexact',0.0001)>=0.0001 local pformat: disp %9.4f `pexact'
	else if round(`pexact',0.00001)>=0.00001 local pformat: disp %9.5f `pexact'
	else if round(`pexact',0.000001)>=0.000001 local pformat: disp %9.6f `pexact'

	
	if `pexact'<`minimum' {
		local pformat="<0`minimum'"
		local fullpformat="p"+trim("`pformat'")
	}
	else local fullpformat="p="+trim("`pformat'")

	if trim("`p'")=="" disp trim("`pformat'")
	else disp trim("`fullpformat'")
	
	return local pformat=trim("`pformat'")
	return local fullpformat=trim("`fullpformat'")

	if trim("`local'")!="" {
		if trim("`p'")=="" c_local `local'=trim("`pformat'")
		else  c_local `local'=trim("`fullpformat'")
	}

end

