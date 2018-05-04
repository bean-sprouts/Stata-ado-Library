/*
PROGRAM: anovadiff.ado
DESCRIPTION: Outputs difference along with 95% CI and p-value continuous variables
PROGRAMMER: Daniel
DATE: 5/7/2013
*/

capture program drop anovadiff
program anovadiff, rclass
        syntax varlist(min=1 max=1 numeric) [if] [in], by(varname) [covariates(varlist) DECimal(integer 0) SDDECimal(integer 0) DIFFDECimal(string) label(string) header REVerse rtf(string)]
        preserve
        capture keep `if' `in'
        *dropping missing by variables
        qui drop if mi(`by')
		
		*if dffdecimal not specifed, then making it the same as decimal
				if "`diffdecimal'"=="" local diffdecimal=`decimal'

        
        *creating a zero-one version of the by variable
                tempvar bygroup
                qui egen `bygroup'=group(`by')
                qui replace `bygroup'=`bygroup'-1
				
		*1/18/2017: Emily added "reverse" option to flip groups
				if trim("`reverse'")!="" {
					tempvar bygroup2
					qui gen `bygroup2' = 0 if `bygroup'==1
					qui replace `bygroup2' = 1 if `bygroup'==0
					qui replace `bygroup' = `bygroup2'
				}        
        
        *outcome variable must be numeric
                if upper(substr("`:type `by''",1,3)) == "STR" {
                        disp as error "Outcome must be numeric."
                        exit
                }
        
        *if label not specifed making it the variable label or var name
                if "`label'"=="" local label: variable label `varlist'
                if "`label'"=="" local label="`varlist'"
        
        *saving copy of the analysis dataset
                tempfile binarydiff_data
                qui save `binarydiff_data'

        * displaying header row 
                if trim("`header'")!="" {
                        *calculating n in by group
                        qui statsby byN=r(N), by(`bygroup' `by'):  count
                        *if header requested, then making it now.
						*1/18/2017: Emily added code to save out header
						*2/22/2017: Emily added code to use variable and value labels in header, if available.
						*variable label, if available
						local bylabel : var label `by'
						if "`bylabel'"=="" local bylabel = "`by'"
						
						*value label (otherwise just uses value)
						*if upper(substr("`:type `by''",1,3))=="STR" local byval=`by'[`i']
						local trueval1 = `by'[1]
						local trueval2 = `by'[2]
						local byval1 : label (`by') `trueval1'
						local byval2 : label (`by') `trueval2'
						
						local header: disp "&`bylabel'=`byval1' (N=" byN[1] ")&`bylabel'=`byval2' (N=" byN[2] ")&Difference&95% CI&p-value"
						noi disp "`header'"
						return local N0 = byN[1]
						return local N1 = byN[2]
                       * disp "&`by'=" `by'[1] " (N=" byN[1] ")&`by'=" `by'[2] " (N=" byN[2] ")&Difference&95% CI&p-value"
                }
        
        
        **dropping missing variables (as they will not be included in the analysis)
                use `binarydiff_data', clear
                qui drop if mi(`varlist')
                qui save `binarydiff_data', replace

        *ensuring there are still remaining variables
                qui sum `bygroup'
                if `r(max)'>=2 {
                        disp as error "By variable cannot be more than 2 levels"
                        exit
                }
                else if `r(max)'==0 {
                        disp as error "By variable must have at least 2 levels."
                        exit
                }
                
        *calculating difference statistics
                qui regress `varlist' `bygroup' `covariates'
                *saving results as locals
                return local N=e(N)
                return local diff=-_b[`bygroup']
                return local se=_se[`bygroup']
				*1/18/2017: Emily split these into creating locals and returning locals so we can use locals below.                
				local lb=_b[`bygroup']*-1-invttail(e(df_r),.025)*_se[`bygroup']
				return local lb=`lb'
                local ub=_b[`bygroup']*-1+invttail(e(df_r),.025)*_se[`bygroup']
				return local ub=`ub'
                local betafmt=string(real(string(_b[`bygroup']*-1,"%9.`diffdecimal'f")),"%9.`diffdecimal'f")
 				return local betaformat = "`betafmt'"
                local ci=string(real(string(_b[`bygroup']*-1-invttail(e(df_r),.025)*_se[`bygroup'],"%9.`diffdecimal'f")),"%9.`diffdecimal'f")+", "+string(real(string(_b[`bygroup']*-1+invttail(e(df_r),.025)*_se[`bygroup'],"%9.`diffdecimal'f")),"%9.`diffdecimal'f")
 				return local ciformat = "`ci'"
                capture return local pexact=ttail(e(df_r),abs(_b[`bygroup']/_se[`bygroup']))*2
                capture pdisplay ttail(e(df_r),abs(_b[`bygroup']/_se[`bygroup']))*2, local(p)
                return local pformat="`p'"
                *if the variance is not calculable the p-value on CI are NA
                if ttail(e(df_r),abs(_b[`bygroup']/_se[`bygroup']))*2==. {
                        local ci="NA"
                        local p="NA"
                }

        *counting N in each group
                use `binarydiff_data', clear
                qui statsby mean=r(mean) sd=r(sd), by(`bygroup' `by'):  sum `varlist'
                
        *formating stat to display
        qui g stat=string(real(string(mean,"%9.`decimal'f")),"%9.`decimal'f") + " (" + string(sd,"%9.`sddecimal'f") + ")" 
						
		*1/18/2017: Emily added in code to save out information on each group already formatted.
		sort `bygroup'
		return local stat0 = stat[1]
		return local stat1 = stat[2]

        *adding in  difference statistics
        qui g betafmt="`betafmt'"
        qui g ci="`ci'"
        qui g p="`p'"
        qui g label="`label'"

        keep `bygroup' stat betafmt ci p label
        *transposing to one line from one line per by group
        qui reshape wide stat, i(betafmt ci p label) j(`bygroup')
        
        *displaying results
        listtex label stat0 stat1 betafmt ci p, type
		
		*1/18/2017: Emily added code for option to export to RTF
		if trim("`rtf'")!="" {
			rtfrstyle label stat0 stat1 betafmt ci p, cwidths(2500) local(b d e) cdadd("\clbrdrt\brdrs\clbrdrb\brdrs\clbrdrr\brdrs\clbrdrl\brdrs")
			if "`header'"!="" local headerrtf = "`b'\ql " + subinstr("`header'","&","`d'\qc ",.) + "`e'"
			listtab label stat0 stat1 betafmt ci p, handle("`rtf'") begin("`b'\ql ") end("`e'") delim("`d'\qc ") head("`headerrtf'")
		}
		
		
		

end


