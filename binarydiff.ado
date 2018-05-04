/*
PROGRAM: binarydiff.ado
DESCRIPTION: Outputs difference along with 95% CI and p-value binary variables
PROGRAMMER: Daniel
DATE: 5/7/2013
*/

capture program drop binarydiff
program binarydiff, rclass
        syntax varlist(min=1 max=1 numeric) [if] [in], by(varname) [label(string) header exact REVerse rtf(string)]
        preserve
        capture keep `if' `in'
        *dropping missing by variables
        qui drop if mi(`by')
        
        *creating a zero-one version of the by variable
                tempvar bygroup
                qui egen `bygroup'=group(`by')
                qui replace `bygroup'=`bygroup'-1
                qui sum `bygroup'
				
		*1/18/2017: EV added option REVerse which allows you to flip by groups (rather than having to flip the signs afterward)
				if trim("`reverse'")!="" {
					tempvar bygroup2
					qui g `bygroup2' = 0 if `bygroup'==1
					qui replace `bygroup2' = 1 if `bygroup'==0
					qui replace `bygroup' = `bygroup2'
				}
        
        *outcome variable must be coded as 0-1
                qui levelsof `varlist'
                local outcomelevels "`r(levels)'"
                capture assert inlist(trim("`outcomelevels'"),"0 1","1","0")
                if _rc>0 {
                        disp as error "Outcome must be coded as 0 or 1."
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
                        disp "`header'"
						return local N0 = byN[1]
						return local N1 = byN[2]
						*disp "&`by'=" `by'[1] " (N=" byN[1] ")&`by'=" `by'[2] " (N=" byN[2] ")&Difference&95% CI&p-value"
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
                *only performing calculation if there is varaiation in the outcome.
                if inlist("`outcomelevels'","0","1") local rd="0.0%"
                qui else {
                        cs `varlist' `bygroup', `exact'
                        *saving results as locals
                        foreach stat in lb_rd ub_rd rd {
                                if round(abs(`r(`stat')'),0.1)>=0.1 local `stat'=string(`r(`stat')'*-100,"%9.0f")+"%"
                                else local `stat'=string(`r(`stat')'*-100,"%9.1f")+"%"
                        }
                        *return local pexact=`r(p)'
                        if trim("`exact'")=="" qui pdisplay r(p), local(p)
                        else qui pdisplay r(p_exact), local(p)
                        local ci="`ub_rd', `lb_rd'"
                }
        
        *calculating n in by group
                qui statsby byN=r(N), by(`bygroup' `by'):  count
                tempfile byfreqs
                qui save `byfreqs'
                
                *adding a line for varlist=1 to ensure there will be a line in the final dataset
                g `varlist'=1
                tempfile byfreqsvar
                qui save `byfreqsvar'
                
        *counting N in each group
        use `binarydiff_data', clear
        qui statsby N=r(N), by(`bygroup' `by' `varlist'):  count
        
        * merging byvariable counts to within varlist counts to calculate frequencies
        qui merge m:1 `bygroup' `by' using `byfreqs', nogen
        qui merge m:1 `bygroup' `by' `varlist' using `byfreqsvar', 
        qui replace N=0 if mi(N)
        
        *keeping one line per by group to report (where varlist==1)
        qui keep if `varlist'==1
        qui g stat=string(N)+" ("+string(N/byN*100,"%9.0f")+"%)" if round(N/byN*100,1)>=10
        qui replace stat=string(N)+" ("+string(N/byN*100,"%9.1f")+"%)" if !(round(N/byN*100,1)>=10)
        
        keep `bygroup' `varlist' stat
        *transposing to one line from one line per by group
        qui reshape wide stat, i(`varlist') j(`bygroup')
        *adding in risk difference statistics
        qui g rd="`rd'"
        qui g ci="`ci'"
        qui g p="`p'"
        qui g label="`label'"
        
        *if everyone in the cohort, no one in the cohort, had the outcome, the p-value and CI are not calculable.
        qui if inlist("`outcomelevels'","0","1") {
                replace p="NA"
                replace ci="NA"
        }
        
        *displaying results
        listtex label stat0 stat1 rd ci p, type
		
		*1/18/2017: Emily added code for option to export to RTF
		if trim("`rtf'")!="" {
			rtfrstyle label stat0 stat1 rd ci p, cwidths(2500) local(b d e) cdadd("\clbrdrt\brdrs\clbrdrb\brdrs\clbrdrr\brdrs\clbrdrl\brdrs")
			if "`header'"!="" local headerrtf = "`b'\ql " + subinstr("`header'","&","`d'\qc ",.) + "`e'"
			listtab label stat0 stat1 rd ci p, handle("`rtf'") begin("`b'\ql ") end("`e'") delim("`d'\qc ") head("`headerrtf'")
		}
        
        return local label=label
        return local stat0=stat0
        return local stat1=stat1
        return local rd=rd
        return local ci=ci
        return local p=p

end

