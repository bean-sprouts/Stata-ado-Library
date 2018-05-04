
// Updated 6/29/2015 by Jessica Poon Emily Vertosick: "by" option now lists levels of "by" variables as first column in table,
// instead of listing separate table for each. Also, "available" option now allows us to see pattern of missing values.
// 3/13/2017: Updated by Emily to allow for tables from widetab to be exported directly to RTF.

capture program drop widetab
program widetab
        qui {
        preserve
        syntax varlist [if] [in] [, table NOSTATistics group * by(varlist) AVAIlable rtf(string)] 
        
        /*Not sure why the marksample command is not working*/
        /*marksample touse*/
        tempvar touse
        g `touse'=0
        replace `touse'=1 `if' `in'
        keep if `touse'

        quietly count 
        local ntot=`r(N)'
        if `ntot'!=0 {
        
                *if we just want to know whether a variable is available, replace the contents with missing or not
                if trim("`available'")!="" {
                        sort `by' `varlist'
                        local tempvarlist ""
                        *for each variable in varlist entered
                        foreach var in `varlist' {
                                tempvar temp`var'
                                *create a binary for missing or not
                                gen `temp`var''= "." if mi(`var')
                                replace `temp`var'' = "available" if !mi(`var')
                                local tempvarlist = "`tempvarlist' `temp`var''"
                                label var `temp`var'' "`var'"
                                drop `var'
                                rename `temp`var'' `var'
                        }
                }
                
                *if by variable isn't specified, report for all
                if trim("`by'")=="" {
                        g _Total = _N
                        local levels = 1
                        local sepopt ""
                        local bytable ""
                }
                
                *if by variable is specified, report missing/non-missing by each level of by
                *local sepopt sets option to list with separators between levels of `by' (instead of every 5 obs)
                *local bytable is set so that if printing using "table" option, the name of the "by" variable will print correctly in the header.
                if trim("`by'")!="" {
                        levelsof `by', local(levels) missing
                        bysort `by': g _Total = _N
                        local sepopt "sepby(`by')"
                        local bytable "`by'&"
                }

                *calculating the stats for each by level
                sort `by' `varlist'
                by `by' `varlist': g _Freq=_N 
                by `by' `varlist': keep if _n==_N
                g _Percentn=_Freq/_Total*100
                format _Percentn %9.1f
                
                *trimming the percentage
                g _Percent=trim(string(_Percentn,"%9.0f"))+"%" if round(_Percentn,1)
                replace _Percent=trim(string(_Percentn,"%9.1f"))+"%" if round(_Percentn,0.1)<10
                replace _Percent="<0.1%" if round(_Percentn,0.1)<0.1
                
                *simplifying values if group specified
                if trim("`group'")!="" {
                        foreach v of varlist `varlist' {
                                local cumvarlist="`cumvarlist' `v'"
                                capture by `by' `cumvarlist': replace `v'="" if _n!=1
                                capture by `by' `cumvarlist': replace `v'=.a if _n!=1
                        }
                }
                
                *if user wants frequencies reported naming the variables here
                if trim("`nostatistics'")=="" local stats=" _Freq _Percent"
                
                *if table not requested, print as normal Stata window table
                noi if trim("`table'")=="" list `by' `varlist' `stats', noobs `options' `sepopt' 
                
                *otherwise, print the table with ampersands
                else {
                        local header=trim("`=subinstr("`bytable'`varlist'`stats'"," ","&",.)'")
                        noi listtex `by' `varlist' `stats', type headlines(`header')
                }
				
				*RTF export: if RTF option is specified, export to RTF
				if trim("`rtf'")!="" {
				
					// First, save out variable labels for header
					local headerlab = ""
					local count = 0
					foreach v in `varlist' {
						
						// Extract variable label, if it exists
						local varlab : var label `v'
						
						// If no variable label, use variable name
						if trim("`varlab'")=="" local varlab = "`v'"
						
						if `count'==0 local varlab = "`varlab'"
							// For the first column, don't put "&" before variable name because this messes up format of table in RTF
						else if `count' > 0 local varlab = "&`varlab'"

						// Add to list of variable labels for header
						local headerlab = "`headerlab'`varlab'"
						local ++count
					}
					
					// If there is a by variable, save out the variable label for the by variable too.
					if trim("`by'")!="" {
						
						// Extract variable label, if it exists
						local bylab : var label `by'
						local bylab = "`bylab'&"
						
						// If no variable label, use variable name
						if trim("`bylab'")=="" local bylab = "`by'"
						
					}
					
					else if trim("`by'")=="" local bylab = ""
					file write `rtf' "{\pard\par}"
					rtfrstyle `by' `varlist' `stats', cwidths(2500) local(b d e) cdadd("\clbrdrt\brdrs\clbrdrb\brdrs\clbrdrr\brdrs\clbrdrl\brdrs")
					local statsheader=trim("`=subinstr("`stats'"," ","&",.)'")
					local headerrtf="`bylab'`headerlab'`statsheader'"
					local headerrtf = subinstr("`headerrtf'","&","`d'\qc ",.)
					local headerrtf = "`b'\ql " + "`headerrtf'" + "`e'"
					listtab `by' `varlist' `stats', handle("`rtf'") begin("`b'\ql ") end("`e'") delim("`d'\qc ") head("`headerrtf'") missnum("Missing")
				}
				
                }
        restore
        }
end
exit
