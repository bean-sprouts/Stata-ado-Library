/*
Programmer: Jessica Poon
Date: 2/15/2013
*/

// Updated 5/23/2014 by Emily Vertosick: "by" option now lists levels of "by" variable as first column in table,
        // instead of listing separate table for each.
		
// 1/12/2017: Emily fixed issue with percentages not displaying correctly.

*drop any previous programs with the same name
capture program drop widecount
*define new program
program widecount
        *set syntax
        syntax varlist [if] [in] [, table NOSTATistics * by(varlist)]

        preserve
        
        *this allows for the if option
        tempvar touse
        quietly g `touse'=0
        quietly replace `touse'=1 `if' `in' 
        quietly keep if `touse' 
        
        * if there are observations to use, continue.
        quietly count 
        local ntot=`r(N)'
        if `ntot'!=0 {
		
				*if we just want to know whether a variable is available, replace the contents with missing or not
				quietly sort `by' `varlist'
                local tempvarlist ""
                *for each variable in varlist entered
                foreach var in `varlist' {
                        tempvar temp`var'
                        *create a binary for missing or not
                        qui gen `temp`var''= "." if mi(`var')
                        qui replace `temp`var'' = "available" if !mi(`var')
                        local tempvarlist = "`tempvarlist' `temp`var''"
                        label var `temp`var'' "`var'"
                        drop `var'
                        rename `temp`var'' `var'
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
                        qui levelsof `by', local(levels) missing
                        bysort `by': g _Total = _N
                        local sepopt "sepby(`by')"
                        local bytable "`by'&"
                }
                
                *sort and capture frequency of each varlist combination
                quietly sort `by' `varlist'
                quietly by `by' `varlist': g _Freq=_N 
                quietly by `by' `varlist': keep if _n==_N
                                                
                qui g _Percentn=_Freq/_Total*100
                format _Percentn %9.1f
			
                *trimming the percentage
                qui g _Percent=trim(string(_Percentn,"%9.0f"))+"%" if round(_Percentn,1)
                qui replace _Percent=trim(string(_Percentn,"%9.1f"))+"%" if round(_Percentn,0.1)<10
                qui replace _Percent="<0.1%" if round(_Percentn,0.1)<0.1
                
                *if user wants frequencies reported naming the variables here
                if trim("`nostatistics'")=="" local stats=" _Freq _Percent"
                
                *if table not requested, print as normal Stata window table
                if trim("`table'")=="" { 
                        list `by' `varlist' `stats', noobs `options' labvar(`label') `sepopt'
                }
                
                *otherwise, print with ampersands for table request
                else {
                        local header=trim("`=subinstr("`bytable'`varlist'`stats'"," ","&",.)'")
                        listtex `by' `varlist' `stats', type headlines(`header')
                }
        }
                restore

end
exit
