
capture program drop tableone
program tableone, rclass
		version 12
		syntax varlist [if] [in] ,  [type(string) by(varname) label(string) header Missing ///
						stattype(integer 1) DECimal(integer 1) test(string) DELIMiter(string) ///
						VARVALue(string asis) statlabel FREQuency rtf(string)]		
		preserve
		
		* Default delimiter &;
		if `"`delimiter'"'=="" local delimiter "&"

		/*Subsetting dataset on selected observations*/
		capture keep `if' `in'

		/*Identifying type of statistics to calculate*/
		local type=trim(upper(substr("`type'",1,3)))
		if !inlist(trim("`type'"),"DIC","CON","CAT","") {
				disp as err "`v': Type must be Dichotomous, Categorical, or Continuous."
				exit
		}

		/*creating a temporary "by" variable if none is indicated*/             
		tempvar by_ord
		if trim("`by'")=="" g `by_ord'=1
		else egen `by_ord'=group(`by')

		
		/*Stats are only calculated for non-missing by variables*/
		qui keep if !mi(`by_ord')
		quietly count
				local N=r(N)    
				return local N=`N'
				*saving a copy of the analysis dataset
				tempfile analysisdata
				qui save `analysisdata'

		
		/*creating dataset that is one line per byvariable level*/
		qui {
				tempfile byNds
				tempvar byN idone
				g `idone'=1
				statsby `byN'=r(N), by(`idone' `by_ord' `by') saving(`byNds'): count
		}
		
		*displaying header row
		if trim("`header'")!="" {
				if "`statlabel'"!="" local byline="&Statistic"
				qui use `byNds', clear
				tempvar totN
				qui egen `totN'=sum(`byN')
				if trim("`by'")=="" local byline="`byline'&N="+string(`totN')
						else {
						foreach i of numlist 1/`=_N' {
								// 2/22/2017: Emily added code to use variable and value labels in header, if available.
								// Variable label, if available
								local bylabel : var label `by'
								if "`bylabel'"=="" local bylabel = "`by'"
								
								// Value label (otherwise just uses value)
								if upper(substr("`:type `by''",1,3))=="STR" local byval=`by'[`i']
								local trueval = `by'[`i']
								else local byval : label (`by') `trueval'
								*else local byval=string(`by'[`i'])
								
								if round(`byN'[`i']/`totN'*100,1)>=10 local byline="`byline'&`bylabel'=`byval' (N="+string(`byN'[`i'])+"; "+string(`byN'[`i']/`totN'*100,"%9.0f")+"%)"
								else if round(`byN'[`i']/`totN'*100,1)<10 local byline="`byline'&`bylabel'=`byval' (N="+string(`byN'[`i'])+"; "+string(`byN'[`i']/`totN'*100,"%9.1f")+"%)"
								return local N`i' = `byN'[`i']
						}
				}
				disp "`byline'"
				local headerrtf = "`byline'"
		}
		if trim("`header'")=="" local headerrtf = ""		

		*looping over every variable in list
		foreach v in `varlist' {
		
				*restoring analytic data
				use `analysisdata', clear

				/*assigning the variable label if not specified*/
				local vlabel=`"`label'"'
				if trim("`vlabel'")=="" local vlabel: variable label `v'
				if trim("`vlabel'")=="" local vlabel="`v'"
				
				/*if type is not specified, then guessing type of variable*/
				if trim("`type'")=="" {
						qui levelsof `v'
						local nlevels: word count `r(levels)'
						if upper(substr("`:type `v''",1,3))=="STR" local vtype="CAT"
						else if inlist(trim("`r(levels)'"),"0 1","1","0") local vtype="DIC"
						else if `nlevels'<10 local vtype="CAT"
						else local vtype="CON"                  
				}
				else local vtype="`type'"
				
		
				****************************************************
				**************  CONTINUOUS VARIABLES  **************
				****************************************************
				if "`vtype'"=="CON" {
								/*calculating p-values if requested*/
										if trim("`test'")=="ranksum" {
												qui ranksum `v', by(`by_ord')
														qui pdisplay 2*(1-normal(abs(r(z)))), local(temp_p)
										}
										else if trim("`test'")=="ttest" {
												qui ttest `v', by(`by_ord')
														qui pdisplay r(p), local(temp_p)
										}
										else if trim("`test'")=="kwallis" {
												qui kwallis `v', by(`by_ord')
														qui pdisplay 1-chi2(`r(df)',`r(chi2_adj)'), local(temp_p)
										}
										else if trim("`test'")=="regress" {
												qui regress `v' i.`by_ord'
														qui pdisplay 1-F(`e(df_m)',`e(df_r)',`e(F)'), local(temp_p)
										}
										else if trim("`test'")!="" {
												disp as err `"Test must be "ranksum", "ttest", "kwallis", or "regress" for continuous variables."'
										}       
								
								
								/*calculating stats by group*/
								qui {
										tempfile bystats
										tempvar N mean var sd min max p25 p50 p75 byNtot Ntot
										statsby `N'=r(N) `mean'=r(mean) `var'=r(Var) `sd'=r(sd) `min'=r(min) `max'=r(max) `p25'=r(p25) `p50'=r(p50)`p75'=r(p75), by(`by_ord' `by') saving(`bystats'): sum `v', d
										use `bystats'
										merge 1:1 `by_ord' `by' using `byNds', nogen
										
										egen `byNtot'=sum(`byN')
										egen `Ntot'=sum(`N')
								}
								
								
								/*formatting stats*/
								tempvar stat vstatlabel
								local format="%9.`decimal'f"
								if `stattype'==1 {
										qui g `vstatlabel'="Median (IQR)"
										qui g `stat'=string(`p50',"`format'")+" ("+ ///
													 string(`p25',"`format'")+", "+ ///
													 string(`p75',"`format'")+")"
										qui replace `stat'="NA" if inlist(`N',0,.)
								}
								else if `stattype'==2 {
										qui g `vstatlabel'="Mean (SD)"
										qui g `stat'=string(`mean',"`format'")+" ("+string(`sd',"`format'")+")"
										qui replace `stat'=string(`mean',"`format'")+" (NA)" if `N'==1
										qui replace `stat'="NA" if inlist(`N',0,.)
								}
								else if `stattype'==3 {
										qui g `vstatlabel'="Median (IQR) Min:Max"
										qui g `stat'=string(`p50',"`format'")+" ("+ ///
													 string(`p25',"`format'")+", "+ ///
													 string(`p75',"`format'")+") "+ ///
													 string(`min',"`format'")+":"+ ///
													 string(`max',"`format'")											
										qui replace `stat'="NA" if inlist(`N',0,.)
								}
								else if `stattype'==4 {
										qui g `vstatlabel'="Coefficient of Variation"
										qui g `stat'=string(`sd'/`mean'*100,"`format'")
										qui replace `stat'="NA" if inlist(`N',0,.) | `N'==1
								}
								else {
										disp as error "stattype not found"
										exit
								}
								
								/*assigning variable label*/
								tempvar labelvar
								qui g `labelvar'="`vlabel'"
								qui replace `labelvar'=`labelvar'+" (N="+string(`Ntot')+")" if `byNtot'!=`Ntot'
								
								*reshaping to one line to print
								keep `labelvar' `vstatlabel' `stat' `by_ord'
								qui reshape wide `stat', j(`by_ord') i(`labelvar' `vstatlabel') 
								
								/*displaying results*/								
								if trim("`temp_p'")!="" {
										tempvar p
										qui g `p'="`temp_p'" in 1
										if "`statlabel'"=="" listtex `labelvar' `stat'* `p', type  delimiter("`delimiter'")
										else listtex `labelvar' `vstatlabel' `stat'* `p', type  delimiter("`delimiter'")

								}
								else if "`statlabel'"==""  listtex `labelvar' `stat'*, type  delimiter("`delimiter'")
								else listtex `labelvar' `vstatlabel' `stat'*, type  delimiter("`delimiter'")
								
								// If exporting to RTF, do that here.
								if trim("`rtf'")!="" {
									if trim("`temp_p'")!="" {
										tempvar p
										qui g `p'="`temp_p'" in 1
										if "`statlabel'"=="" {
											rtfrstyle `labelvar' `stat'* `p', cwidths(2500) local(b d e) cdadd("\clbrdrt\brdrs\clbrdrb\brdrs\clbrdrr\brdrs\clbrdrl\brdrs")
											if "`header'"!="" local headerrtf = "`b'\ql " + subinstr("`byline'","&","`d'\qc ",.) + "`d'\qc " + "p value" + "`e'"
											listtab `labelvar' `stat'* `p', handle("`rtf'") begin("`b'\ql ") end("`e'") delim("`d'\qc ") head("`headerrtf'")
										}
										else {
											rtfrstyle `labelvar' `vstatlabel' `stat'* `p', cwidths(2500) local(b d e) cdadd("\clbrdrt\brdrs\clbrdrb\brdrs\clbrdrr\brdrs\clbrdrl\brdrs")
											if "`header'"!="" local headerrtf = "`b'\ql " + subinstr("`byline'","&","`d'\qc ",.) + "`d'\qc " + "p value" + "`e'"
											listtab `labelvar' `vstatlabel' `stat'* `p', handle("`rtf'") begin("`b'\ql ") end("`e'") delim("`d'\qc ") head("`headerrtf'")
										}
									}
									else if "`statlabel'"=="" {
										rtfrstyle `labelvar' `stat'*, cwidths(2500) local(b d e) cdadd("\clbrdrt\brdrs\clbrdrb\brdrs\clbrdrr\brdrs\clbrdrl\brdrs")
										if "`header'"!="" local headerrtf = "`b'\ql " + subinstr("`byline'","&","`d'\qc ",.) + "`e'"
										listtab `labelvar' `stat'*, handle("`rtf'") begin("`b'\ql ") end("`e'") delim("`d'\qc ") head("`headerrtf'")
									}
									else {
										rtfrstyle `labelvar' `vstatlabel' `stat'*, cwidths(2500) local(b d e) cdadd("\clbrdrt\brdrs\clbrdrb\brdrs\clbrdrr\brdrs\clbrdrl\brdrs")
										if "`header'"!="" local headerrtf = "`b'\ql " + subinstr("`byline'","&","`d'\qc ",.) + "`e'"
										listtab `labelvar' `vstatlabel' `stat'*, handle("`rtf'") begin("`b'\ql ") end("`e'") delim("`d'\qc ") head("`headerrtf'")
									}									
								}

								*returning results as locals
								foreach rown of numlist 1/`=_N' {
										local coln=0
										foreach returnv of varlist `stat'* {
												local ++coln
												return local `v'_`rown'_`coln'=`returnv'[`rown']
										}
								}
								if trim("`temp_p'")!="" return local p="`temp_p'"                                       
								
				} /* end of continuous variable*/
				
				*****************************************************
				**************  CATEGORICAL VARIABLES  **************
				*****************************************************
				if inlist("`vtype'","CAT","DIC") {
								/*calculating p-value if requested*/
								if trim("`test'")=="chi2" {
										qui tab `v' `by_ord', chi2
												qui pdisplay r(p), local(temp_p)
								}
								else if trim("`test'")=="exact" {
										qui tab `v' `by_ord', exact
												qui pdisplay r(p_exact), local(temp_p)
								}
								else if trim("`test'")!="" {
										disp as err `"Test must be "chi2" or "exact" for categorical variables."'
								}
				
								/*creating label that is displayed for levels of variable*/
								/*getting the value label*/
								local labval: value label `v'
								tempvar vdisp
								if upper(substr("`:type `v''",1,3)) == "STR" {
										qui g `vdisp'="    "+`v'
								}
								else {
										if trim("`labval'")=="" {
												qui g `vdisp'="    "+string(`v')
										}
										else {
												qui g `vdisp'=""
												qui levelsof `v'
												foreach i of numlist `r(levels)' {
														local labval`i': label (`v') `i'
														qui replace `vdisp'="    "+trim("`labval`i''") if `v'==`i'                      
												}
										}
								}
								if trim("`missing'")!="" qui replace `vdisp'="    Unknown" if missing(`v')
								label var `vdisp' "vdisp"

								/*calculating stats by group*/
								qui {
										tempfile bystats
										tempvar N byNtot byNtot0 Ntot v_ord vmissing vtotal 
										g `idone'=1
										*dropping missing values of v if not requested in tabulation
										if "`missing'"=="" drop if mi(`v')
										
										g `vmissing'=missing(`v')
										egen `v_ord'=group(`vmissing' `v'), `missing'
										bysort `v': g `vtotal'=_N
										
										if "`missing'"!="" local missingwcomma=", missing"
										statsby `N'=r(N), by(`idone' `by_ord' `by' `vmissing' `vtotal' `v' `v_ord' `vdisp' `missingwcomma') saving(`bystats'): count
										*saving dataset that is one line per variable per by variable
										use `bystats', clear
										keep `idone' `vmissing' `vtotal' `v' `v_ord' `vdisp'
										duplicates drop
										
										joinby `idone' using `byNds'
										tempfile allcombinations
										save `allcombinations'
										
										
										*merging counts with every possible combination of v and byvar
										use `bystats', clear
										merge m:1 `idone' `by_ord' `by' `vmissing' `vtotal' `v' `v_ord' `vdisp' using `allcombinations', nogen
										replace `N'=0 if mi(`N')
										bysort `by_ord' `by': egen `byNtot'=sum(`N')

										egen `Ntot'=sum(`N')
										label var `byN' "(byN) Number of Observations within by group (including any missing values in `v')"
										label var `N' "(N) Number of Observations within by group and `v' (from statsby)"
										label var `byNtot' "(byNtot) Number of observed values with `v' groups"
										label var `Ntot' "(Ntot) Number of observed values"

								}
								
										
								**  if missing values, adding n to variable label
								capture assert `byNtot'==`byN'
								if _rc>0 local vlabel=`"`vlabel' (N="'+string(`Ntot')+")"								
								
								/*formatting stats*/
								tempvar stat vstatlabel
								if `stattype'==1 {
										qui g `vstatlabel'="n (%)"
										qui g `stat'=string(`N')+" (<0.1%)"
										qui replace `stat'=string(`N')+" ("+string(`N'/`byNtot'*100,"%9.1f")+"%)" if `N'/`byNtot'*100>=0.1
										qui replace `stat'=string(`N')+" ("+string(`N'/`byNtot'*100,"%9.0f")+"%)" if round(`N'/`byNtot'*100,1)>=10
										qui replace `stat'=string(`N')+" (0%)" if `N'==0
								}
								else if `stattype'==2 {
										qui g `vstatlabel'="N"
										qui g `stat'=string(`N')
								}
								else if `stattype'==3 {
										qui g `vstatlabel'="%"
										qui g `stat'="<0.1%"
										qui replace `stat'=string(`N'/`byNtot'*100,"%9.1f")+"%" if `N'/`byNtot'*100>=0.1
										qui replace `stat'=string(`N'/`byNtot'*100,"%9.0f")+"%" if round(`N'/`byNtot'*100,1)>=10
										qui replace `stat'=string(`N')+"0%" if `N'==0
								}
								else if `stattype'==4 {
										qui g `vstatlabel'="n %"
										qui g `stat'=string(`N')+" <0.1%"
										qui replace `stat'=string(`N')+" "+string(`N'/`byNtot'*100,"%9.1f")+"%" if `N'/`byNtot'*100>=0.1
										qui replace `stat'=string(`N')+" "+string(`N'/`byNtot'*100,"%9.0f")+"%" if round(`N'/`byNtot'*100,1)>=10
										qui replace `stat'=string(`N')+" 0%" if `N'==0
								}
								else if `stattype'==5 {
										qui g `vstatlabel'="n/N (%)"
										qui g `stat'=string(`N')+"/"+string(`byNtot')+" (<0.1%)"
										qui replace `stat'=string(`N')+"/"+string(`byNtot')+" ("+string(`N'/`byNtot'*100,"%9.1f")+"%)" if `N'/`byNtot'*100>=0.1
										qui replace `stat'=string(`N')+"/"+string(`byNtot')+" ("+string(`N'/`byNtot'*100,"%9.0f")+"%)" if round(`N'/`byNtot'*100,1)>=10
										qui replace `stat'=string(`N')+"/"+string(`byNtot')+" (0%)" if `N'==0
								}
								else {
										disp as error "stattype not found"
										exit
								}
								

								*changing to a wide format.
								keep `by_ord' `vmissing' `vtotal' `v' `v_ord' `vdisp' `vstatlabel' `stat'
								
								qui reshape wide `stat', j(`by_ord') i(`vstatlabel' `vmissing' `vtotal' `v_ord' `v' `vdisp') 
								if "`frequency'"=="" sort `v_ord'
								else if "`frequency'"!="" gsort `vmissing' -`vtotal' `v_ord'

								
								***  for categorical variables displaying row for label
								if "`vtype'"=="CAT" & "`varvalue'"=="" disp `"`vlabel'"'
								qui else {
										*if dichotomous type, then keeping line where variable is 1
										if "`vtype'"=="DIC" & "`varvalue'"=="" local varvalue=1
										if upper(substr("`:type `v''",1,3))=="STR" {
												keep if inlist(`v',"`varvalue'") | mi(`v')
												replace `vdisp'=`"`vlabel'"' if inlist(`v',"`varvalue'")
										}
										else {
												keep if inlist(`v',`varvalue') | mi(`v')
												replace `vdisp'=`"`vlabel'"' if inlist(`v',`varvalue')
										}
								}

								/*displaying results*/
								if trim("`temp_p'")!="" {
										tempvar p
										qui g `p'="`temp_p'" in 1
										if "`statlabel'"=="" listtex `vdisp' `stat'* `p', type  delimiter("`delimiter'")
										else listtex `vdisp' `vstatlabel' `stat'* `p', type  delimiter("`delimiter'")

								}
								else if "`statlabel'"=="" listtex `vdisp' `stat'*, type  delimiter("`delimiter'")
								else listtex `vdisp' `vstatlabel' `stat'*, type  delimiter("`delimiter'")
								
								// If exporting as RTF
								if trim("`rtf'")!="" {

									// When exporting to RTF, we need to save out the variable name/label in vdisp or it won't display otherwise.
									if trim("`rtf'")!="" & "`vtype'"=="CAT" {									
										local obs = `r(N)' + 1
										qui set obs `obs'
										gsort `vdisp', mfirst										
										qui replace `vdisp' = `"`vlabel'"' if mi(`vdisp')
										qui replace `v_ord' = 0 if mi(`v_ord')
										gsort `v_ord'
									}
									
									if trim("`temp_p'")!="" {
										tempvar p
										qui g `p' = "`temp_p'" in 1
										if "`statlabel'"=="" {
											rtfrstyle `vdisp' `stat'* `p', cwidths(2500) local(b d e) cdadd("\clbrdrt\brdrs\clbrdrb\brdrs\clbrdrr\brdrs\clbrdrl\brdrs")
											if "`header'"!="" local headerrtf = "`b'\ql " + subinstr("`byline'","&","`d'\qc ",.)+ "`d'\qc " + "p value" + "`e'"
											listtab `vdisp' `stat'* `p', handle("`rtf'") begin("`b'\ql ") end("`e'") delim("`d'\qc ") head("`headerrtf'")
										}
										else {
											rtfrstyle `vdisp' `vstatlabel' `stat'* `p', cwidths(2500) local(b d e) cdadd("\clbrdrt\brdrs\clbrdrb\brdrs\clbrdrr\brdrs\clbrdrl\brdrs")
											if "`header'"!="" local headerrtf = "`b'\ql " + subinstr("`byline'","&","`d'\qc ",.) + "`d'\qc " + "p value" + "`e'"
											listtab `vdisp' `vstatlabel' `stat'* `p', handle("`rtf'") begin("`b'\ql ") end("`e'") delim("`d'\qc ") head("`headerrtf'")
										}
									}
									else if "`statlabel'"=="" {
										rtfrstyle `vdisp' `stat'*, cwidths(2500) local(b d e) cdadd("\clbrdrt\brdrs\clbrdrb\brdrs\clbrdrr\brdrs\clbrdrl\brdrs")
										if "`header'"!="" local headerrtf = "`b'\ql " + subinstr("`byline'","&","`d'\qc ",.) + "`e'"								
										listtab `vdisp' `stat'*, handle("`rtf'") begin("`b'\ql ") end("`e'") delim("`d'\qc ") head("`headerrtf'")
									}
									
									else {
										rtfrstyle `vdisp' `vstatlabel' `stat'*, cwidths(2500) local(b d e) cdadd("\clbrdrt\brdrs\clbrdrb\brdrs\clbrdrr\brdrs\clbrdrl\brdrs")
										if "`header'"!="" local headerrtf = "`b'\ql " + subinstr("`byline'","&","`d'\qc ",.) + "`e'"
										listtab `vdisp' `vstatlabel' `stat'*, handle("`rtf'") begin("`b'\ql ") end("`e'") delim("`d'\qc ") head("`headerrtf'")
									}
									
								}
								
								*returning results as locals
								foreach rown of numlist 1/`=_N' {
										local coln=0
										foreach returnv of varlist `stat'* {
												local ++coln
												return local `v'_`rown'_`coln'=`returnv'[`rown']
										}
										return local `v'_row`rown'=`v'[`rown']
										if trim("`labval'")!="" {
											local val`labval'`rown' : label (`v') `=`v'[`rown']'
											return local `v'_rowlbl`rown' = `"`val`labval'`rown''"'
										}
										
								}
								if trim("`temp_p'")!="" return local p="`temp_p'"
								
				} /* end of categorical variable*/
				
		} /*end of varlist loop*/
				
		*restoring data
		restore 
end






