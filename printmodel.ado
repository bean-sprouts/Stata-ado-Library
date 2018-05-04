capture program drop printmodel
capture program drop getcovars

program printmodel, rclass sortpreserve
version 13
syntax [anything] [, DECimal(integer 2) text HEADer non _cons saving(string asis) outsheet rtf(string)]

        tempname coefmatrix
        mat `coefmatrix'=r(table)

        * eform=0 means linear predictor, eform=1 means exponentiated
        local eform=`coefmatrix'[rownumb(matrix(`coefmatrix'),"eform"),1]
                if `eform'==0 local esttype="xb"
                else if `eform'==1 local esttype="exp"
                capture assert inlist(`eform',0,1)
                        if _rc>0 {
                                disp as error "eform from r(table) is a value other than 0 and 1.  Program needs to be updated to account for this!  I am not sure what it means!"
                                exit 1
                        }
        
        ***********************
        **  REGRESSION TYPE  **
        ***********************
                *creating local with command name*
                if trim("`e(cmd2)'")!=""  local cmd=trim("`e(cmd2)'")
                else local cmd=trim("`e(cmd)'")

                *if we are looking at a model done with imputed values, store the fact that this was an imputation estimation in a precommand local
                if "`cmd'"=="mi estimate" local cmd=trim("`e(cmd_mi)'")

        ********************
        **  N FROM MODEL  **
        ********************
                local modeln=`e(N)'
                        return local N=`modeln'
                        
        ******************************************************
        **  Creating variable list for covariates to print  **
        ******************************************************
                qui if "`anything'"!="" {
                        ds `anything'
                        local covarlist "`r(varlist)'"
                        *if constant requested, adding it from variable list
                        if "`_cons'"!="" local covarlist="`covarlist' _cons"
                }
                *if varlist not specified, then printing all vars in model
                else {
                        local covarlist: colnames `coefmatrix'
                        *if constant not requested, removing it from variable list
                        if "`_cons'"=="" local covarlist=reverse(subinstr(reverse("`covarlist'"),"snoc_","",1))
                        *printing error if xi: prefix not used for factor variables
                        foreach v in `covarlist' {
                                if index("`v'",".") & substr("`v'",1,2)!="o." {
                                        disp as err "`v': Factor variables and time-series operators not allowed.  For factor variables, use the xi: prefix."
                                        exit 1
                                }
                        }
                }
                

        *********************************
        **  saving results into table  **
        *********************************
        if `"`saving'"'!="" {
                tempfile out
                tempname memhold
                postfile `memhold' str50(label coef ci p)  using `out'
                
        }
        
        ****************************************************
        **  Print Options, OR, HR, XB, etc., and headers  **
        ****************************************************
                *putting out a header if requested
                if "`n'"=="" disp "N=`modeln'"
                if trim("`header'")!="" {
                        *tagging the type of beta estimate reported depending on the command
                        if trim("`esttype'")=="exp" {
                                if inlist("`cmd'","nbreg","poisson") local betatype="Incidence Rate Ratio"
                                else if inlist("`cmd'","stcox") local betatype="Hazard Ratio"
                                else if "`cmd'"=="stcrreg" local betatype="Subhazard Ratio"
                                else if inlist("`cmd'","logit","logistic","clogit") local betatype="Odds Ratio"
                                else local betatype="exp(beta)"
                        }       
                        *if estimated not exponentiated, estimate is a coefficnent, tag that
                        else local betatype="Coefficient"
                        if trim("`text'")!="" {
                                disp "Variable: `betatype' (95% Confidence Interval; p-value)"
                                local headerrtf "Variable: `betatype' (95% Confidence Interval; p-value)"
								return local header headerrtf
                        }
                        else {
                                disp "Variable&`betatype'&95% Confidence Interval&p-value"
                                local headerrtf "Variable&`betatype'&95% Confidence Interval&p-value"
								return local header headerrtf
                                if `"`saving'"'!="" post `memhold' ("Variable") ("`betatype'") ("95% Confidence Interval") ("p-value")
                        }
                }
                        
                
        ********************************************************
        **  GETTING COEFFICIENT ESTIMATES FOR EACH COVARIATE  **
        ********************************************************
				
				// Set up variables to save out results for RTF
				if trim("`rtf'")!="" {
					tempvar cov_rtf coef_rtf ci_rtf p_rtf
					qui g `cov_rtf' = ""
					qui g `coef_rtf' = ""
					qui g `ci_rtf' = ""
					qui g `p_rtf' = ""
				}
				local rtfnum = 0
		
                * Calculate estimates for each covariate
                local matcolnames: colnames `coefmatrix'
                foreach v in `covarlist' {
						// Add to "rtfnum" for each loop for RTF export. If not exporting, this is fine to run anyway as it doesn't do anything.
						local ++rtfnum
				
						*checking that the specifed variable is in the coef matrix
						local vinlist: list v in matcolnames
								if `vinlist'==0 {
										disp as err "`v' not a variable in model."
										exit 1
								}
						*getting the position of the varname in the matrix
						local covarn: list posof "`v'" in matcolnames
						 
						*getting beta coefficient, se, ub, lb, and p-value 
						if substr(word("`matcolnames'",`covarn'),1,2)!="o." {
								local `v'_beta=`coefmatrix'[rownumb(matrix(`coefmatrix'),"b"),`covarn']
								local se`v'=`coefmatrix'[rownumb(matrix(`coefmatrix'),"se"),`covarn']
										return local `v'_se=`se`v''
								local `v'_pexact=`coefmatrix'[rownumb(matrix(`coefmatrix'),"pvalue"),`covarn']
								local `v'_betalb=`coefmatrix'[rownumb(matrix(`coefmatrix'),"ll"),`covarn']
								local `v'_betaub=`coefmatrix'[rownumb(matrix(`coefmatrix'),"ul"),`covarn']

								*formatting p-value
								if trim("`text'")!="" qui pdisplay ``v'_pexact', local(`v'_pfmt) p
								else qui pdisplay ``v'_pexact', local(`v'_pfmt) 
										*returning p-values
										return local `v'_pexact="``v'_pexact'"
										return local `v'_pfmt="``v'_pfmt'"
								
								*rounding estimates
								foreach estimate in `v'_beta `v'_betalb `v'_betaub {
												*formating to display appropriate precision
												local `estimate'fmt=trim(string(``estimate'',"%9.`decimal'f"))
														*return results
														return local `estimate'="``estimate''"
														return local `estimate'fmt="``estimate'fmt'"
								}

								*creating var for CI
								local `v'_ci=trim("``v'_betalbfmt'")+", "+trim("``v'_betaubfmt'")
								return local `v'_ci=trim("``v'_ci'")
								
								*extract variable label
								if "`v'"=="_cons" local `v'_label="_cons"
								else {
										local `v'_label: var label `v'
										if trim("``v'_label'")=="" local `v'_label="`v'"
								}
								
								// For RTF export, save out variables for variable label/name, coefficient, CI, p value
								if trim("`rtf'")!="" {
								
									// If it's a categorical variable, include two extra spaces for variable label and reference group before first level
									if substr("`v'",1,2)=="_I" {
										
										// Original variable name
										local truename = substr("``v'_label'",1,index("``v'_label'","=")-1)
										
										// Variable label, if available
										local truelabel : var label `truename'
										if "`truelabel'"=="" local truelabel = "`truename'"
										
										// Variable value and value label, if available
										// For string, use actual value
										if upper(substr("`:type `truename''",1,3))=="STR" local trueval = substr("``v'_label'",index("``v'_label'","=")+2,.)
										
										// If not, use value label if available
										else {
											local ival = substr("`v'",-1,.)										
											local trueval : label (`truename') `ival'
										}
										
										// For the first time, save out a blank line with varname and another with reference group
										levelsof `cov_rtf', clean
										if strpos("`r(levels)'","`truelabel'")==0 {
											
											// Reference group label
											if upper(substr("`:type `truename''",1,3))=="STR" {
												preserve
												sort `truename'
												local refval = `truename'[1]
												restore
											}
											else local refval: label (`truename') `=`ival'-1'
											
											// Blank line with variable name
											qui replace `cov_rtf' = "`truelabel'" in `rtfnum'
											local ++rtfnum
											
											// Next line with reference group value/value label and p value
											qui replace `cov_rtf' = "    `refval'" in `rtfnum'
											qui replace `coef_rtf' = "Reference" in `rtfnum'
											qui replace `ci_rtf' = "-" in `rtfnum'
											local ++rtfnum
										}
										
										// For non-reference levels of categorical variable (including for first level)
											
										// Label with value/value label (if available) and save out coefficients
										qui replace `cov_rtf' = "    `trueval'" in `rtfnum'
										qui replace `coef_rtf' = "``v'_betafmt'" in `rtfnum'
										qui replace `ci_rtf' = "``v'_ci'" in `rtfnum'
										qui replace `p_rtf' = "``v'_pfmt'" in `rtfnum'											
										
									}
									
									// If not a categorical/factor variable, save out
									else {
										qui replace `cov_rtf' = "``v'_label'" in `rtfnum'
										qui replace `coef_rtf' = "``v'_betafmt'" in `rtfnum'
										qui replace `ci_rtf' = "``v'_ci'" in `rtfnum'
										qui replace `p_rtf' = "``v'_pfmt'" in `rtfnum'
									}
									
								}

						}
						
						*for omitted variables print NA
						else if index("`v'","o.") {
								*assigning variable label
								local var_label: var label `=subinstr("`v'","o.","",1)'
								if trim("`var_label'")=="" local var_label="`v'"
								
								local v=subinstr("`v'","o.","o_",1)
								
								local `v'_label="`var_label'"
								local `v'_beta="NA"
								local se`v'="NA"
								local `v'_pexact="NA"
								local `v'_betalb="NA"
								local `v'_betaub="NA"
								local `v'_pfmt="NA"
								local `v'_ci="NA"
								
								return local `v'_pexact="NA"
								return local `v'_pfmt="NA"
								return local `v'_ci="NA"
								foreach estimate in `v'_beta `v'_betalb `v'_betaub {
												*formatting to display appropriate precision
												local `estimate'fmt="NA"
														*return results
														return local `estimate'="NA"
														return local `estimate'fmt="NA"
								}
						}       
						*this loop should never be reached.  if it does then there is situation we did not anticiapte
						else {
								print as err "`v': Not sure how to deal with this variable."
								exit 1
						}
								
						*returning variable label
						return local `v'_label="`v'"
																		
						/*displaying result to be copied into text of manuscript*/
						if trim("`text'")!="" {
										disp  trim("``v'_label'") " " trim("``v'_betafmt'") " (95% CI " trim("``v'_ci'") "; " trim("``v'_pfmt'") ")"
										return local `v'_stats=trim("``v'_betafmt'")+" (95% CI "+trim("``v'_ci'")+"; "+trim("``v'_pfmt'")+")"
										return local `v'_all=trim("``v'_label'")+" "+trim("``v'_betafmt'")+"(95% CI "+trim("``v'_ci'")+"; "+trim("``v'_pfmt'")+")"
						}
						/*displaying result to be copied into table*/
						else {
										disp trim("``v'_label'") "&" trim("``v'_betafmt'") "&" trim("``v'_ci'") "&" trim("``v'_pfmt'") 
										return local `v'_stats=trim("``v'_betafmt'")+"&"+trim("``v'_ci'")+"&"+trim("``v'_pfmt'") 
										return local `v'_all=trim("``v'_label'")+"&"+trim("``v'_betafmt'")+"&"+trim("``v'_ci'")+"&"+trim("``v'_pfmt'") 
										if `"`saving'"'!="" post `memhold' (`"`=trim("``v'_label'")'"') (`"`=trim("``v'_betafmt'")'"') (`"`=trim("``v'_ci'")'"') (`"`=trim("``v'_pfmt'")'"')
						}
				
				} /*end covar loop*/
				
				*
				if `"`saving'"'!="" {
						preserve
						
						postclose `memhold'
						use `out', clear
						qui compress
						
						if "`outsheet'"!="" outsheet using `saving'
						else save `saving'
						
						restore
				}


        *return matrix of coefficients
        return matrix table=`coefmatrix'
				
		// Save out to RTF file
		if trim("`rtf'")!="" {
			file write `rtf' "{\pard N=`modeln'\par}"
			rtfrstyle `cov_rtf' `coef_rtf' `ci_rtf' `p_rtf', cwidths(2500) local(b d e) cdadd("\qc\clbrdrt\brdrs\clbrdrb\brdrs\clbrdrr\brdrs\clbrdrl\brdrs")
			local headerrtf = "`b'\ql " + subinstr("`headerrtf'","&","`d'\qc ",.) + "`e'"
			if "`header'"=="" local headerrtf = ""
			listtab `cov_rtf' `coef_rtf' `ci_rtf' `p_rtf' if !mi(`cov_rtf'), handle("`rtf'") begin("`b'\ql ") delim("`d'\qc ") end("`e'") head("`headerrtf'")
            file write `rtf' "{\pard \par}"
		}

end

program getcovars, rclass
syntax anything [if] [in] [pweight fweight aweight] [, *]

        return local varlist="`anything'"
                return local modeloptions="`options'"
end

