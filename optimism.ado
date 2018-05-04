/*
PROGRAM: optimism.ado
PROGRAMMER: Melissa + Dan
DATE: 1/26/2015
DESCRIPTION: This program generates Optimism-corrected statistics using bootstrap methods
*/


//OPTIMISM PROGRAM
capture program drop optimism
program optimism, rclass

	syntax anything(equalok) [if] [in], train(string) apply(string)[reps(integer 200) seed(integer -1) strata(passthru) DECimal(integer 2) saving(string asis)]

qui {
	
	*keeping observations for analyses
	tempfile dataorig
	save `dataorig'
	capture keep `if' `in'
	
	*if seed requested setting it now.
	if `seed'>=0 set 
	
	*seperating commands into different locals
	local traini=0
	while length(trim("`train'"))>0 {
		local ++traini
		if index("`train'",";") local train`traini'=substr("`train'",1,index("`train'",";")-1)
		else local train`traini'="`train'"
		
		local train=subinstr("`train'","`train`traini''","",1)
		
		*if first character is a semicolon deleteing it
		if substr("`train'",1,1)==";" local train=subinstr("`train'",";","",1)
		disp "`train`traini''"
	}
	local applyi=0
	while length(trim("`apply'"))>0 {
		local ++applyi
		if index("`apply'",";") local apply`applyi'=substr("`apply'",1,index("`apply'",";")-1)
		else local apply`applyi'="`apply'"
		
		local apply=subinstr("`apply'","`apply`applyi''","",1)
		
		*if first character is a semicolon deleteing it
		if substr("`apply'",1,1)==";" local apply=subinstr("`apply'",";","",1)
		disp "`apply`applyi''"
		
	}
			
			
	//Return local variables
	//We wish to return variables that specify their labels and if there is no label the order in which the variable is listed will be used
	//Separate the label from the statistic if there is an = sign
	//this section also generates a variable called stat`numvars' which will be used in the return statement
	*generate a new local which includes on the statistical commands portion of anything		
	local runstats=""
	local numvars=0
	foreach statistic in `anything' {
		local ++numvars
								
		*identifying those with an = sign
		if index("`statistic'","=")!=0 {
			*create the label for the returned locals
			local lab`numvars'=word(subinstr("`statistic'","="," ",.),1)
			*create the label for the table
			local tablelab`numvars'=word(subinstr("`statistic'","="," ",.),1)
			*create the statistic of interest
			local stat`numvars'=word(subinstr("`statistic'","="," ",.),2)
		}              
		else {
			*create the statistic of interest
			local stat`numvars'="`statistic'"
			*create the label for those without a label for the local list
			local lab`numvars'="stat`numvars'"
			*create the label for the table for those without a label
			local tablelab`numvars'="`statistic'"
		}
		local runstats="`runstats'" + " `stat`numvars''"
	}
		


	*estimating optimism through bootstraps
	tempname memhold
	tempfile results
	postfile `memhold' boot str20(dataset statname) stat using `results'

	foreach boot of numlist 1/`reps' {
		disp "boot `boot'"
		use `dataorig', clear
		capture keep `if' `in'

		bsample, `strata'
	
		*obtaining estimate on bootstrap sampled dataset
			qui foreach i of numlist 1/`traini' {
				`train`i''
			}
			qui foreach i of numlist 1/`applyi' {
				`apply`i''
			}		
			
		*saving out stats estimates from bootstrapped dataset
		foreach statistic in `runstats' {
			post `memhold' (`boot') ("bootstrapped") ("`statistic'") (``statistic'')
		}
		
		*applying to the original dataset
			use `dataorig', clear
			capture keep `if' `in'
			
			disp "`applyi'"
			foreach i of numlist 1/`applyi' {
				count
				disp "i=`i'"
				disp "`apply`i''"
				`apply`i''
			}	
		
		
		*saving out stats estimates from original dataset
		foreach statistic in `runstats' {
			post `memhold' (`boot') ("original") ("`statistic'") (``statistic'')
		}
		
		
		*Calculating the Apparent Value
			use `dataorig', clear
			capture keep `if' `in'
					
			*obtaining estimate apparent (uncorrected)
			qui foreach i of numlist 1/`traini' {
				`train`i''
			}
			qui foreach i of numlist 1/`applyi' {
				`apply`i''
			}
			
		*saving out stats estimates from original dataset
		foreach statistic in `runstats' {
			post `memhold' (`boot') ("apparent") ("`statistic'") (``statistic'')
		}	

	}
	postclose `memhold'
	use `results', clear
	
	*saving out output of results if requested
	if "`saving'"!="" {
		save `saving'
	}
	
	
	//reshaping data and outputting results
	*so we needed to drop the duplicates in order to reshape
		duplicates drop
		*reshape to wide format
		reshape wide stat, i(boot statname) j(dataset) string
		
		
	//Format results based on format option
		*setting format for AUC and CI
		local format="%9.`decimal'f"
		
			
	*Generate the n-statistic set for each statistic requested
	local numvars=0
	foreach statistic in `runstats' {
		local ++numvars
		preserve
		
		*keeping only results pertinant to the statistic
		keep if statname=="`statistic'"

		*mean c-index of the 200 bootstrapped models on the original dataset
		sum statoriginal
		g original=`r(mean)'
		local original=string(`r(mean)',"`format'")
		
		*calculating the optimism
		g opt = statbootstrapped - statoriginal
		sum opt 
		g optimism=`r(mean)'
		local optimism=string(`r(mean)',"`format'")
		
		*mean c-index of the 200 bootstrapped models on the bootsrapped dataset
		sum statbootstrapped 
		g bootstrapped=`r(mean)'
		local bootstrapped=string(`r(mean)',"`format'")
	
		*apparent value of the statistic
		sum statapparent
		assert `r(sd)'==0
		g apparent=`r(mean)'
		local apparent=string(`r(mean)',"`format'")
		
		*Calculate the bootstrap-corrected value
		g bc = apparent - optimism
		sum bc
		local bscorrected=string(`r(mean)',"`format'")
		
		//Outputting apparent statistic, optimism, and corrected statistic
		noi disp "`tablelab`numvars''"
	
		*outputting using locals
		noi disp "Apparent Value&Optimism&Bootstrapped-Corrected Value"
		noi disp "`apparent'&`optimism'&`bscorrected'"
		

	//Return values as locals	
	return local optimism`lab`numvars''="`optimism'"
	return local apparent`lab`numvars''="`apparent'"
	return local bscorrected`lab`numvars''="`bscorrected'"
	return local avoriginal`lab`numvars''="`original'"
	return local avbootstrapped`lab`numvars''="`bootstrapped'"
	
	restore
	}	

	//Return original dataset
	use `dataorig', clear
}

end

