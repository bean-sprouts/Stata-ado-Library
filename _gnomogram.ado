

/*
THE rotterdam FUNCTION GENERATES A NEW VARIABLE 
THAT IS THE PREDICTION FOR PCa ON BX BASED ON THE 
ERSPC ROTTERDAM COHORT
PROGRAMMER: AM
DATE: 12/20/2012
*/

*
capture drop nomorisk
capture program drop _gnomogram
program define _gnomogram
	version 12

	gettoken type 0 : 0
	gettoken g 0 : 0
	gettoken eqs 0 : 0
	gettoken lparen 0 : 0, parse("(")
	gettoken rparen 0 : 0, parse(")")
	
	local cmd="`=substr(trim("`0'"),2,.)'"
	
	
	if inlist(trim(upper("`rparen'")),"BX VICKERS 2008","ROTTERDAM") {
		disp as result _newline "`rparen'"
		rotterdam `cmd'
		g `type' `g' = nomogram_tempvarname798465
		drop nomogram_tempvarname798465	
		}
	else if inlist(trim(upper("`rparen'")),"ROTTERDAMHG") {
		disp as result _newline "`rparen'"
		rotterdamhg `cmd'
		g `type' `g' = nomogram_tempvarname798465
		drop nomogram_tempvarname798465	
		}
	else if trim(upper("`rparen'"))=="BCR STEPHENSON 2005" {
		disp as result _newline "`rparen'"
		stephenson `cmd'
		g `type' `g' = nomogram_tempvarname798465
		drop nomogram_tempvarname798465	
		}		
	else if trim(upper("`rparen'"))=="KATTAN" {
		disp as result _newline "`rparen'"
		kattan `cmd'
		g `type' `g' = nomogram_tempvarname798465
		drop nomogram_tempvarname798465	
		}
	else if trim(upper("`rparen'"))=="PROTECT" {
		disp as result _newline "`rparen'"
		protect `cmd'
		g `type' `g' = nomogram_tempvarname798465
		drop nomogram_tempvarname798465	
		}
	else if trim(upper("`rparen'"))=="PROTECTHG" {
		disp as result _newline "`rparen'"
		protecthg `cmd'
		g `type' `g' = nomogram_tempvarname798465
		drop nomogram_tempvarname798465	
		}
	else if trim(upper("`rparen'"))=="MSKCC PREOP BCR" {
		disp as result _newline "`rparen'"
		preopbcr `cmd'
		g `type' `g' = nomogram_tempvarname798465
		drop nomogram_tempvarname798465	
		}
	else if trim(upper("`rparen'"))=="MSKCC CORES PREOP BCR" {
		disp as result _newline "`rparen'"
		preopbcrcores `cmd'
		g `type' `g' = nomogram_tempvarname798465
		drop nomogram_tempvarname798465	
		}	
	else if trim(upper("`rparen'"))=="MSKCC POSTOP BCR" {
		disp as result _newline "`rparen'"
		postopbcr `cmd'
		g `type' `g' = nomogram_tempvarname798465
		drop nomogram_tempvarname798465	
		}
	* if rparen is empty - list out possible nomogram choices
	else {
		disp as result _newline "Nomogram `rparen' does not exist" _newline
		disp "Available Nomograms Include: " _newline _newline _column(5) "BX VICKERS 2008" _newline "BCR STEPHENSON 2005" /// 
		_newline  "KATTAN" _newline  "PROTECT" _newline  "PROTECTHG" ///
		_newline  "MSKCC PREOP BCR" _newline  "MSKCC CORES PREOP BCR" _newline  "MSKCC POSTOP BCR" 
		error  
		}
		
end
	
capture program drop protect
program define protect
	version 12
	
	syntax [if] [in], [age(varname) tpsa(varname) fpsa(varname) ipsa(varname) hk2(varname) dreneg(varname) drepos(varname) xb g(string)]
	marksample touse	
	
	disp as result "PROTECT ANYGRADE NOMOGRAM"

	qui foreach v in age tpsa fpsa ipsa hk2 {
		if trim("``v''") == "" {
			if "`v'"=="age" noi disp as text "`v' assumed to be variable `v' in years"
			else if inlist("`v'","tpsa","fpsa","ipsa","hk2") noi disp "`v' assumed to be variable `v' in ng/ml"
		
			*marking the observations to exclude due to missing values
			replace `touse'=. if mi(`v')
			
			local `v'="`v'"
		}
	}

	
	qui {
		tempvar sptpsa1 sptpsa2 spfpsa1 spfpsa2 x1 x2 x3 x4 drenegvar dreposvar xb_full phat_full
		
		*Outputting error message if the clinical stage does not satisfy the inclusion criteria
		capture assert `tpsa'>=0 if `touse' 
		if _rc>0 {
			disp as err "tPSA must be >= 0"
			exit
		}
		
		capture assert inrange(`age',0,130) if `touse' & !mi(`age')
		if _rc>0 {
			disp as err "Age must be between 0 and 130 years"
			exit
		}

		*output a message to notify coders of dre specifications
		if (trim("`drepos'")!="" & trim("`dreneg'")=="") | (trim("`dreneg'")!="" & trim("`drepos'")=="") {
			disp as err "Specify neither or both dreneg and drepos"
			exit
		}
		if trim("`dreneg'")=="" {
			g `drenegvar'=0
			disp as err "Specify dreneg to include DRE in prediction model"
		}
		else g `drenegvar'=`dreneg'
		
		if trim("`drepos'")=="" {
			g `dreposvar'=0
			disp as err "Specify drepos to include DRE in prediction model"
		}
		else g `dreposvar'=`drepos'
		
		
		g `sptpsa1'=max(`tpsa'-0,0)^3 - ((265-0)/(265-5.61))*max(`tpsa'-5.61,0)^3 + ((5.61-0)/(265-5.61))*max(`tpsa'-265,0)^3
		g `sptpsa2'=max(`tpsa'-3.92,0)^3 - ((265-3.92)/(265-5.61))*max(`tpsa'-5.61,0)^3 + ((5.61-3.92)/(265-5.61))*max(`tpsa'-265,0)^3
		g `spfpsa1'=max(`fpsa'-0,0)^3 - ((24-0)/(24-1.21))*max(`fpsa'-1.21,0)^3 + ((1.21-0)/(24-1.21))*max(`fpsa'-24,0)^3
		g `spfpsa2'=max(`fpsa'-0.82,0)^3 - ((24-0.82)/(24-1.21))*max(`fpsa'-1.21,0)^3 + ((1.21-0.82)/(24-1.21))*max(`fpsa'-24,0)^3


		g `xb_full'=-2.4265 + 0.0339*`age' + 0.4001*`tpsa' + -0.0002*`sptpsa1' + ///
					  -0.0011*`sptpsa2' + -3.4534*`fpsa' + 0.6256*`spfpsa1' + ///
					  -1.7992*`spfpsa2' + 0.0495*`ipsa' + 9.5970*`hk2'+-0.2*`drenegvar'+0.8*`dreposvar' if `tpsa'<=25
		
		replace `xb_full'=-1.1200 + 0.1118*`tpsa' if `tpsa'>25
		
		g `phat_full'=  exp(`xb_full')/(1+ exp(`xb_full') )

	}
		
		/*NEW VARIABLE IS EQUAL TO THE FIRST NONMISSING VARIABLE IN THE VARLIST */
		if trim("`xb'")=="" {
			g nomogram_tempvarname798465=`phat_full' if `touse'
			*g ${EGEN_Varname} =`phat_full' if `touse'
			*label var ${EGEN_Varname} "Risk of PCa on Bx from ProtecT Model"
		}
		else if trim("`xb'")!="" {
			g nomogram_tempvarname798465 =`xb_full' if `touse'
			*label var ${EGEN_Varname} "Risk of PCa on Bx from ProtecT Model (Linear Prediction)"
		}
			
	
	* list out the paper info: title, first author, journal and PMID 
	
end

capture program drop protecthg
program define protecthg
	version 12
	
	syntax [if] [in], [age(varname) tpsa(varname) fpsa(varname) ipsa(varname) hk2(varname) dreneg(varname) drepos(varname) xb g(string)]
	marksample touse
		
	disp as result "PROTECT HIGHGRADE NOMOGRAM"

	qui foreach v in age tpsa fpsa ipsa hk2 {
		if trim("``v''") == "" {
			if "`v'"=="age" noi disp as text "`v' assumed to be varaible `v' in years"
			else if inlist("`v'","tpsa","fpsa","ipsa","hk2") noi disp "`v' assumed to be varaible `v' in ng/ml"
		
			*marking the observations to exclude due to missing values
			replace `touse'=. if mi(`v')
		
			local `v'="`v'"
		}
	}
	
	qui {
		tempvar sptpsa1 sptpsa2 spfpsa1 spfpsa2 x1 x2 x3 x4 drenegvar dreposvar xb_full phat_full 
		
		*Outputting error message if the clinical stage does not satisfy the inclusion criteria
		capture assert `tpsa'>=0 if `touse' 
		if _rc>0 {
			disp as err "tPSA must be >= 0"
			exit
		}
		
		capture assert inrange(`age',0,130) if `touse' & !mi(`age')
		if _rc>0 {
			disp as err "Age must be between 0 and 130 years"
			exit
		}
		
		/*if trim("`dreneg'")=="" g `drenegvar'=0
		else g `drenegvar'=`dreneg'
		if trim("`drepos'")=="" g `dreposvar'=0
		else g `dreposvar'=`drepos'*/	
		*output a message to notify coders of dre specifications
		if (trim("`drepos'")!="" & trim("`dreneg'")=="") | (trim("`dreneg'")!="" & trim("`drepos'")=="") {
			disp as err "Specify neither or both dreneg and drepos"
			exit
		}
		if trim("`dreneg'")=="" {
			g `drenegvar'=0
			disp as err "Specify dreneg to include DRE in prediction model"
		}
		else g `drenegvar'=`dreneg'
		
		if trim("`drepos'")=="" {
			g `dreposvar'=0
			disp as err "Specify drepos to include DRE in prediction model"
		}
		else g `dreposvar'=`drepos'
		
		
		g `sptpsa1'=max(`tpsa'-0,0)^3 - ((265-0)/(265-5.61))*max(`tpsa'-5.61,0)^3 + ((5.61-0)/(265-5.61))*max(`tpsa'-265,0)^3
		g `sptpsa2'=max(`tpsa'-3.92,0)^3 - ((265-3.92)/(265-5.61))*max(`tpsa'-5.61,0)^3 + ((5.61-3.92)/(265-5.61))*max(`tpsa'-265,0)^3
		g `spfpsa1'=max(`fpsa'-0,0)^3 - ((24-0)/(24-1.21))*max(`fpsa'-1.21,0)^3 + ((1.21-0)/(24-1.21))*max(`fpsa'-24,0)^3
		g `spfpsa2'=max(`fpsa'-0.82,0)^3 - ((24-0.82)/(24-1.21))*max(`fpsa'-1.21,0)^3 + ((1.21-0.82)/(24-1.21))*max(`fpsa'-24,0)^3


		g `xb_full'=-6.4469 + 0.0570*`age' + 0.8755*`tpsa' + -0.0054*`sptpsa1' + 0.0144*`sptpsa2' + ///
					-5.4028*`fpsa' + 0.8872*`spfpsa1' + -2.5798*`spfpsa2' + ///
					2.1820*`ipsa' + 6.9644*`hk2'+-0.3*`drenegvar'+1.2*`dreposvar' if `tpsa'<=25
		
		replace `xb_full'=0.8701 + 0.0136*`tpsa' if `tpsa'>25
		
		g `phat_full'=  exp(`xb_full')/(1+ exp(`xb_full') )

	}
		
		/*NEW VARIABLE IS EQUAL TO THE FIRST NONMISSING VARIABLE IN THE VARLIST */
		if trim("`xb'")=="" {
			g nomogram_tempvarname798465=`phat_full' if `touse'
			*g ${EGEN_Varname} =`phat_full' if `touse'
			*label var ${EGEN_Varname} "Risk of High Grade PCa on Bx from ProtecT Model"
		}
		else if trim("`xb'")!="" {
			g nomogram_tempvarname798465 =`xb_full' if `touse'
			*label var ${EGEN_Varname} "Risk of High Grade PCa on Bx from ProtecT Model (Linear Prediction)"
		}
			
	
	* list out the paper info: title, first author, journal and PMID 
	
end


capture program drop rotterdam
program define rotterdam
	version 12
	
	syntax [if] [in], [age(varname) tpsa(varname) fpsa(varname) ipsa(varname) hk2(varname) dre(varname) xb g(string)]
	marksample touse
	
	disp as result "ROTTERDAM NOMOGRAM"

	qui foreach v in age tpsa fpsa ipsa hk2 {
		if trim("``v''") == "" {
			if "`v'"=="age" noi disp as text "`v' assumed to be varaible `v' in years"
			else if inlist("`v'","tpsa","fpsa","ipsa","hk2") noi disp "`v' assumed to be varaible `v' in ng/ml"
		
			*marking the observations to exclude due to missing values
			replace `touse'=. if mi(`v')
		
			local `v'="`v'"
		}
	}
	
	
	qui {
		tempvar sptpsa1 sptpsa2 spfpsa1 spfpsa2 x1 x2 x3 x4 xb_full phat_full 
		
		*Outputting error message if the clinical stage does not satisfy the inclusion criteria
		capture assert `tpsa'>=0 if `touse' 
		if _rc>0 {
			disp as err "tPSA must be >= 0"
			exit
		}
		
		capture assert inrange(`age',0,130) if `touse' & !mi(`age')
		if _rc>0 {
			disp as err "Age must be between 0 and 130 years"
			exit
		}
		
		gen `sptpsa1' = -(162 - 4.4503)/(162 - 3)*(`tpsa'-3)^3 + max(`tpsa'-4.4503, 0)^3	
		gen `sptpsa2' = -(162 - 6.4406)/(162 - 3)*(`tpsa'-3)^3 + max(`tpsa'-6.4406, 0)^3
		 
		gen `spfpsa1' = -(11.8 - 0.84)/(11.8 - 0.25)*(`fpsa'-0.25)^3 + max(`fpsa'-0.84, 0)^3	
		replace `spfpsa1' = (11.8 - 0.84)*(0.84 - 0.25)*(11.8 + 0.84 + 0.25 - 3*`fpsa') if `fpsa'>11.8

		gen `spfpsa2' = -(11.8 - 1.29)/(11.8 - 0.25)*(`fpsa' - 0.25)^3 + max(`fpsa' - 1.29, 0)^3
		replace `spfpsa2' = (11.8 - 1.29)*(1.29 - 0.25)*(11.8 + 1.29 + 0.25 - 3*`fpsa') if `fpsa'>11.8

		**  WITHOUT DRE  **		
		if trim("`dre'")=="" {
			disp as err "Specify DRE variable to include in model"
			/* DEFINE THE VARIOUS INPUTS TO THE MODEL */
			g `x1' = 0.0846726*`tpsa' - 0.0211959*`sptpsa1' + 0.0092731*`sptpsa2'
			g `x2' = -3.717517*`fpsa' - 0.6000171*`spfpsa1' + 0.275367 *`spfpsa2'
			g `x3' = 3.968052*`ipsa' 
			g `x4' = 4.508231*`hk2' 

			/* AND CALCULATE LINEAR PREDICTION */
			g `xb_full' = -1.735529 + 0.0172287*`age' + `x1' + `x2' + `x3' + `x4'
		}
		** WITH DRE  **
		else {
			g `x1' = 0.0637121*`tpsa' - 0.0199247*`sptpsa1' + 0.0087081*`sptpsa2'
			g `x2' = -3.460508*`fpsa' - 0.4361686*`spfpsa1' + 0.1801519*`spfpsa2'
			g `x3' = 4.014925*`ipsa' 
			g `x4' = 3.523849*`hk2' 

			/* AND CALCULATE LINEAR PREDICTION */
			g `xb_full' = -1.373544 + 0.9661025*`dre' + 0.0070077*`age' + `x1' + `x2' + `x3' + `x4'
		}
		* and for those with PSA>25, seperate calculation
		replace `xb_full' = 0.0733628*`tpsa' - 1.377984 if `tpsa'>25

		* and get prediction
		g `phat_full' = exp(`xb_full')/(1 + exp(`xb_full'))

	}
		
		/*NEW VARIABLE IS EQUAL TO THE FIRST NONMISSING VARIABLE IN THE VARLIST */
		if trim("`xb'")=="" {
			g nomogram_tempvarname798465=`phat_full' if `touse'
			*g ${EGEN_Varname} =`phat_full' if `touse'
			*label var ${EGEN_Varname} "Risk of PCa on Bx from Rotterdam Model"
		}
		else if trim("`xb'")!="" {
			g nomogram_tempvarname798465=`xb_full' if `touse'
			*label var ${EGEN_Varname} "Risk of PCa on Bx from Rotterdam Model (Linear Prediction)"
		}
			
	
	* list out the paper info: title, first author, journal and PMID 
	
end

capture program drop rotterdamhg
program define rotterdamhg
	version 12
	
	syntax [if] [in], [age(varname) tpsa(varname) fpsa(varname) ipsa(varname) hk2(varname) dre(varname) xb g(string)]
	marksample touse
	
	disp as result "ROTTERDAM NOMOGRAM (high grade)"

	qui foreach v in age tpsa fpsa ipsa hk2 {
		if trim("``v''") == "" {
			if "`v'"=="age" noi disp as text "`v' assumed to be varaible `v' in years"
			else if inlist("`v'","tpsa","fpsa","ipsa","hk2") noi disp "`v' assumed to be varaible `v' in ng/ml"
		
			*marking the observations to exclude due to missing values
			replace `touse'=. if mi(`v')
		
			local `v'="`v'"
		}
	}
	
		
	qui {
		tempvar sptpsa1 sptpsa2 spfpsa1 spfpsa2 xb_full phat_full
		
		*Outputting error message if the clinical stage does not satisfy the inclusion criteria
		capture assert `tpsa'>=0 if `touse' 
		if _rc>0 {
			disp as err "tPSA must be >= 0"
			exit
		}
		
		capture assert inrange(`age',0,130) if `touse' & !mi(`age')
		if _rc>0 {
			disp as err "Age must be between 0 and 130 years"
			exit
		}
		
		local k0tpsa=0.0100
		local kNtpsa=245.0000 
		local k1tpsa=4.2400 
		local k2tpsa=6.2000
		local k0fpsa=0.0000
		local kNfpsa=17.3000 
		local k1fpsa=0.8300 
		local k2fpsa=1.2700
		
		gen `sptpsa1' = -((scalar(`kNtpsa') - scalar(`k1tpsa')) / (scalar(`kNtpsa') - scalar(`k0tpsa'))) * (`tpsa' - scalar(`k0tpsa'))^3 
			replace `sptpsa1' = `sptpsa1' + (`tpsa' - scalar(`k1tpsa'))^3 if `tpsa' > scalar(`k1tpsa') & scalar(`k1tpsa') > scalar(`k0tpsa')

		gen `sptpsa2' = -((scalar(`kNtpsa') - scalar(`k2tpsa')) / (scalar(`kNtpsa') - scalar(`k0tpsa'))) * (`tpsa' - scalar(`k0tpsa'))^3 
			replace `sptpsa2' = `sptpsa2' + (`tpsa' - scalar(`k2tpsa'))^3 if `tpsa' > scalar(`k2tpsa') & scalar(`k2tpsa') > scalar(`k0tpsa')

		
		gen `spfpsa1' = -((scalar(`kNfpsa') - scalar(`k1fpsa')) / (scalar(`kNfpsa') - scalar(`k0fpsa'))) * (`fpsa' - scalar(`k0fpsa'))^3 
			replace `spfpsa1' = `spfpsa1' + (`fpsa' - scalar(`k1fpsa'))^3 if `fpsa' > scalar(`k1fpsa') & scalar(`k1fpsa') > scalar(`k0fpsa')

		gen `spfpsa2' = -((scalar(`kNfpsa') - scalar(`k2fpsa')) / (scalar(`kNfpsa') - scalar(`k0fpsa'))) * (`fpsa' - scalar(`k0fpsa'))^3 
			replace `spfpsa2' = `spfpsa2' + (`fpsa' - scalar(`k2fpsa'))^3 if `fpsa' > scalar(`k2fpsa') & scalar(`k2fpsa') > scalar(`k0fpsa')
			
		***********************
		****  WITHOUT DRE  ****
		***********************
		/* AND CALCULATE LINEAR PREDICTION */
		if trim("`dre'")=="" {
				disp as err "Specify DRE variable to include in model"
				g `xb_full' = -9.114314686532467 + `age'*.0826815 + `tpsa'*.4849174 + `sptpsa1'*-.0051319 + `sptpsa2'*.0044963 + ///
							`fpsa'*-3.760412 + `spfpsa1'*-.3237043 + `spfpsa2'*.1404233 + `ipsa'*2.538381 + `hk2'*4.114796
		}
							
		********************
		****  WITH DRE  ****
		********************
		else g `xb_full' =-9.003410369690776 + `age'*.0696931 + `tpsa'*.4852348 + `sptpsa1'*-.003319 + `sptpsa2'*.0032038 + ///
							`fpsa'*-3.173072 + `spfpsa1'*-.0727489 + `spfpsa2'*-.0171745 + `ipsa'*2.265643 + `hk2'*4.117794 + `dre'*1.308971

							
		* and for those with PSA>25, seperate calculation
		replace `xb_full' = -.9755556108334966 + `tpsa'*.0373936 if `tpsa'>25



		
		* and get prediction
		g `phat_full' = exp(`xb_full')/(1 + exp(`xb_full'))

	}
		
		/*NEW VARIABLE IS EQUAL TO THE FIRST NONMISSING VARIABLE IN THE VARLIST */
		if trim("`xb'")=="" {
			g nomogram_tempvarname798465=`phat_full' if `touse'
			*g ${EGEN_Varname} =`phat_full' if `touse'
			*label var ${EGEN_Varname} "Risk of High Grade PCa on Bx from Rotterdam Model"
		}
		else if trim("`xb'")!="" {
			g nomogram_tempvarname798465=`xb_full' if `touse'
			*label var ${EGEN_Varname} "Risk of High Grade PCa on Bx from Rotterdam Model (Linear Prediction)"
		}
			
	
	* list out the paper info: title, first author, journal and PMID 
	
	/* AS this is an unpublished model, the code to create it is below.
			cd "O:\Outcomes\Andrew\Analytic Projects\1 Active\Lilja ERSPC 2005\_Submitted or Inactive\_Rotterdam\paper 1 replication of Goteborg round 1 biopsy"	
			use "training round 1 replication.dta", clear
			append using "validation round 1 replication.dta", gen(vali)

			* spline all markers
			local markers="tpsa fpsa"
			foreach v of local markers {

				* save out the location of knots and minimum and maximum values in the training set so that they
				* can be used when computing splines in the validation set
				quietly centile `v', c(33 66)
				local k1`v' = r(c_1)
				local k2`v' = r(c_2)
					quietly sum `v'
					local k0`v' = r(min)
					local kN`v' = r(max)
				disp "`v'"  "&" %9.4f `k0`v''  "&" %9.4f `kN`v'' "&"  %9.4f `k1`v'' "&"  %9.4f `k2`v''

						gen sp`v'1 = -((scalar(`kN`v'') - scalar(`k1`v'')) / (scalar(`kN`v'') - scalar(`k0`v''))) * (`v' - scalar(`k0`v''))^3 
						replace sp`v'1 = sp`v'1 + (`v' - scalar(`k1`v''))^3 if `v' > scalar(`k1`v'') & scalar(`k1`v'') > scalar(`k0`v'')
				
						gen sp`v'2 = -((scalar(`kN`v'') - scalar(`k2`v'')) / (scalar(`kN`v'') - scalar(`k0`v''))) * (`v' - scalar(`k0`v''))^3 
						replace sp`v'2 = sp`v'2 + (`v' - scalar(`k2`v''))^3 if `v' > scalar(`k2`v'') & scalar(`k2`v'') > scalar(`k0`v'')

			}
			
			logit higrade age tpsa sptpsa1 sptpsa2 fpsa spfpsa1 spfpsa2 ipsa hk2 if tpsa<25, nolog
			logit higrade age tpsa sptpsa1 sptpsa2 fpsa spfpsa1 spfpsa2 ipsa hk2 dre if tpsa<25, nolog
			logit higrade tpsa if tpsa>25
*/
	
end


capture program drop stephenson 
program define stephenson 
	version 12
	
	syntax [if] [in], [tpsa(varname) pathgg1(varname) pathgg2(varname) yos(varname) ece(varname) svi(varname) lni(varname) sms(varname) xb g(string)]
	marksample touse
	
	disp as result "STEPHENSON NOMOGRAM"

	qui foreach v in tpsa pathgg1 pathgg2 yos ece svi lni sms {
	if trim("``v''") == "" {        
			if "`v'"=="yos" noi disp as text "`v' assumed to be varaible `v' in years"
			else if inlist("`v'","pretxpsa") noi disp as text "`v' assumed to be varaible `v' in ng/mL"
			else if inlist("`v'","pathgg1", "pathgg2") noi disp as text "`v' assumed to be varaible `v' in numeric format 0-5"
			else if inlist("`v'","ece", "svi", "lni", "sms") noi disp as text "`v' assumed to be varaible `v' in binary format"
		
			*marking the observations to exclude due to missing values
			replace `touse'=. if mi(`v')
		
			local `v'="`v'"
		}
	}


/* 		* 5-YEAR NOMOGRAM PROBABILITY
		Variables required for the calculation are:
			Pretreatment PSA													* tpsa  
			Primary and secondary Gleason scores from pathology specimen        * pathgg1 & pathgg2
			Year of radical prostatectomy                                       * yos  
			Extracapsular extension (yes or no)                                 * ece
			Seminal vesicle invasion (yes or no)                                * svi
			Lymph node involvement (yes or no)                                  * lni
			Positive surgical margins (yes or no)                               * sms
			
		* NOTE: The nomogram is valid only for patients who did not receive any neoadjuvant therapy.	
		
		*/    
		
	qui {
		tempvar xb_full phat_full badtpsa
		
		*Outputting error message if the clinical stage does not satisfy the inclusion criteria
		capture assert `tpsa'>=0 if `touse' 
		if _rc>0 {
			disp as err "tPSA must be >= 0"
			exit
		}
		
	* calculate the linear prediction for each patient based on their predictor values
		qui g `xb_full' = 	110.52302 ///
			+ 0.39237815 * log(`tpsa') ///
			- 0.032117323 * max(log(`tpsa') - 1.0650555, 0) ^ 3 ///
			+ 0.059594888 * max(log(`tpsa') - 1.8407084, 0) ^ 3 ///
			- 0.027477565 * max(log(`tpsa') - 2.747335, 0) ^ 3  ///
			+ 1.3792909 * (`pathgg1' > 3) ///
			+ 0.57039863 * (`pathgg2' > 3) ///
			+ 0.3931648 * (`sms'==1 ) ///
			+ 0.94251512 * (`svi'==1) ///
			+ 0.43389125 * (`lni'==1) ///
			+ 0.8748934 * (`ece'==1) ///
			- 0.056105835 * min(`yos', 2005)

	* "scalar basesurv" is the 5-year baseline survival probability from the lookup table
	scalar basesurv = 0.9442806 

	* compute the 5-year survival probability for each patient based on their predictor values
	qui g `phat_full' = basesurv ^ exp(`xb_full') 
	
		if trim("`xb'")=="" {
			qui g nomogram_tempvarname798465=`phat_full' if `touse'
		}
		else if trim("`xb'")!="" {
			qui replace `g' =`xb' if `touse'
		}
		}		
		
		* list out the paper info: title, first author, journal and PMID 
		disp _newline "Nomogram Publication Information:"
		disp _column(5) "Postoperative Nomogram Predicting the 10-Year Probability of Prostate Cancer Recurrence After Radical Prostatectomy"
		disp _column(5) "Andrew J. Stephenson,Peter T. Scardino,James A. Eastham,Fernando J. Bianco Jr,Zohar A. Dotan,Christopher J. DiBlasio,Alwyn Reuther,Eric A. Klein and Michael W. Kattan"
		disp _column(5) "Journal of Clinical Oncology"
		disp _column(5) "doi: 10.1200/JCO.2005.01.867 JCO October 1, 2005 vol. 23 no. 28 7005-7012 PMID 16192588"

	end	
	
capture program drop kattan 
program define kattan
	version 12
	
	syntax [if] [in], [tpsa(varname) clinstage(varname) bxgg1(varname) bxgg2(varname) xb g(string)]
	marksample touse
	
	disp as result "KATTAN NOMOGRAM"

	qui foreach v in tpsa clinstage bxgg1 bxgg2 {
	if trim("``v''") == "" {        
			if "`v'"=="tpsa" noi disp as text "`v' assumed to be varaible `v' in ng/mL"
			else if inlist("`v'","clinstage") noi disp as text "`v' assumed to be varaible `v'"
			else if inlist("`v'","bxgg1", "bxgg2") noi disp as text "`v' assumed to be varaible `v' numeric form 5-10"
		
			*marking the observations to exclude due to missing values
			replace `touse'=. if mi(`v')
		
			local `v'="`v'"
		}
	}

	qui {
		tempvar logtpsa linpred preopnomogram
		tempname basesurv
		
		*Outputting error message if the clinical stage does not satisfy the inclusion criteria
		capture assert inlist(trim(upper(`clinstage')),"T1A","T1B","T1C","T2A","T2B","T2C","T3A","T3B","") if `touse'
		if _rc>0 {
			disp as err "Clinical stage must contain the following values only: T1A, T1B, T1C, T2A, T2B, T2C, T3A, T3B"
			exit
		}
		capture assert `tpsa'>=0 if `touse' 
		if _rc>0 {
			disp as err "tPSA must be >= 0"
			exit
		}
		
		
		g `logtpsa' = log(`tpsa')
		g `linpred' = -1.547257	///
			+ 0.31586671 * `logtpsa' 	///
			- 0.028514152 * max(`logtpsa' - 0.18255172, 0) ^ 3 ///
			+ 1.5861471 * max(`logtpsa' - 1.5260563, 0) ^ 3 	///
			- 3.426208 * max(`logtpsa' - 1.9169226, 0) ^ 3 	///
			+ 2.1286969 * max(`logtpsa' - 2.360854, 0) ^ 3 	///
			- 0.26012183 * max(`logtpsa' - 3.3565481, 0) ^ 3	///
			+ (upper(`clinstage')=="T1C") * -0.70223659 	///
			+ (upper(`clinstage')=="T2A") * -0.40668753 	///
			+ (upper(`clinstage')=="T2B") * 0.22358979 		///
			+ (upper(`clinstage')=="T2C") * 0.021192634		///
			+ (upper(`clinstage')=="T3A" | upper(`clinstage')=="T3B") * 0.7597205 ///
			+ (`bxgg1' < 3 & `bxgg2'==3) * 0.71327666	///
			+ (`bxgg1'==3 & `bxgg2' < 3) * 0.80371467	///
			+ (`bxgg1'==3 & `bxgg2'==3) * 0.77753355	///
			+ (`bxgg1' < 4 & `bxgg2' > 3) * 1.4273016 	///
			+ (`bxgg1' > 3 & `bxgg2' > 0) * 1.4957297			
		
		* "scalar basesurv" is the 5-year baseline survival probability from the lookup table
		scalar `basesurv' = 0.7800056
		
		* compute the 5-year survival probability for each patient based on their predictor values
		g `preopnomogram' = 1 - `basesurv'^exp(`linpred')	
	}

		if trim("`xb'")=="" {
			qui g nomogram_tempvarname798465=`preopnomogram' if `touse'
			*label var `g' "Risk of BCR 5-years after RP from Kattan Model"
		}
		else if trim("`xb'")!="" {
			qui g nomogram_tempvarname798465=`linpred' if `touse'
			*label var `g' "Risk of BCR 5-years after RP from Kattan Model (Linear Prediction)"
		}
		
		* list out the paper info: title, first author, journal and PMID 
		disp _newline "Nomogram Publication Information:"
		disp _column(5) "Title: A Preoperative Nomogram for Disease Recurrence Following Radical Prostatectomy for Prostate Cancer"
		disp _column(5) "Michael W. Kattan, James A. Eastham, Alan M. F. Stapleton, Thomas M. Wheeler, Peter T. Scardino"
		disp _column(5) "Journal: Journal of the National Cancer Institute"
		disp _column(5) "doi: 10.1093/jnci/90.10.766 JNCI 1998 May 20;90(10):766-71. PMID 9605647"

	end			

/*
*set trace on 
clear
set obs 100
drawnorm age tpsa fpsa ipsa hk2 pretxpsa
replace pretxpsa = abs(pretxpsa)
g double rp_gg1 = round((5-1)*runiform() + 1)
g double rp_gg2 = round((5-1)*runiform() + 1)
g double svi = round((0-1)*runiform() + 1)
g double sms = round((0-1)*runiform() + 1)
g double lni = round((0-1)*runiform() + 1)
g double ece = round((0-1)*runiform() + 1)
g double yos = round((10-1)*runiform() + 1)
replace yos = 1990+yos
replace age = . in 1
replace tpsa = . in 3
g clinstage="T2C"
egen nomopostop=nomogram(MSKCC POSTOP BCR), pathgg1(rp_gg1) pathgg2(rp_gg2)
egen nomopreop=nomogram(MSKCC PREOP BCR), bxgg1(rp_gg1) bxgg2(rp_gg2)
g poscores=round(10*runiform())
g negcores=round(10*runiform())
egen nomopreopcores=nomogram(MSKCC CORES PREOP BCR), bxgg1(rp_gg1) bxgg2(rp_gg2)
egen nomovickers=nomogram(BX VICKERS 2008)
egen nomostephenson=nomogram(BCR STEPHENSON 2005), pathgg1(rp_gg1) pathgg2(rp_gg2)
egen nomokattan=nomogram(KATTAN), bxgg1(rp_gg1) bxgg2(rp_gg2)
sum nomo*
*/


/*
PROGRAM: _gnomogram
PROGRAMMER: Melissa
DATE: 4/4/2016
DESCRIPTION: Added MSKCC dynamic nomograms and outputting edits to the ado
*/

capture program drop preopbcr 
program define preopbcr
	version 12
	
	syntax [if] [in], [tpsa(varname) clinstage(varname) bxgg1(varname) bxgg2(varname) xb g(string)]
	marksample touse
	
	disp as result "MSKCC PREOPERATIVE BCR NOMOGRAM"

	qui foreach v in tpsa clinstage bxgg1 bxgg2 {
	if trim("``v''") == "" {        
			if "`v'"=="tpsa" noi disp as text "`v' assumed to be variable `v' in ng/mL"
			else if inlist("`v'","clinstage") noi disp as text "`v' assumed to be variable `v'"
			else if inlist("`v'","bxgg1", "bxgg2") noi disp as text "`v' assumed to be variable `v' numeric form 5-10"
			
			*marking the observations to exclude due to missing values
			replace `touse'=. if mi(`v')
		
			local `v'="`v'"

		}
	}

	qui {
		tempvar sptpsa1 sptpsa2 linpred preopnomogram 
		tempname scaling


		*Outputting error message if the clinical stage does not satisfy the inclusion criteria
		capture assert inlist(trim(upper(`clinstage')),"T1C","T2A","T2B","T2C","T3A","T3B","") if `touse'
		if _rc>0 {
			disp as err "Clinical stage must contain the following values only: T1C, T2A, T2B, T2C, T3A, T3B"
			exit
		}
		capture assert `tpsa'>=0 if `touse' 
		if _rc>0 {
			disp as err "tPSA must be >= 0"
			exit
		}
		
		
		*Generating spline terms for PSA
		local k1tpsa=0.21
		local k2tpsa=4.6
		local k3tpsa=7
		local k4tpsa=96.53
		
		gen `sptpsa1'= (max(`tpsa'-scalar(`k1tpsa'),0))^3 - ((max(`tpsa'-scalar(`k3tpsa'),0))^3*((scalar(`k4tpsa')-scalar(`k1tpsa'))/(scalar(`k4tpsa')-scalar(`k3tpsa')))) + ((max(`tpsa'-scalar(`k4tpsa'),0))^3*((scalar(`k3tpsa')-scalar(`k1tpsa'))/(scalar(`k4tpsa')-scalar(`k3tpsa'))))
		gen `sptpsa2'= (max(`tpsa'-scalar(`k2tpsa'),0))^3 - ((max(`tpsa'-scalar(`k3tpsa'),0))^3*((scalar(`k4tpsa')-scalar(`k2tpsa'))/(scalar(`k4tpsa')-scalar(`k3tpsa')))) + ((max(`tpsa'-scalar(`k4tpsa'),0))^3*((scalar(`k3tpsa')-scalar(`k2tpsa'))/(scalar(`k4tpsa')-scalar(`k3tpsa'))))

		*generating the linear predictor
		g `linpred' = 6.87501087	///
			+ -0.47546268 * `tpsa' 	///
			+ 0.00397062 * `sptpsa1' ///
			+ -0.01100004 * `sptpsa2'	///
			+ -2.17451547 * (`bxgg1'>3)  ///
			+ -0.88163218 * (`bxgg2'>3) ///
			+ (upper(`clinstage')=="T2A") * -0.53530558 ///
			+ (upper(`clinstage')=="T2B") * -1.02585827 ///
			+ (upper(`clinstage')=="T2C") * -0.73522864	///
			+ (upper(`clinstage')=="T3A" | upper(`clinstage')=="T3B") * -1.15236139 
			
			
		* "scaling parater for log- logistic model
		scalar `scaling' = 1.15475454
		
		* compute the 5-year survival probability for each patient based on their predictor values
		g `preopnomogram' = (1 + (5*exp(-1*`linpred'))^(1/`scaling'))^-1
	}

		if trim("`xb'")=="" {
			qui g nomogram_tempvarname798465=`preopnomogram' if `touse'
			*label var `g' "MSKCC PreOp 5 year Risk of BCR"
		}
		else if trim("`xb'")!="" {
			qui g nomogram_tempvarname798465=`linpred' if `touse'
			*label var `g' "MSKCC Pre-Op 5 year Risk of BCR (Linear Prediction)"
		}
		
		* list out the paper info: title, first author, journal and PMID 
		disp _newline "Nomogram Publication Information:"
		disp _column(5) "MSKCC Dynamic Preoperative BCR Nomogram"
		disp _column(5) "Website: https://www.mskcc.org/nomograms/prostate/pre-op/coefficients"
		disp _column(5) "Date of Download: 3/15/2016"

	end					
			

			
capture program drop preopbcrcores 
program define preopbcrcores
	version 12
	
	syntax [if] [in], [tpsa(varname) clinstage(varname) bxgg1(varname) bxgg2(varname) poscores(varname) negcores(varname) xb g(string)]
	marksample touse
	
	disp as result "MSKCC PREOPERATIVE BCR NOMOGRAM WITH CORES"

	qui foreach v in tpsa clinstage bxgg1 bxgg2 poscores negcores {
	if trim("``v''") == "" {        
			if "`v'"=="tpsa" noi disp as text "`v' assumed to be variable `v' in ng/mL"
			else if inlist("`v'","clinstage") noi disp as text "`v' assumed to be variable `v'"
			else if inlist("`v'","bxgg1", "bxgg2") noi disp as text "`v' assumed to be variable `v' numeric form 5-10"
			else if inlist("`v'","poscores", "negcores") noi disp as text "`v' assumed to be variable `v' numeric form"
		
			*marking the observations to exclude due to missing values
			replace `touse'=. if mi(`v')
		
			local `v'="`v'"
		}
	}

	qui {
		tempvar sptpsa1 sptpsa2 linpred preopnomogram   
		tempname scaling


		*Outputting error message if the clinical stage does not satisfy the inclusion criteria
		capture assert inlist(trim(upper(`clinstage')),"T1C","T2A","T2B","T2C","T3A","T3B","") if `touse'
		if _rc>0 {
			disp as err "Clinical stage must contain the following values only: T1C, T2A, T2B, T2C, T3A, T3B"
			exit
		}
		capture assert `tpsa'>=0 if `touse' 
		if _rc>0 {
			disp as err "tPSA must be >= 0"
			exit
		}
		
		
		*Generating spline terms for PSA
		local k1tpsa=0.21
		local k2tpsa=4.6
		local k3tpsa=7
		local k4tpsa=96.53
		
		gen `sptpsa1'= (max(`tpsa'-scalar(`k1tpsa'),0))^3 - ((max(`tpsa'-scalar(`k3tpsa'),0))^3*((scalar(`k4tpsa')-scalar(`k1tpsa'))/(scalar(`k4tpsa')-scalar(`k3tpsa')))) + ((max(`tpsa'-scalar(`k4tpsa'),0))^3*((scalar(`k3tpsa')-scalar(`k1tpsa'))/(scalar(`k4tpsa')-scalar(`k3tpsa'))))
		gen `sptpsa2'= (max(`tpsa'-scalar(`k2tpsa'),0))^3 - ((max(`tpsa'-scalar(`k3tpsa'),0))^3*((scalar(`k4tpsa')-scalar(`k2tpsa'))/(scalar(`k4tpsa')-scalar(`k3tpsa')))) + ((max(`tpsa'-scalar(`k4tpsa'),0))^3*((scalar(`k3tpsa')-scalar(`k2tpsa'))/(scalar(`k4tpsa')-scalar(`k3tpsa'))))

		*generating the linear predictor
		g `linpred' = 7.01406658	///
			+ -0.69296928 * `tpsa' 	///
			+ 0.00655692 * `sptpsa1' ///
			+ -0.01835119 * `sptpsa2'	///
			+ -2.13085826 * (`bxgg1'>3)  ///
			+ -0.83564678 * (`bxgg2'>3) ///
			+ (upper(`clinstage')=="T2A") * -0.31238428 ///
			+ (upper(`clinstage')=="T2B") * -0.69059638 ///
			+ (upper(`clinstage')=="T2C") * -0.38502 ///
			+ (upper(`clinstage')=="T3A" | upper(`clinstage')=="T3B") * -0.66371692 ///
			+ 0.06702081 * `negcores' ///
			+ -0.01624338 * `poscores'
			
		
		* "scaling parater for log- logistic model
		scalar `scaling' = 1.10758467
		
		* compute the 5-year survival probability for each patient based on their predictor values
		g `preopnomogram' = (1 + (5*exp(-1*`linpred'))^(1/`scaling'))^-1
	}

		if trim("`xb'")=="" {
			qui g nomogram_tempvarname798465=`preopnomogram' if `touse'
			*label var `g' "MSKCC Pre-Op (Cores) 5-year Risk of BCR"
		}
		else if trim("`xb'")!="" {
			qui g nomogram_tempvarname798465=`linpred' if `touse'
			*label var `g' "MSKCC Pre-Op (Cores) 5-year Risk of BCR (Linear Prediction)"
		}
		
		* list out the paper info: title, first author, journal and PMID 
		disp _newline "Nomogram Publication Information:"
		disp _column(5) "MSKCC Dynamic Preoperative BCR Nomogram with Cores"
		disp _column(5) "Website: https://www.mskcc.org/nomograms/prostate/pre-op/coefficients"
		disp _column(5) "Date of Download: 3/15/2016"
		

	end					
		
capture program drop postopbcr 
program define postopbcr 
	version 12
	
	syntax [if] [in], [tpsa(varname) pathgg1(varname) pathgg2(varname) ece(varname) svi(varname) lni(varname) sms(varname) xb g(string)]
	marksample touse
	
	disp as result "MSKCC POSTOPERATIVE BCR NOMOGRAM"

	qui foreach v in tpsa pathgg1 pathgg2 ece svi lni sms {
	if trim("``v''") == "" {        
			if inlist("`v'","tpsa") noi disp as text "`v' assumed to be variable `v' in ng/mL"
			else if inlist("`v'","pathgg1", "pathgg2") noi disp as text "`v' assumed to be variable `v' in numeric format 0-5"
			else if inlist("`v'","ece", "svi", "lni", "sms") noi disp as text "`v' assumed to be variable `v' in binary format"
		
			*marking the observations to exclude due to missing values
			replace `touse'=. if mi(`v')
		
			local `v'="`v'"
		}
	}
	
	qui {
		tempvar sptpsa1 sptpsa2 linpred preopnomogram 
			tempname scaling
		
		*Outputting error message if the clinical stage does not satisfy the inclusion criteria
		capture assert `tpsa'>=0 if `touse' 
		if _rc>0 {
			disp as err "tPSA must be >= 0"
			exit
		}
		
		
		*Generating spline terms for PSA
		local k1tpsa=0.21
		local k2tpsa=4.6
		local k3tpsa=7
		local k4tpsa=96.53

		gen `sptpsa1'= (max(`tpsa'-scalar(`k1tpsa'),0))^3 - ((max(`tpsa'-scalar(`k3tpsa'),0))^3*((scalar(`k4tpsa')-scalar(`k1tpsa'))/(scalar(`k4tpsa')-scalar(`k3tpsa')))) + ((max(`tpsa'-scalar(`k4tpsa'),0))^3*((scalar(`k3tpsa')-scalar(`k1tpsa'))/(scalar(`k4tpsa')-scalar(`k3tpsa'))))
		gen `sptpsa2'= (max(`tpsa'-scalar(`k2tpsa'),0))^3 - ((max(`tpsa'-scalar(`k3tpsa'),0))^3*((scalar(`k4tpsa')-scalar(`k2tpsa'))/(scalar(`k4tpsa')-scalar(`k3tpsa')))) + ((max(`tpsa'-scalar(`k4tpsa'),0))^3*((scalar(`k3tpsa')-scalar(`k2tpsa'))/(scalar(`k4tpsa')-scalar(`k3tpsa'))))
	

/* 		* 5-YEAR NOMOGRAM PROBABILITY
		Variables required for the calculation are:
			Pretreatment PSA													* tpsa  
			Primary and secondary Gleason scores from pathology specimen        * pathgg1 & pathgg2
			Extracapsular extension (yes or no)                                 * ece
			Seminal vesicle invasion (yes or no)                                * svi
			Lymph node involvement (yes or no)                                  * lni
			Positive surgical margins (yes or no)                               * sms
		*/    
		
	*generating the linear predictor
		g `linpred' = 6.39812614	///
			+ -0.29687539 * `tpsa' 	///
			+ 0.00271688 * `sptpsa1' ///
			+ -0.00757855 * `sptpsa2'	///
			+ -2.28938706 * (`pathgg1'>3)  ///
			+ -0.86724004 * (`pathgg2'>3) ///
			+ -0.87752933 * (`sms'==1 ) ///
			+ -0.62873088 * (`svi'==1) ///
			+ -1.20264841 * (`lni'==1) ///
			+ -0.50294333 * (`ece'==1) 
			
			
		* "scaling parater for log- logistic model
		scalar `scaling' = 1.0151869
		
		* compute the 5-year survival probability for each patient based on their predictor values
		g `preopnomogram' = (1 + (5*exp(-1*`linpred'))^(1/`scaling'))^-1
	}

		if trim("`xb'")=="" {
			qui g nomogram_tempvarname798465=`preopnomogram' if `touse'
			*label var `g' "MSKCC Post-Op 5-year Risk of BCR"
		}
		else if trim("`xb'")!="" {
			qui g nomogram_tempvarname798465=`linpred' if `touse'
			*label var `g' "MSKCC Post-Op 5-year Risk of BCR (Linear Prediction)"
		}
			
		
		* list out the paper info: title, first author, journal and PMID 
		disp _newline "Nomogram Publication Information:"
		disp _column(5) "MSKCC Dynamic Postoperative BCR Nomogram"
		disp _column(5) "Website: https://www.mskcc.org/nomograms/prostate/pre-op/coefficients"
		disp _column(5) "Date of Download: 3/15/2016"
	

	end	

