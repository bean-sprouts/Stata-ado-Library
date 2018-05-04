/*
Name: printlabel
Author: Jessica Poon
Date: 4/23/2013
Update by Melissa on 11/16/2017 to incorporate standard JRSC variables
*/

capture program drop printlabel
program printlabel
	version 12
	
	*syntax requires a location to save the created file and the disease we are interested in, only replace existing file if specified
	syntax , SAVing(string asis) [FORDISease(string) replace]
	preserve
	
	*give an error code if disease not specified as one of the disease types available
	*though technically an optional specification, this allows us to give a more specific error message.
	if trim("`fordisease'")=="" {
		disp as error "Disease must be specified as prostate, bladder, kidney or JRSC"
		exit
	}
	
	qui {
		*import the excel of all the common variables we use along with their associated labels and formats
		import excel using "O:\Outcomes\Andrew\_UROLOGY Biostatistics\SOPs, guidelines etc\New Project Template\Derived Variables - Template.xlsx", firstrow clear sheet("Derived Variables")
		*only keep those that match our disease of interest or basic stuff that applies to all diseases
		keep if index(upper(disease),upper("`fordisease'")) | index(upper(disease),upper("all"))
		*now we no longer need the disease column
		drop disease nts format
		gen codeadded=""
		*this is the excel equation that shows the user what format they end up with
		*because of how ugly it is, we are only showing the format in the first line: the user will have to pull the format down through the column to use it
		gen format=`"=IF(ISERROR(VLOOKUP(D"'+trim(string(_n+1))+`",Formats!A:B,2,FALSE))," ",VLOOKUP(D"'+trim(string(_n+1))+`",Formats!A:B,2,FALSE))"' in 1
		order codeadded varname label formatname format

		*export according to the file specification in a sheet called Derived Variables, including first row variables 
		capture export excel using `saving', `replace' sheet("Derived Variables") firstrow(variables)
		*decided to give our own error message for this because the Stata error message tells the user that the SHEET is already in use,
		*no matter what file we want to save it as and also suggests that we use sheetmodify or sheetreplace. This is not what we want the user to do.
		if _rc==602 {
			disp as error "This filename is already in use in current directory, must specify replace option if you want to overwrite this file."
			exit
		}
		
		//cut down on the list of formats for the second page
		*only keep the formatnames that we actually used for this format
		keep formatname
		bysort formatname: keep if _n==1 & trim(formatname)!=""
		tempfile formats
		save `formats'
			
		*import the second page of the template excel
		import excel using "O:\Outcomes\Andrew\_UROLOGY Biostatistics\SOPs, guidelines etc\New Project Template\Derived Variables - Template.xlsx", firstrow clear sheet("Formats")
		*merge to see which formats we actually use
		merge 1:1 formatname using `formats'
		*there shouldn't be any format names on the first page that wasn't on the second page
		assert _merge!=2
		*only keep those that showed up on the first page and the second page
		keep if _merge==3
		drop _merge
	}
	*now export this as a second page to the excel we had made with all the standard variable names
	export excel using `saving', sheet("Formats") sheetmodify firstrow(variables)
	
	restore
end
