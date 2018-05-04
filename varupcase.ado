/*
PROGRAM: varupcase
PROGRAMMER: Daniel
DATE: 5/4/2011
This program turns all variable names into uppercase, and all character strings into uppercase as well.
*/

program varupcase
	version 11.0
	syntax varlist

	foreach v of varlist `varlist' {
		quietly capture replace `v' = upper(`v')
		quietly capture rename `v' `=upper("`v'")'
	}

end
