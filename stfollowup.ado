* DATA NEEDS TO BE STSET BEFORE RUNNING THIS PROGRAM

* 9/16/2015: Updated to provide IQR in additional to median.

* this ado file will give the number of patients and events and the median followup time for survivors
capture program drop stfollowup
program stfollowup, byable(recall)
	quietly{
		* "touse" is the subset indicator when the estimates are desired by a particular variable
		marksample touse
				
		count if _t0==0 & `touse'
		local total=r(N)
		count if _d==1 & `touse'
		local event=r(N)
		
		sum _t if _d==0 & `touse', d
		local followup=r(p50)
		local p25=r(p25)
		local p75=r(p75)
	}
	disp "No. patients " `total'
	disp "No. events " `event'
	disp "Median followup time for survivors is " round(`followup', 0.01) " (IQR " round(`p25', 0.01) ", " round(`p75', 0.01) ")"
end


