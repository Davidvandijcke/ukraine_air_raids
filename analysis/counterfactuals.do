**~~~~~~~~~~~~~~~~**~~~~~~~~~~~~~~~~**~~~~~~~~~~~~~~~~**~~~~~~~~~~~~~~~~**
* 0. Set the working directory below *
**~~~~~~~~~~~~~~~~**~~~~~~~~~~~~~~~~**~~~~~~~~~~~~~~~~**~~~~~~~~~~~~~~~~**
	** cd "~/repo/"
	cd "/Users/austinlw/Harris-Public-Policy Dropbox/Austin Wright/Ukraine_AirRaids/repo/"

**~~~~~~~~~~~~~~~~**~~~~~~~~~~~~~~~~**~~~~~~~~~~~~~~~~**~~~~~~~~~~~~~~~~**
**~~~~~~~~~~~~~~~~**~~~~~~~~~~~~~~~~**~~~~~~~~~~~~~~~~**~~~~~~~~~~~~~~~~**
** 1. Storing scalars for counterfactuals  
**~~~~~~~~~~~~~~~~**~~~~~~~~~~~~~~~~**~~~~~~~~~~~~~~~~**~~~~~~~~~~~~~~~~**
**~~~~~~~~~~~~~~~~**~~~~~~~~~~~~~~~~**~~~~~~~~~~~~~~~~**~~~~~~~~~~~~~~~~**

** Scalar for correlation between civilian casualty events in VIINA and post-alert movement 
	scalar mvmt_coef = -0.0000674
	
** Scalar for movement summary statistics at beginning + end of sample period
	scalar mvmt_start = 1938
	scalar mvmt_end = 421.85
	
**~~~~~~~~~~~~~~~~**~~~~~~~~~~~~~~~~**~~~~~~~~~~~~~~~~**~~~~~~~~~~~~~~~~**
**~~~~~~~~~~~~~~~~**~~~~~~~~~~~~~~~~**~~~~~~~~~~~~~~~~**~~~~~~~~~~~~~~~~**
** 2. Processing ACLED data on civilian casualties
**~~~~~~~~~~~~~~~~**~~~~~~~~~~~~~~~~**~~~~~~~~~~~~~~~~**~~~~~~~~~~~~~~~~**
**~~~~~~~~~~~~~~~~**~~~~~~~~~~~~~~~~**~~~~~~~~~~~~~~~~**~~~~~~~~~~~~~~~~**
	
	** Load ACLED data 
	import excel "data/raw/acled/Ukraine_Black_Sea_2020_2022.xlsx", sheet("Sheet1") firstrow clear

	** Build civilian harm measures
	gen civcas=FATALITIES if ACTOR2=="Civilians (Ukraine)"
	gen civcas_event=1 if ACTOR2=="Civilians (Ukraine)"&FATALITIES>=1&FATALITIES!=.

	collapse (sum) civcas civcas_event, by(EVENT_DATE)

	** Set time series structure of aggregated ACLED data 
	tsset EVENT_DATE
		gen year=year(EVENT_DATE)
		gen month=month(EVENT_DATE)
		gen day=day(EVENT_DATE)
	
	** Generate civilian casualty elasticities (day-specific + moving average)
	gen civcas_elasticity=civcas/civcas_event
	gen cc_elas_7dayMovAvg= (civcas_elasticity + L.civcas_elasticity + L2.civcas_elasticity + L3.civcas_elasticity + L4.civcas_elasticity + L5.civcas_elasticity + L6.civcas_elasticity)/7

	** Reproduce time series from movement data 
	keep if year>=2022&month>=3
	drop if month==3&day<=14
	keep if month<10
	
**~~~~~~~~~~~~~~~~**~~~~~~~~~~~~~~~~**~~~~~~~~~~~~~~~~**~~~~~~~~~~~~~~~~**
**~~~~~~~~~~~~~~~~**~~~~~~~~~~~~~~~~**~~~~~~~~~~~~~~~~**~~~~~~~~~~~~~~~~**
** 3. Combine time series data on movement and conflict activity + harm
**~~~~~~~~~~~~~~~~**~~~~~~~~~~~~~~~~**~~~~~~~~~~~~~~~~**~~~~~~~~~~~~~~~~**
**~~~~~~~~~~~~~~~~**~~~~~~~~~~~~~~~~**~~~~~~~~~~~~~~~~**~~~~~~~~~~~~~~~~**

	** Generate matching date 
	gen smart_date = mdy(month, day, year)
		sort smart_date 
	
	merge smart_date using "data/final/alert_countXday.dta"

	** Construct period of conflict indicator
	gen period_of_conflict=1 if month==3|month==4
		replace period_of_conflict=2 if month==5|month==6
		replace period_of_conflict=3 if month>=7
	
	** Construct interpolated movement 
	gen mvmt_avg = mvmt_start if month==3&day==15
		replace mvmt_avg = mvmt_end if month==9&day==30
		ipolate mvmt_avg EVENT_DATE, generate(mvmt_trend)

**~~~~~~~~~~~~~~~~**~~~~~~~~~~~~~~~~**~~~~~~~~~~~~~~~~**~~~~~~~~~~~~~~~~**
**~~~~~~~~~~~~~~~~**~~~~~~~~~~~~~~~~**~~~~~~~~~~~~~~~~**~~~~~~~~~~~~~~~~**
** 4. Counterfactuals: avoided harm (full sample) and excess harm (periods 2 and 3)
**~~~~~~~~~~~~~~~~**~~~~~~~~~~~~~~~~**~~~~~~~~~~~~~~~~**~~~~~~~~~~~~~~~~**
**~~~~~~~~~~~~~~~~**~~~~~~~~~~~~~~~~**~~~~~~~~~~~~~~~~**~~~~~~~~~~~~~~~~**
	
	****************************
	** Avoided harm exercises ** 
	****************************
	
	** Avoided harm, scenario 1: no movement 
	gen deaths_avoided_1 		= alert	* civcas_elasticity * (0 - movement) * mvmt_coef 
	** Avoided harm, scenario 2: minimal movement
	gen deaths_avoided_2 	= alert	* civcas_elasticity * (421.8466 - movement) * mvmt_coef 
		
	***************************
	** Excess harm exercises ** 
	***************************

	** Trend relative movement during first period (high compliance)
	gen mvmt_trend_delta = mvmt_trend - mvmt_start
	** Day-specific movement relative movement during first period (high compliance)
	gen mvmt_actual_delta = movement - mvmt_start

	** Excess harm, scenario 1: movement trend, day-specific civilian harm elasticity 
	gen excess_deaths_1 		= alert * civcas_elasticity * mvmt_trend_delta * mvmt_coef 
	** Excess harm, scenario 2: day-specific movement, day-specific civilian harm elasticity 
	gen excess_deaths_2 	= alert * civcas_elasticity * mvmt_actual_delta * mvmt_coef 
	** Excess harm, scenario 3: movement trend, moving average of civilian harm elasticity 
	gen excess_deaths_3 	= alert * cc_elas_7dayMovAvg * mvmt_trend_delta * mvmt_coef 
	
		foreach v in excess_deaths_1 excess_deaths_2 excess_deaths_3{
			replace `v' =. if period_of_conflict==1
		}
		
	keep deaths_avoided_1 deaths_avoided_2 excess_deaths_1 excess_deaths_2 excess_deaths_3 civcas smart_date month period_of_conflict
	
	save "data/final/counterfactuals.dta", replace 

\\

**~~~~~~~~~~~~~~~~**~~~~~~~~~~~~~~~~**~~~~~~~~~~~~~~~~**~~~~~~~~~~~~~~~~**
**~~~~~~~~~~~~~~~~**~~~~~~~~~~~~~~~~**~~~~~~~~~~~~~~~~**~~~~~~~~~~~~~~~~**
** 5. Avoided harm exercises 
**~~~~~~~~~~~~~~~~**~~~~~~~~~~~~~~~~**~~~~~~~~~~~~~~~~**~~~~~~~~~~~~~~~~**
**~~~~~~~~~~~~~~~~**~~~~~~~~~~~~~~~~**~~~~~~~~~~~~~~~~**~~~~~~~~~~~~~~~~**
	
	use "data/final/counterfactuals.dta", clear 

	preserve 
	collapse (sum) deaths_avoided_1 deaths_avoided_2 civcas
	// deaths_avoided_1		deaths_avoided_2		civcas
	// 1617.007				1269.027				3617
	restore 	

	** set sample
	keep if deaths_avoided_1!=.

**~~~~~~~~~~~~~~~~**~~~~~~~~~~~~~~~~**~~~~~~~~~~~~~~~~**~~~~~~~~~~~~~~~~**~~~~~~~~~~~~~~~~**
** BOOTSTRAP WITH REPLACEMENT
**~~~~~~~~~~~~~~~~**~~~~~~~~~~~~~~~~**~~~~~~~~~~~~~~~~**~~~~~~~~~~~~~~~~**~~~~~~~~~~~~~~~~**

preserve 

set matsize 10000

matrix store_shuffle = J(10000, 4, .) 

forval i=1/10000{
	
	preserve 
	
	set seed `i'
	
	matrix store_shuffle[`i', 1] = `i'

	quietly bsample
		
	collapse (sum) deaths_avoided_1 deaths_avoided_2 civcas
	
	matrix store_shuffle[`i', 2] = deaths_avoided_1[1]
	matrix store_shuffle[`i', 3] = deaths_avoided_2[1]
	matrix store_shuffle[`i', 4] = civcas[1]
	
	restore 
}

svmat store_shuffle, names(shuf_deaths)

 	gen deaths_avoided_1_rate=shuf_deaths2/shuf_deaths4
 	gen deaths_avoided_2_rate=shuf_deaths3/shuf_deaths4

	summarize deaths_avoided_1_rate	
	summarize deaths_avoided_2_rate
		
restore 
		
**~~~~~~~~~~~~~~~~**~~~~~~~~~~~~~~~~**~~~~~~~~~~~~~~~~**~~~~~~~~~~~~~~~~**~~~~~~~~~~~~~~~~**		
** BOOTSTRAP WITH REPLACEMENT + PERIOD CLUSTERS
**~~~~~~~~~~~~~~~~**~~~~~~~~~~~~~~~~**~~~~~~~~~~~~~~~~**~~~~~~~~~~~~~~~~**~~~~~~~~~~~~~~~~**

preserve 

set matsize 10000

matrix store_shuffle = J(10000, 4, .) 

forval i=1/10000{
	
	preserve 
	
	set seed `i'
	
	matrix store_shuffle[`i', 1] = `i'

	quietly bsample, cluster(period_of_conflict)
		
	collapse (sum) deaths_avoided_1 deaths_avoided_2 civcas
	
	matrix store_shuffle[`i', 2] = deaths_avoided_1[1]
	matrix store_shuffle[`i', 3] = deaths_avoided_2[1]
	matrix store_shuffle[`i', 4] = civcas[1]
	
	restore 
}

svmat store_shuffle, names(shuf_deaths)

	gen deaths_avoided_1_rate=shuf_deaths2/shuf_deaths4
 	gen deaths_avoided_2_rate=shuf_deaths3/shuf_deaths4

	summarize deaths_avoided_1_rate	
	summarize deaths_avoided_2_rate
	
restore 
	
**~~~~~~~~~~~~~~~~**~~~~~~~~~~~~~~~~**~~~~~~~~~~~~~~~~**~~~~~~~~~~~~~~~~**~~~~~~~~~~~~~~~~**		
** BOOTSTRAP WITH REPLACEMENT + MONTHLY CLUSTERS
**~~~~~~~~~~~~~~~~**~~~~~~~~~~~~~~~~**~~~~~~~~~~~~~~~~**~~~~~~~~~~~~~~~~**~~~~~~~~~~~~~~~~**

preserve 

set matsize 10000

matrix store_shuffle = J(10000, 4, .) 

forval i=1/10000{
	
	preserve 
	
	set seed `i'
	
	matrix store_shuffle[`i', 1] = `i'

	quietly bsample, cluster(month)
		
	collapse (sum) deaths_avoided_1 deaths_avoided_2 civcas
	
	matrix store_shuffle[`i', 2] = deaths_avoided_1[1]
	matrix store_shuffle[`i', 3] = deaths_avoided_2[1]
	matrix store_shuffle[`i', 4] = civcas[1]
	
	restore 
}

svmat store_shuffle, names(shuf_deaths)

	gen deaths_avoided_1_rate=shuf_deaths2/shuf_deaths4
 	gen deaths_avoided_2_rate=shuf_deaths3/shuf_deaths4

	summarize deaths_avoided_1_rate	
	summarize deaths_avoided_2_rate

restore 
	
**~~~~~~~~~~~~~~~~**~~~~~~~~~~~~~~~~**~~~~~~~~~~~~~~~~**~~~~~~~~~~~~~~~~**~~~~~~~~~~~~~~~~**	
** PERMUTED SAMPLING WITHOUT REPLACEMENT
**~~~~~~~~~~~~~~~~**~~~~~~~~~~~~~~~~**~~~~~~~~~~~~~~~~**~~~~~~~~~~~~~~~~**~~~~~~~~~~~~~~~~**
 
preserve 

set matsize 10000

// set seed 1234

matrix store_shuffle = J(10000, 4, .) 

// drop if deaths_avoided_ce==.

forval i=1/10000{
	
	preserve 
	
	set seed `i'
	
	matrix store_shuffle[`i', 1] = `i'

	quietly gen random = runiform()
		sort random 
		quietly keep if _n<=90
		
	collapse (sum) deaths_avoided_1 deaths_avoided_2 civcas
	
	matrix store_shuffle[`i', 2] = deaths_avoided_1[1]
	matrix store_shuffle[`i', 3] = deaths_avoided_2[1]
	matrix store_shuffle[`i', 4] = civcas[1]
	
	restore 
}

svmat store_shuffle, names(shuf_deaths)

	gen deaths_avoided_1_rate=shuf_deaths2/shuf_deaths4
 	gen deaths_avoided_2_rate=shuf_deaths3/shuf_deaths4

	summarize deaths_avoided_1_rate	
	summarize deaths_avoided_2_rate
		
restore 
		
**~~~~~~~~~~~~~~~~**~~~~~~~~~~~~~~~~**~~~~~~~~~~~~~~~~**~~~~~~~~~~~~~~~~**
**~~~~~~~~~~~~~~~~**~~~~~~~~~~~~~~~~**~~~~~~~~~~~~~~~~**~~~~~~~~~~~~~~~~**
** 6. Excess harm exercises 
**~~~~~~~~~~~~~~~~**~~~~~~~~~~~~~~~~**~~~~~~~~~~~~~~~~**~~~~~~~~~~~~~~~~**
**~~~~~~~~~~~~~~~~**~~~~~~~~~~~~~~~~**~~~~~~~~~~~~~~~~**~~~~~~~~~~~~~~~~**
	
	use "data/final/counterfactuals.dta", clear

	preserve 
		collapse (sum) excess_deaths_1 excess_deaths_2 excess_deaths_3 civcas if period_of_conflict>=2
		// excess_deaths_1	excess_deaths_3		excess_deaths_2	civcas
		// 276.5217			216.4506			149.7302			1817
	restore 

	** set sample
	keep if period_of_conflict>=2

	format smart_date %tdMon_dd
	tsset smart_date

tw (tsline excess_deaths_1, lcolor(green)) (tsline excess_deaths_2, lcolor(red)) (tsline excess_deaths_3, lcolor(blue)), ///
title("") xtitle("") ytitle("Excess deaths (total)") ///
	legend(size(vsmall) pos(12) ring(0) col(1) order(1 "Mvmt Trend + Casualty Raw       = 15% Excess Deaths" 2 "Mvmt Raw + Casualty Raw         =  8%  Excess Deaths" 3 "Mvmt Trend + Casualty Trend     = 12% Excess Deaths" )) 
	
	graph export "output/figs/ts_excessdeaths.png", replace
	
**~~~~~~~~~~~~~~~~**~~~~~~~~~~~~~~~~**~~~~~~~~~~~~~~~~**~~~~~~~~~~~~~~~~**~~~~~~~~~~~~~~~~**		
** BOOTSTRAP WITH REPLACEMENT + MONTHLY CLUSTERS
**~~~~~~~~~~~~~~~~**~~~~~~~~~~~~~~~~**~~~~~~~~~~~~~~~~**~~~~~~~~~~~~~~~~**~~~~~~~~~~~~~~~~**

set matsize 10000

matrix store_shuffle = J(10000, 5, .) 

forval i=1/10000{
	
	preserve 
	
	set seed `i'
	
	matrix store_shuffle[`i', 1] = `i'

	quietly bsample, cluster(month)
		
	collapse (sum)  excess_deaths_1 excess_deaths_2 excess_deaths_3 civcas
	
	matrix store_shuffle[`i', 2] = excess_deaths_1[1]
	matrix store_shuffle[`i', 3] = excess_deaths_2[1]
	matrix store_shuffle[`i', 4] = excess_deaths_3[1]
	matrix store_shuffle[`i', 5] = civcas[1]
	
	restore 
}

svmat store_shuffle, names(shuf_deaths)

 	gen excess_deaths_1_rate=shuf_deaths2/shuf_deaths5
 	gen excess_deaths_2_rate=shuf_deaths3/shuf_deaths5
 	gen excess_deaths_3_rate=shuf_deaths4/shuf_deaths5

	summarize excess_deaths_1_rate
	summarize excess_deaths_2_rate	
	summarize excess_deaths_3_rate	

	
