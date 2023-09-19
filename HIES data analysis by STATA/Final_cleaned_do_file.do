clear all
set more off
cd "C:\Users\USER\OneDrive\Documents\487-paper" // Working folder 


***------------------------------------------------------------------------***

// ----------------     Merging all the HH information     -----------------// 

***------------------------------------------------------------------------***


* generating unique id by hhold - hhid value && keeping only relevant variables *

// HH_sec_1A.dta  --  Household Information Roster
use "C:\Users\USER\Downloads\HIES2016 24_7_2019-20230807T192137Z-001\HIES2016 24_7_2019\HH_SEC_1A.dta" , clear
drop _merge s1aq08 s1aq09 s1aq10 s1aq11 s1aq12 s1aq13 s1aq14 s1aq15 s1aq16 s1aq17
gen unid = string(hhold) + string(hhid)
save HH_SEC_1A__updated, replace

//HH_SEC_2A  --  Section 2 (Education), Part A (Literacy and Education Attainment)
use "C:\Users\USER\Downloads\HIES2016 24_7_2019-20230807T192137Z-001\HIES2016 24_7_2019\HH_SEC_2A.dta" , clear
drop _merge s2aq01 s2aq02 s2aq05 s2aq06 s2aq07
gen unid = string(hhold) + string(hhid)
save HH_SEC_2A__updated, replace

//HH_SEC_2B  --  Section 2 (Education), Part B (Current Enrollment)
use "C:\Users\USER\Downloads\HIES2016 24_7_2019-20230807T192137Z-001\HIES2016 24_7_2019\HH_SEC_2B.dta" , clear
drop s2bq02 s2bq0a s2bq04 s2bq05 s2bq06 s2bq07 s2bq08a s2bq08b s2bq08c s2bq08d s2bq08e s2bq08f s2bq08g s2bq08h s2bq08i s2bq08j s2bq08k s2bq08l s2bq08m s2bq08n s2bq08o s2bq08p s2bq08q _merge
gen unid = string(hhold) + string(hhid)
save HH_SEC_2B__updated, replace

//HH_SEC_6B  --  Section 2 (Education), Part B (Current Enrollment)
use "C:\Users\USER\Downloads\HIES2016 24_7_2019-20230807T192137Z-001\HIES2016 24_7_2019\HH_SEC_6B.dta" , clear
drop s2bq02 s2bq0a s2bq04 s2bq05 s2bq06 s2bq07 s2bq08a s2bq08b s2bq08c s2bq08d s2bq08e s2bq08f s2bq08g s2bq08h s2bq08i s2bq08j s2bq08k s2bq08l s2bq08m s2bq08n s2bq08o s2bq08p s2bq08q _merge
gen unid = string(hhold) + string(hhid)
save HH_SEC_6B__updated, replace

//HH_SEC_7A  --  Section 7 (Agriculture), Part A (Landholding)
use "C:\Users\USER\Downloads\HIES2016 24_7_2019-20230807T192137Z-001\HIES2016 24_7_2019\HH_SEC_7A.dta" , clear
drop _merge s7aq06 s7aq05 s7aq04
gen unid = string(hhold) + string(hhid)
save HH_SEC_7A__updated, replace

//HH_SEC_8C  --  Section 8 (Other Assets amd Income), Part C (Migration and Remittance)
use "C:\Users\USER\Downloads\HIES2016 24_7_2019-20230807T192137Z-001\HIES2016 24_7_2019\HH_SEC_8C.dta", clear
drop s8cq00 s8cq03 s8cq04 s8cq05a s8cq05b s8cq06 s8cq07 s8cq08 s8cq09 s8cq12 s8cq13 s8cq14 s8cq15 s8cq16a s8cq16b s8cq16c s8cq17 s8cq18a1 s8cq18a2 s8cq18b1 s8cq18c1 s8cq18c2 s8cq18d1 s8cq18d2 s8cq17_r s8cq18d2_r s8cq18a2_r s8cq18b2_r s8cq18c2_r _merge s8cq18b2
gen unid = string(hhold) + string(hhid)
save HH_SEC_8C__updated, replace

//HH_SEC_6A  --  Slum area
use "C:\Users\USER\Downloads\HIES2016 24_7_2019-20230807T192137Z-001\HIES2016 24_7_2019\HH_SEC_6A.dta", clear
keep hhold s6aq25
rename s6aq25 slum_area
save HH_6A_slum_area, replace

//hh_sec_a  --  HouseHold information
use "C:\Users\USER\Downloads\HIES2016 24_7_2019-20230807T192137Z-001\HIES2016 24_7_2019\hh_sec_a", clear 
keep hhold hhid id_01_name s1aq01 s1aq02 s1aq03 s1aq04 s1aq05 s1aq07 s2aq04 s7aq01
gen unid = string(hhold) + string(hhid)
save HH, replace


*** // merging all the HH level data

merge 1:1 hhold using "C:\Users\USER\OneDrive\Documents\487-paper\HH_6A_slum_area.dta"
drop if _merge == 2 | _merge == 1
drop _merge
merge 1:m unid using "C:\Users\USER\OneDrive\Documents\487-paper\HH_SEC_1A__updated.dta"
drop if _merge == 2 | _merge == 1
drop _merge
merge 1:1 unid using "C:\Users\USER\OneDrive\Documents\487-paper\HH_SEC_2A__updated.dta"
drop if _merge == 2 | _merge == 1
drop _merge
merge 1:1 unid using "C:\Users\USER\OneDrive\Documents\487-paper\HH_SEC_2B__updated.dta"
drop if _merge == 2 | _merge == 1
drop _merge
merge 1:1 unid using "C:\Users\USER\OneDrive\Documents\487-paper\HH_SEC_6B__updated.dta"
drop if _merge == 2 | _merge == 1
drop _merge
merge 1:1 unid using "C:\Users\USER\OneDrive\Documents\487-paper\HH_SEC_7A__updated.dta"
drop if _merge == 2 | _merge == 1
drop _merge
merge 1:1 unid using "C:\Users\USER\OneDrive\Documents\487-paper\HH_SEC_8C__updated.dta"
drop if _merge == 2 | _merge == 1
drop _merge

// just keeping the HH head observations and renaming the variables
drop if s1aq02 != 1

rename s1aq05 maritial_status
rename s2aq04 edulev_HHead
rename s7aq01 cult_land
rename s1aq03 age
sort hhold
save HH, replace



***------------------------------------------------------------------------***

// -------------     Collapsing some variables by hhold    -----------------// 

***------------------------------------------------------------------------***



// HH maximum education

use "C:\Users\USER\OneDrive\Documents\487-paper\HH_SEC_2A.dta" , clear
collapse (max) s2aq04 , by(hhold)
save highest_household_Education, replace 


****//---   Calculating total income in different section by hhold   ---\\****
//==========================================================================//


// Social safety net program benefit
use "C:\Users\USER\Downloads\HIES2016 24_7_2019-20230807T192137Z-001\HIES2016 24_7_2019\HH_SEC_1C.dta", clear
replace s1cq10a=0 if s1cq10a==.
replace s1cq101d=0 if s1cq101d==.
replace s1cq102d=0 if s1cq102d==.
gen ssnincome= s1cq10a+ s1cq101d+ s1cq102d // Total SSN benefit
save HH_SEC_1C_m.dta, replace
collapse (sum) ssnincome, by (hhold)
save "ssnben.dta", replace

// Wage employment
use "C:\Users\USER\Downloads\HIES2016 24_7_2019-20230807T192137Z-001\HIES2016 24_7_2019\HH_SEC_4B.dta", clear
drop if s4bq00==.
drop if s4bq01==0
rename s4bq00 s4aq00
rename s4bq0a s4aq0a
save "HH_SEC_4B_updated.dta", replace
use "C:\Users\USER\Downloads\HIES2016 24_7_2019-20230807T192137Z-001\HIES2016 24_7_2019\HH_SEC_4A.dta", clear
merge m:n hhold hhid s4aq0a s4aq00 using HH_SEC_4B_updated.dta, nogen
replace s4aq02=0 if s4aq02==.
replace s4aq03=0 if s4aq03==.
replace s4bq02c=0 if s4bq02c==.
replace s4bq05b=0 if s4bq05b==.
replace s4bq08=0 if s4bq08==.
replace s4bq09=0 if s4bq09==.
gen swincome= (s4aq02*s4aq03*s4bq02c)+ (s4aq02*s4aq03*s4bq05b)+ (s4bq08*s4aq02)+ s4bq09 
collapse (sum) swincome, by (hhold)
save salary_wage_inc.dta, replace

// Non-Agricultural Enterprise Income
use "C:\Users\USER\Downloads\HIES2016 24_7_2019-20230807T192137Z-001\HIES2016 24_7_2019\HH_SEC_05.dta", clear
gen naeincome= s5q07* s5q20/100
save "HH_SEC_05_m.dta", replace
collapse (sum) naeincome, by(hhold)
save "naeinc.dta", replace


****//---   Calculating Agriculture Income   ---\\****

// Crop revenue
use "C:\Users\USER\Downloads\HIES2016 24_7_2019-20230807T192137Z-001\HIES2016 24_7_2019\HH_SEC_7B.dta", clear
save HH_SEC_7B__updated.dta, replace
replace s7bq04b=0 if s7bq04b==.
replace s7bq05=0 if s7bq05==.
replace s7bq06=0 if s7bq06==.
replace s7bq07=0 if s7bq07==.
gen croprev= (s7bq05+ s7bq06+ s7bq07)*s7bq04b // crop revenue = (consumption+ sell+ net stock)* Unit price
collapse (sum) croprev, by (hhold)
save croprev.dta, replace

// livestock & poultry revenue
use "C:\Users\USER\Downloads\HIES2016 24_7_2019-20230807T192137Z-001\HIES2016 24_7_2019\HH_SEC_7C2.dta", clear
save HH_SEC_7C2_m.dta, replace
drop if s7c2q00==220
replace s7c2q07b=0 if s7c2q07b==.
replace s7c2q08b=0 if s7c2q08b==.
gen lsprorev= s7c2q07b+ s7c2q08b
save "HH_SEC_7C2__updated.dta", replace
collapse (sum) lsprorev, by( hhold)
save lsprorev.dta, replace

***fishing & fish farming revenue (Section 7C3)***
use "C:\Users\USER\Downloads\HIES2016 24_7_2019-20230807T192137Z-001\HIES2016 24_7_2019\HH_SEC_7C3.dta", clear
save HH_SEC_7C3_m.dta, replace
drop if s7c3q00==230
replace s7c3q11b=0 if s7c3q11b==.
replace s7c3q12b=0 if s7c3q12b==.
gen fishrev= s7c3q11b+ s7c3q12b // fishing revenue= sell+ consume
save "HH_SEC_7C3_mp.dta", replace
collapse (sum) fishrev, by( hhold)
save fishrev.dta, replace // house hold fishing & fish firming revenue

***forestry revenue (Section 7C4)***
use "C:\Users\USER\Downloads\HIES2016 24_7_2019-20230807T192137Z-001\HIES2016 24_7_2019\HH_SEC_7C4.dta", clear
save HH_SEC_7C4_m.dta, replace
drop if s7c4q00==240
replace s7c4q15=0 if s7c4q15==.
replace s7c4q16=0 if s7c4q16==.
gen forrev= s7c4q15+ s7c4q16  // forestry revenue= sell+ consume
save "HH_SEC_7C4_mp.dta", replace
collapse (sum) forrev, by( hhold)
save forrev.dta, replace // house hold forestry revenue

***Expenses on agricultural input (Section 7D)***
use "C:\Users\USER\Downloads\HIES2016 24_7_2019-20230807T192137Z-001\HIES2016 24_7_2019\HH_SEC_7D.dta", clear
save HH_SEC_7D_m.dta, replace
drop if s7dq00==330
gen agroinexp=s7dq02b
replace agroinexp=0 if agroinexp==.
save "HH_SEC_7D_m.dta", replace
collapse (sum) agroinexp, by(hhold)
save agroinexp.dta, replace

// Income/rent from agricultural assets (Section 7E)
use "C:\Users\USER\Downloads\HIES2016 24_7_2019-20230807T192137Z-001\HIES2016 24_7_2019\HH_SEC_7E.dta", clear
save HH_SEC_7E_m.dta, replace
drop if s7eq00==420
gen aaincome=s7eq04
drop if s7eq00==420
replace aaincome=0 if aaincome==.
save "HH_SEC_7E_m.dta", replace
collapse (sum) aaincome, by(hhold)
save aaincome.dta, replace //house hold agricultural asset income/rent

// Energy Consumption (Section 9C)
use "C:\Users\USER\Downloads\HIES2016 24_7_2019-20230807T192137Z-001\HIES2016 24_7_2019\HH_SEC_9C.dta", clear
drop psu hhid s9cq00 s9cq03 s9cq02
save HH_SEC_9C_m.dta, replace
collapse (sum) s9cq01, by(hhold)
rename s9cq01 energy_exp
save energy_exp.dta, replace //house hold energy consumption expenditure

//Health expenditure (Section 1C)
use "C:\Users\USER\Downloads\HIES2016 24_7_2019-20230807T192137Z-001\HIES2016 24_7_2019\HH_SEC_3A.dta", clear
keep hhold s3aq20j s3aq14e
collapse (sum) s3aq14e (sum) s3aq20j, by( hhold)
gen tothelexp = s3aq14e + s3aq20j
drop s3aq14e s3aq20j
save health_exp, replace //house hold total health expenditure

***Agricultual income calculation***
use "C:\Users\USER\Downloads\HIES2016 24_7_2019-20230807T192137Z-001\HIES2016 24_7_2019\hh_sec_a", clear
save agroinc.dta, replace
merge 1:1 hhold using croprev.dta, nogen
merge 1:1 hhold using lsprorev.dta, nogen
merge 1:1 hhold using fishrev.dta, nogen
merge 1:1 hhold using forrev.dta, nogen
merge 1:1 hhold using agroinexp.dta, nogen
merge 1:1 hhold using aaincome.dta, nogen
gen agroincome= croprev+ lsprorev+ fishrev+ forrev- agroinexp+ aaincome 
save agroincex.dta, replace 
keep hhold agroincome
save "agroinc.dta", replace //house hold wise agricultural income


// Income from others assets (Section 8B)
use "C:\Users\USER\Downloads\HIES2016 24_7_2019-20230807T192137Z-001\HIES2016 24_7_2019\HH_SEC_8B.dta", clear
replace s8bq01=0 if s8bq01==.
replace s8bq02=0 if s8bq02==.
replace s8bq03a=0 if s8bq03a==.
replace s8bq03b=0 if s8bq03b==.
replace s8bq03c=0 if s8bq03c==.
replace s8bq04=0 if s8bq04==.
replace s8bq05=0 if s8bq05==.
replace s8bq06=0 if s8bq06==.
replace s8bq07=0 if s8bq07==.
replace s8bq08=0 if s8bq08==.
replace s8bq09=0 if s8bq09==.
replace s8bq11=0 if s8bq11==.
replace s8bq12=0 if s8bq12==.
replace s8bq13=0 if s8bq13==.
gen otrastincome= s8bq01+ s8bq02+ s8bq03a+ s8bq03b+ s8bq03c+ s8bq04+ s8bq05+ s8bq06+ s8bq07+ s8bq08+ s8bq09+ s8bq11+ s8bq12+ s8bq13 //total income from differnt other asset category
save "HH_SEC_8B_m.dta", replace
keep hhold otrastincome
save otrastinc.dta, replace // House hold income from other assets

***House hold Income Calulation***
use "C:\Users\USER\Downloads\HIES2016 24_7_2019-20230807T192137Z-001\HIES2016 24_7_2019\hh_sec_a", clear
save hholdid.dta, replace
merge 1:1 hhold using ssnben.dta, nogen
merge 1:1 hhold using swinc.dta, nogen
merge 1:1 hhold using naeinc.dta, nogen
merge 1:1 hhold using agroinc.dta, nogen
merge 1:1 hhold using otrastinc.dta, nogen
gen hhincome= ssnincome+ swincome+ naeincome+ agroincome+ otrastincome //all income from section 1, 4, 5, 7 & 8
save hholdincome.dta, replace
keep hhold hhincome
save hhinc.dta, replace // house hold income summary file
rename hhincome Total_Income
save Total_Income, replace 


****//--------------------   Cleaning 6B datasest   --------------------\\****
//==========================================================================//

use "C:\Users\USER\OneDrive\Documents\487-paper\HH_SEC_6B__updated", clear
drop hhid s6bq03 s6bq04 s6bq05a s6bq05b s6bq05c s6bq05d unid id_01_name
drop if s6bq00 != 101
drop s6bq00
rename s6bq02 drought
replace drought = 0 if drought == 2
label variable drought "experience in drought shock in the past 12 months"
save 6B_drought, replace

use "C:\Users\USER\OneDrive\Documents\487-paper\HH_SEC_6B__updated", clear
drop hhid s6bq03 s6bq04 s6bq05a s6bq05b s6bq05c s6bq05d unid id_01_name
drop if s6bq00 != 102
drop s6bq00
rename s6bq02 floods
replace floods = 0 if floods == 2
label variable floods "experience in floods shock in the past 12 months"
save 6B_floods, replace

use "C:\Users\USER\OneDrive\Documents\487-paper\HH_SEC_6B__updated", clear
drop hhid s6bq03 s6bq04 s6bq05a s6bq05b s6bq05c s6bq05d unid id_01_name
drop if s6bq00 != 103
drop s6bq00
rename s6bq02 erosion
replace erosion = 0 if erosion == 2
label variable erosion "experience in erosion shock in the past 12 months"
save 6B_erosion, replace

use "C:\Users\USER\OneDrive\Documents\487-paper\HH_SEC_6B__updated", clear
drop hhid s6bq03 s6bq04 s6bq05a s6bq05b s6bq05c s6bq05d unid id_01_name
drop if s6bq00 != 118
drop s6bq00
rename s6bq02 torEQFire
replace torEQFire = 0 if torEQFire == 2
label variable drought "experience in tornado/EQ/Fire shock in the past 12 months"
save 6B_torEQFire, replace
save 6B_erosion, replace



****//--------------------   Cleaning 9C datasest   --------------------\\****
//==========================================================================//

use "C:\Users\USER\Downloads\HIES2016 24_7_2019-20230807T192137Z-001\HIES2016 24_7_2019\HH_SEC_9C.dta", clear
keep if s9cq00 == 241 | s9cq00 == 242 | s9cq00 == 243 | s9cq00 == 244 | s9cq00 == 245 | s9cq00 == 246 | s9cq00 == 247 | s9cq00 == 248 | s9cq00 == 249  //just keeping HH energy consumption data

* reshaping data from long to wide format by hhold
sort hhold hhid s9cq00
drop psu hhid id s9cq02 s9cq03
gen jn = s9cq00
reshape wide s9cq00 s9cq01, i(hhold) j( jn )
drop jn

// renaming and cleaning the variables
rename s9cq01241 firewood
rename s9cq01242 cow_dung_bhushi
rename s9cq01243 jute_stick
rename s9cq01244 kerosene
rename s9cq01245 agri_products_fuel
rename s9cq01246 gas_LPG
rename s9cq01247 electricity
rename s9cq01248 coal
rename s9cq01249 other_fuels-candles
rename s9cq01249 other_fuels__candles
drop s9cq00241 s9cq00242 s9cq00243 s9cq00244 s9cq00245 s9cq00246 s9cq00247 s9cq00248 s9cq00249
replace firewood = 0 if firewood == .
replace cow_dung_bhushi = 0 if cow_dung_bhushi == .
replace jute_stick = 0 if jute_stick == .
replace kerosene = 0 if kerosene == .
replace agri_products_fuel = 0 if agri_products_fuel == .
replace gas_LPG = 0 if gas_LPG == .
replace electricity = 0 if electricity == .
replace coal = 0 if coal == .
replace other_fuels__candles = 0 if other_fuels__candles == .

duplicates report hhold
duplicates drop hhold, force

gen totEnergyCons = firewood + cow_dung_bhushi + jute_stick + kerosene + agri_products_fuel + gas_LPG + electricity + coal + other_fuels__candles
save 9C_wide_format_energy_cons, replace


***------------------------------------------------------------------------***

// -------------          Merging all the dataset          -----------------// 

***------------------------------------------------------------------------***


use "C:\Users\USER\OneDrive\Documents\487-paper\HH", clear

merge 1:1 hhold using "C:\Users\USER\OneDrive\Documents\487-paper\Total_Income" 
drop if _merge == 2
drop _merge
save income_floods, replace
merge 1:1 hhold using "C:\Users\USER\OneDrive\Documents\487-paper\croprev"
drop if _merge == 2
drop _merge
merge 1:1 hhold using "C:\Users\USER\OneDrive\Documents\487-paper\highest_household_education"
drop if _merge == 2
drop _merge
merge 1:1 hhold using "C:\Users\USER\OneDrive\Documents\487-paper\ssninc"
drop if _merge == 2
drop _merge
merge 1:1 hhold using "C:\Users\USER\OneDrive\Documents\487-paper\health_exp"
drop if _merge == 2
drop _merge
merge 1:1 hhold using "C:\Users\USER\OneDrive\Documents\487-paper\6B_dought"
drop if _merge == 2
drop _merge
merge 1:1 hhold using "C:\Users\USER\OneDrive\Documents\487-paper\6B_floods"
drop if _merge == 2
drop _merge
merge 1:1 hhold using "C:\Users\USER\OneDrive\Documents\487-paper\6B_erosion"
drop if _merge == 2
drop _merge
merge 1:1 hhold using "C:\Users\USER\OneDrive\Documents\487-paper\6B_torEQFire"
drop if _merge == 2
drop _merge
merge 1:1 hhold using "C:\Users\USER\OneDrive\Documents\487-paper\6B_tornado_shock"
drop if _merge == 2
drop _merge
merge 1:1 hhold using "C:\Users\USER\OneDrive\Documents\487-paper\highest_household_Education"
drop if _merge == 2
drop _merge
merge 1:1 hhold using "C:\Users\USER\OneDrive\Documents\487-paper\HH_members"
drop if _merge == 2
drop _merge
merge 1:1 hhold using "C:\Users\USER\OneDrive\Documents\487-paper\9C_wide_format_energy_cons"
drop if _merge == 2
drop _merge


*** cleaning dataset -- renameing and labeling

label variable croprev "Revenue from crop production"
label variable ssnincome "Total income from Social Safety Net Programs"
label variable agroinexp "Agricuriculture input expenditure"
label variable hh_members "Total members in the HH"
label variable maxedu_HH "Maximum education in the HH"
label variable edulev_HHead "HH head education level"
label variable totEnergyCons "Total energy consumption in the last 12 months"
label variable Total_Income "Total HH income"
label variable Total_Income "Total HH income/annually"
label variable floods "experience in flood shock in the past 12 months"
label variable firewood "Firewood consumption"
label variable cow_dung_bhushi "Cow dung, cakes, bhushi and wood-powder consumption"
label variable jute_stick "Jute stick consumption"
label variable kerosene "Kerosene consumption"
label variable agri_products_fuel "Agri by products fuel: paddy, hag, pressed sugarcane and dried com plants, etc."
label variable gas_LPG "Gas (natural, bio-gas)/LPG"
label variable electricity "Electricity consumption"
label variable coal "Coal consumption"
label variable coal "Pit coal, char coal, wood coal consumption"
label variable agri_products_fuel "Agri by products fuel: paddy, hag, pressed sugarcane and dried com plants, etc. consumption"
label variable gas_LPG "Gas (natural, bio-gas)/LPG consumption"
label variable other_fuels__candles "Other fuels and light (e.g. matches and candles etc.) consumption"


***------------------------------------------------------------------------***

// -------------              Analysis                     -----------------// 

***------------------------------------------------------------------------***

* Descriptive analysis

mrtab firewood- other_fuels__candles, by(division) row nofreq

tabstat firewood - other_fuels__candles, by(floods) stat(mean sd)
tabstat firewood - other_fuels__candles, by(drought) stat(mean sd)
tabstat firewood - other_fuels__candles, by(erosion) stat(mean sd)
tabstat firewood - other_fuels__candles, by(torEQFire) stat(mean sd)


* Linear regression model for estimating energy consumption based on household characteristics.
//-------------------------------------------------------------------------------------------//
reg lnEnergyCons i.division hh_members maxedu_HH cult_land slum_area lnIncome

estat hettest
vif
estat ovtest

// presence of heteroskedasticity -- using robust standard error
reg lnEnergyCons i.division hh_members maxedu_HH cult_land slum_area lnIncome, robust


* Probit regression model for identifying areas more likely to be affected by disasters.
//-------------------------------------------------------------------------------------------//

probit floods i.division_id slum_area
probit drought i.division_id slum_area
probit erosion i.division_id slum_area
probit torEQFire i.division_id slum_area

* Propensity Score Matching for finding the difference in energy consumption between households affected by disasters and those unaffected
//-------------------------------------------------------------------------------------------//

psmatch2 nd_shock i.division hh_members maxedu_HH cult_land slum_area lnIncome, logit outcome(totEnergyCons)
psgraph
twoway (kdensity _pscore if _treated==1) (kdensity _pscore if _treated==0, lpattern(dash)), legend( label( 1 "affected") label( 2 "unaffected") ) xtitle("propensity score")

// ATE
teffects psmatch (totEnergyCons) (floods i.division hh_members maxedu_HH cult_land slum_area lnIncome, probit)
teffects psmatch (totEnergyCons) (drought i.division hh_members maxedu_HH cult_land slum_area lnIncome, probit)
teffects psmatch (totEnergyCons) (erosion i.division hh_members maxedu_HH cult_land slum_area lnIncome, probit)
teffects psmatch (totEnergyCons) (torEQ_fire i.division hh_members maxedu_HH cult_land slum_area lnIncome, probit)

// ATET
teffects psmatch (totEnergyCons) (floods i.division hh_members maxedu_HH cult_land slum_area lnIncome, probit), atet
teffects psmatch (totEnergyCons) ( drought i.division hh_members maxedu_HH cult_land slum_area lnIncome, probit), atet
teffects psmatch (totEnergyCons) ( erosion i.division hh_members maxedu_HH cult_land slum_area lnIncome, probit), atet
teffects psmatch (totEnergyCons) ( torEQ_fire i.division hh_members maxedu_HH cult_land slum_area lnIncome, probit), atet













