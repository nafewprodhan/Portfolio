
**************Variable Label***************

label variable A3 "Site of Interview (Camp No)"
label variable B1 "Age"
label variable B2 "Sex"
label variable B3 "Do you have difficulty in seeing, even if wearing glasses?"
label variable B4 "Do you have difficulty in hearing, even if using a hearing aid(s)?"
label variable B5 "Do you have difficulty in walking or climbing steps?"
label variable B6 "Do you have difficulty remembering or concentrating?"
label variable B7 "Do you have difficulty with self-care, such as washing all over or dressing?"
label variable B8 "Using your usual language, do you have difficulty in communicating, for example understanding or being understood?"
label variable B9 "Marital Status"
label variable B10_a "GBV"
label variable B10_b "SRH"
label variable B10_c "WASH"
label variable B10_d "DRR"
label variable B11_a "life-skills activities"
label variable B11_b "Shine activities"
label variable B11_c "life-skill and leadership training"
label variable B11_d "life skill training for youths"
label variable B11_e "Others (specify)"
label variable B11_f "None of the above"
label variable C1 "Which following are the forms of GBV"
label variable C2_a "Psychosocial support & counselling "
label variable C2_b "Training on EMAP tools"
label variable C2_c "training on Gender Based Violence Pocket Guide (GBVPG)"
label variable C2_d "Training on GBV prevention & protection"
label variable C2_e "information about GBV prevention & mitigation"
label variable C2_f "Others (specify)"
label variable C3 "Did you received services from WGSS (Women-Girl Safe Space)"
label variable C4 "Do your husband/in-laws/ family members allow you to receive services from WGSS?"
label variable C5 "Do you allow your wife or female family members to receive services from WGSS?"
label variable C6 "Did any community groups have been created and trained up for GBV prevention?"
label variable C7 "Have you ever been a victim of violence in your personal life?"
label variable C10 "In last three years did you experience or see any GBV occurrences in the camp?"
label variable C11_a "Applied the knowledge received from CARE’s sessions"
label variable C11_b "Contacted with community leaders (imams, majhi)"
label variable C11_c "Contacted to CARE outreach member"
label variable C11_d "None of the above"
label variable C12_a "Apply the knowledge received from CARE’s sessions"
label variable C12_b "Contact with community leaders (imams, majhi)"
label variable C12_c "Contact to CARE outreach member"
label variable C12_d "None of the above"
label variable D1_a "ANC"
label variable D1_b "PNC"
label variable D1_c "Family Planning method"
label variable D1_d "Refferal services for Safe delivery"
label variable D1_e "STI/RTI"
label variable D1_f "Menstrual hygiene"
label variable D1_g "Awareness session/information on reproductive health"
label variable D1_h "Couple counselling"
label variable D1_i "None of the above"
label variable D2_a "CARE’s Health Post"
label variable D2_b "Outreach clinics"
label variable D2_c "Women Girl Safe Space"
label variable D2_d "Court-yard meetings"
label variable D2_e "Door-to-door visit"
label variable D3 "Do you decide to take the SRH services yourself or do you need permission from someone in your family?"
label variable D4_a "From my husband"
label variable D4_b "From my mother"
label variable D4_c "From mother-in-law"
label variable D4_d "Other family members"
label variable D5_a "ANC"
label variable D5_b "PNC"
label variable D5_c "Family Planning method"
label variable D5_d "Referral services for Safe delivery"
label variable D5_e "STI/RTI"
label variable D5_f "Menstrual hygiene"
label variable D5_g "Awareness session/information on reproductive health"
label variable D5_h "Couple counselling"
label variable D6 "Which of the following contraceptive methods did you use in the past for family planning?"
label variable D7 "What is your current contraceptive method?"
label variable D8 "Who makes the decision to receive family planning and other sexual health services?"
label variable D11 "Do you think everyone should know about sexual health regardless of men, women, boys, and girls?"
label variable E1 "Do you know what to do if there is a disaster and whom to contact for assistance?"
label variable E2 "Did you participate in any DRR activities as a cash-for-work laborer to construct/repair targeted structures?"
label variable E3_a "slope"
label variable E3_b "drainage system"
label variable E3_c "Durable pathways"
label variable E3_d "staircase with railing"
label variable E3_e "All of the above"
label variable E4 "Did you participated in Block development committee meetings?"
label variable E5 "Did you help CARE team to identify the sites which are vulnerable to disasters?"
label variable E6 "Do you participate in disaster campaigns, or disaster awareness or training sessions?"
label variable F1_a "good sanitation facility"
label variable F1_b "adequate safe drinking water"
label variable F1_c "gender-segregated latrines and bathing cubicles"
label variable F1_d "All of the above"
label variable F1_e "None of the above"
label variable F2_a "Deep tube wells"
label variable F2_b "Shallow tube wells"
label variable F2_c "Tap stand/Pipe water"
label variable F2_d "Water tank"
label variable F2_e "Protected dug well"
label variable F2_f "Unprotected dug well"
label variable F2_g "Others (specify)"
label variable F3 "Do you use latrines in a proper way regularly?"
label variable F4 "Are the WASH facilities have easy access for people with disability?"
label variable F5_a "Availability of safe drinking water"
label variable F5_b "People use latrines regularly"
label variable F5_c "People maintain personal hygiene"
label variable F5_d "Others (specify)"
label variable F5_e "All of the above"
label variable F5_f "None of the above"
label variable G "How out of five would you rate CARE's services in improving health, safety and infrastructure?"

**************Value Label***************

lab def Age 1"10-19" 2"20-49" 3"50 and above"
lab val B1 Age

lab def Sex 1"Male" 2"Female" 3"Others"
lab val B2 Sex

lab def DifSee 1"No difficulty" 2"Some difficulty" 3"A lot of difficulty" 4"Cannot do at all" 5"Refused" 6"Don’t know"
lab val B3 DifSee

lab def DifHear 1"No difficulty" 2"Some difficulty" 3"A lot of difficulty" 4"Cannot do at all" 5"Refused" 6"Don’t know"
lab val B4 DifHear

lab def DifWalk 1"No difficulty" 2"Some difficulty" 3"A lot of difficulty" 4"Cannot do at all" 5"Refused" 6"Don’t know"
lab val B5 DifWalk

lab def DifRemem 1"No difficulty" 2"Some difficulty" 3"A lot of difficulty" 4"Cannot do at all" 5"Refused" 6"Don’t know"
lab val B6 DifRemem

lab def DifComm 1"No difficulty" 2"Some difficulty" 3"A lot of difficulty" 4"Cannot do at all" 5"Refused" 6"Don’t know"
lab val B8 DifComm

lab def Marrital 1"Married" 2"Unmarried" 3"Widow/widower" 4"Divorced"
lab val B9 Marrital

lab def services_a 1"Yes" 2"No"
lab val B10_a services_a

lab def services_b 1"Yes" 2"No"
lab val B10_b services_b

lab def services_c 1"Yes" 2"No"
lab val B10_c services_c

lab def training_a 1"Yes" 2"No"
lab val B11_a training_a

lab def training_b 1"Yes" 2"No"
lab val B11_b training_b

lab def training_c 1"Yes" 2"No"
lab val B11_c training_c

lab def training_d 1"Yes" 2"No"
lab val B11_d training_d

lab def training_e 1"Yes" 2"No"
lab val B11_e training_e

lab def training_f 1"Yes" 2"No"
lab val B11_f training_f

lab def GBV 1"Intimate Partner Violence" 2"restriction of access to services and information" 3"Rape" 4"Early/forced marriage" 5"Physical assault" 6"Verbal/emotional abuse" 7"Others (specify)" 8"All of the above" 9"Don’t Know"
lab val C1 GBV

lab def Services_a 1"Yes" 2"No"
lab val C2_a training_a

lab def Services_b 1"Yes" 2"No"
lab val C2_b training_b

lab def Services_c 1"Yes" 2"No"
lab val C2_c training_c

lab def Services_d 1"Yes" 2"No"
lab val C2_d training_d

lab def Services_e 1"Yes" 2"No"
lab val C2_e training_e

lab def Services_f 1"Yes" 2"No"
lab val C2_f training_f

lab def WGSS 1"Yes" 2"No"
lab val C3 WGSS

lab def HusWGSS 1"Yes" 2"No"
lab val C4 HusWGSS

lab def WifWGSS 1"Yes" 2"No"
lab val C5 WifWGSS

lab def CommGBV 1"Yes" 2"No"
lab val C6 CommGBV

lab def victim 1"Yes" 2"No"
lab val C7 victim

lab def ExpGBV 1"Yes" 2"No"
lab val C10 ExpGBV

lab def Actions_a 1"Yes" 2"No"
lab val C11_a Actions_a

lab def Actions_b 1"Yes" 2"No"
lab val C11_b Actions_b

lab def Actions_c 1"Yes" 2"No"
lab val C11_c Actions_c

lab def Actions_d 1"Yes" 2"No"
lab val C11_d Actions_d

lab def DoGBV_a 1"Yes" 2"No"
lab val C12_a DoGBV_a

lab def DoGBV_b 1"Yes" 2"No"
lab val C12_b DoGBV_b

lab def DoGBV_c 1"Yes" 2"No"
lab val C12_c DoGBV_c

lab def DoGBV_d 1"Yes" 2"No"
lab val C12_d DoGBV_d

lab def SRH_a 1"Yes" 2"No"
lab val D1_a SRH_a

lab def SRH_b 1"Yes" 2"No"
lab val D1_b SRH_b

lab def SRH_c 1"Yes" 2"No"
lab val D1_c SRH_c

lab def SRH_d 1"Yes" 2"No"
lab val D1_d SRH_d

lab def SRH_e 1"Yes" 2"No"
lab val D1_e SRH_e

lab def SRH_f 1"Yes" 2"No"
lab val D1_f SRH_f

lab def SRH_g 1"Yes" 2"No"
lab val D1_g SRH_g

lab def SRH_h 1"Yes" 2"No"
lab val D1_h SRH_h

lab def SRH_i 1"Yes" 2"No"
lab val D1_i SRH_i

lab def SRHfrom_a 1"Yes" 2"No"
lab val D2_a SRHfrom_a

lab def SRHfrom_b 1"Yes" 2"No"
lab val D2_b SRHfrom_b

lab def SRHfrom_c 1"Yes" 2"No"
lab val D2_c SRHfrom_c

lab def SRHfrom_d 1"Yes" 2"No"
lab val D2_d SRHfrom_d

lab def SRHfrom_e 1"Yes" 2"No"
lab val D2_e SRHfrom_e

lab def decide 1"I take my decisions" 2"I need permission"
lab val D3 decide

lab def whom_a 1"Yes" 2"No"
lab val D4_a whom_a

lab def whom_b 1"Yes" 2"No"
lab val D4_b whom_b

lab def whom_c 1"Yes" 2"No"
lab val D4_c whom_c

lab def whom_d 1"Yes" 2"No"
lab val D4_d whom_d

lab def SRHreq_a 1"Yes" 2"No"
lab val D5_a SRHreq_a

lab def SRHreq_b 1"Yes" 2"No"
lab val D5_b SRHreq_b

lab def SRHreq_c 1"Yes" 2"No"
lab val D5_c SRHreq_c

lab def SRHreq_d 1"Yes" 2"No"
lab val D5_d SRHreq_d

lab def SRHreq_e 1"Yes" 2"No"
lab val D5_e SRHreq_e

lab def SRHreq_f 1"Yes" 2"No"
lab val D5_f SRHreq_f

lab def SRHreq_g 1"Yes" 2"No"
lab val D5_g SRHreq_g

lab def SRHreq_h 1"Yes" 2"No"
lab val D5_h SRHreq_h

lab def Contra 1"Condom" 2"Birth Control Pills" 3"Depo Injection" 4"Implant" 5"Others (specify)" 6"Didn’t use any in the past"
lab val D6 Contra

lab def CurrContra 1"Condom" 2"Birth Control Pills" 3"Depo Injection" 4"Implant" 5"Others (specify)" 6"Didn’t use any method currently"
lab val D7 CurrContra

lab def MakDeci 1"Jointly (Husband & wife)" 2"Only my partner (husband or wife)" 3"Myself" 4"Don’t wish to reply" 5"Others (specify)"
lab val D8 MakDeci

lab def everyone 1"Yes" 2"No"
lab val D11 everyone

lab def disaster 1"Yes" 2"No"
lab val E1 disaster

lab def DRR 1"Yes" 2"No"
lab val E2 DRR

lab def infras_a 1"Yes" 2"No"
lab val E3_a infras_a

lab def infras_b 1"Yes" 2"No"
lab val E3_b infras_b

lab def infras_c 1"Yes" 2"No"
lab val E3_c infras_c

lab def infras_d 1"Yes" 2"No"
lab val E3_d infras_d

lab def infras_e 1"Yes" 2"No"
lab val E3_e infras_e

lab def Block 1"Yes" 2"No"
lab val E4 Block

lab def CAREteam 1"Yes" 2"No"
lab val E5 CAREteam

lab def aware 1"Yes" 2"No"
lab val E6 aware

lab def WASH_a 1"Yes" 2"No"
lab val F1_a WASH_a

lab def WASH_b 1"Yes" 2"No"
lab val F1_b WASH_b

lab def WASH_c 1"Yes" 2"No"
lab val F1_c WASH_c

lab def WASH_d 1"Yes" 2"No"
lab val F1_d WASH_d

lab def WASH_e 1"Yes" 2"No"
lab val F1_e WASH_e

lab def drinking_a 1"Yes" 2"No"
lab val F2_a drinking_a

lab def drinking_b 1"Yes" 2"No"
lab val F2_b drinking_b

lab def drinking_c 1"Yes" 2"No"
lab val F2_c drinking_c

lab def drinking_d 1"Yes" 2"No"
lab val F2_d drinking_d

lab def drinking_e 1"Yes" 2"No"
lab val F2_e drinking_e

lab def drinking_f 1"Yes" 2"No"
lab val F2_f drinking_f

lab def drinking_g 1"Yes" 2"No"
lab val F2_g drinking_g

lab def latrines 1"Yes" 2"No" 3"I don’t know how to use latrines"
lab val F3 latrines

lab def WASHeasy 1"Yes" 2"No" 3"There are separate infrastructures for them"
lab val F4 WASHeasy

lab def improve_a 1"Yes" 2"No"
lab val F5_a improve_a

lab def improve_b 1"Yes" 2"No"
lab val F5_b improve_b

lab def improve_c 1"Yes" 2"No"
lab val F5_c improve_c

lab def improve_d 1"Yes" 2"No"
lab val F5_d improve_d

lab def improve_e 1"Yes" 2"No"
lab val F5_e improve_e

lab def improve_f 1"Yes" 2"No"
lab val F5_f improve_f

***************Analysis*****************

