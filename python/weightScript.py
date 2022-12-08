xSection=3.309e-2
kFactor=2.27
brFraction=3*0.58*0.58*0.0227
nTotal=1.5e6
sigmaBr=xSection*kFactor*brFraction
print("xSection : ",xSection , " fb --> ",xSection*1e-3," pb")
print("kFactor : ",kFactor)
print("brFraction : ",brFraction)
print("sigma x br : ",sigmaBr ," fb --> ",sigmaBr*1e-3," pb")
print("sigma x br / hgg : ",sigmaBr/0.0227, " fb --> ", sigmaBr*1e-3/0.0227 ," pb")

print("expected number of events : " ,sigmaBr," / fb-1" )
print("expected number of events for 57 /fb: " ,sigmaBr*57 )

print("Weight per event for nTotal : ",nTotal/1e6,"M = ",sigmaBr/nTotal)
print("Weight per event for nTotal : ",nTotal/1e6,"M in fgg  = ",sigmaBr*1e-3/nTotal)

