signalRule.addSignal <-
function(rule, signal){

  rule$slist = c(rule$slist, list(signal))
  rule$snpIdx = c(rule$snpIdx, signal$var.idx)
  
  return (rule)
}

