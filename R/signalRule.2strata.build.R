signalRule.2strata.build <-
function(sigStr="g9=11 and g13=11", sigType="geno2d", para=c(-5, 1)){

  ## Dec08Change!!!: allow D/R coding
  #if (sigType !="geno1d") stop( paste("sigType=", sigType, " is not implemented.", sep=""))
  if( !is.element( sigType, c("geno2d", "D/R") ) ) stop( paste("Wrong input: sigType=", sigType, ". It is not implemented.", sep=""))
  
  if( length(para)!=2) stop( paste("Wrong length of input para:", length(para), sep=""))
        
  rule = signalRule.contr()
  rule = signalRule.setInter(rule, para[1])
  sig1 = signal.new(para[2], type=sigType, signalStr = sigStr)
  rule = signalRule.addSignal(rule, sig1)

  return (rule)
}

