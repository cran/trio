bkMap.spHap <-
function(bkMap, subjectCt){

  len = length(bkMap$bkLens)

  re = NULL
  people = NULL
  people.in = NULL
  for( i.bk in 1:len ){
    curIn = sample(bkMap$bkLens[i.bk], size = subjectCt, replace = T, prob = bkMap$bks[[i.bk]][,bkMap$probCol])

    curExp = lapply(curIn, FUN=function(item, bkMap, i.bk){
                              as.character(bkMap$bks[[i.bk]][item,bkMap$expCol])
                            }, bkMap = bkMap, i.bk=i.bk)
    people = c(people, list(unlist(curExp)))
    people.in = c(people.in, list(curIn))
  }
  subjects = matrix(unlist(people), ncol=length(people), byrow=F)
  subjects.in = matrix(unlist(people.in), ncol=length(people.in), byrow=F)

  re = list(subjects=subjects, subjectsIn=subjects.in)
  return(re)
}

