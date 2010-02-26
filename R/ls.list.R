ls.list <-
function(path=NULL, objname,  exam.ext = .3){
  
  myobj = NULL
 
  
  if( (!is.null(path)) & (!is.character(objname))){
    assign(objname, NULL)
    load(path)
    myobj = get(objname)
  }else{
    myobj = objname
  }

  if( !is.list(myobj)){
    stop(paste("Object name (", objname, ") is not a list!", sep=""))
  }
  element = NULL
  l.len = length(myobj)
  names.all = names(myobj)
  for( i in 1:l.len){
    #print(i)
    #print(element)
    elm = myobj[[i]]
    # str(elm)
    if( !is.list(elm)){
      element = c(element, paste("Object", i, ") name=", names.all[i], "; value=", elm, ";", sep=""))
    }else{
      # check structure
      #str(elm[[1]])
      if (length(elm)==1){
        element = c(element, paste("Object", i, ") List::name=", names.all[i], "; length=1", sep=""))
      }else if(length(elm)>1){
        is.allSame = qStr.pattern(elm, exam.ext=exam.ext)
        ##print(is.allSame)
        if( is.null(is.allSame)){
          element = c(element, paste("Object", i, ") List::name=", names.all[i], "; length=", length(elm), ";", sep=""))
        }else{
          ## same objects 
          element = c(element, paste("Object", i, ") List::name=", names.all[i], "; length=", length(elm), ";", sep=""))
      
          element = c(element, paste("            structure::", is.allSame, ";", sep=""))

        }
      }else{
        element = c(element, paste("Object", i, ") List::name=", names.all[i], "; length=0", sep=""))
      }

    }
  }
  return(element)
}

