moveDir <-
function(par=getwd(), fromRoot, toRoot, subDir=NULL){


  fListdir = list.files(path = file.path(par, fromRoot))
  dir.index = file.info(file.path(par, fromRoot, fListdir))$isdir
  dirs = fListdir[dir.index]
  for( j in dirs){
    newdir = file.path(par, toRoot, j)
    dir.create(path=newdir, showWarnings = TRUE)
  }
  ## addtional subdirectory
  if(!is.null(subDir))
    dir.create(path=file.path(par, toRoot, subDir[1], subDir[2]), showWarnings = TRUE)
  
  fList = list.files(path = file.path(par, fromRoot), all.files = T,
           full.names = F, recursive = T)

  for( i in fList){
    f1 = file.path( par, fromRoot, i)
    f2 = file.path( par, toRoot, i)
    print(f1)
    print(f2)
    file.copy(from=f1, to=f2, overwrite=T )

  }
}

