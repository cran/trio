grp.palette <-
function(name, width = 500, height =500){
     jpeg(filename = name, width = width, height = height,
          pointsize = 10, quality = 100, bg = "white", res = NA)
     
     
     recL = 1
     colMa = matrix(colors(), ncol = 25, nrow = 27, byrow = T)
     
     my.xlim = c(0, 25)
     my.ylim = c(0, 27)
     plot(my.xlim, my.ylim, type = "n", xlab="", ylab="", axes=F)
     
     
     for(row in 1:27){
       for(col in 1:25){
         rect( (col-1)* recL + .1, (row-1)* recL + .1, col*recL - .1, row*recL - .1, col = colMa[row, col] )
         #print(paste("x=", (col-1)* recL + .1, ", y=", (row-1)* recL + .1))
         #print(paste("col=", col, "row=", row))
         
       }
     
     }
     axis(2, seq(.5,27.5, by=2), seq(0,670, by=50), cex=.5)
     axis(1)
     grid(nx = NULL, ny = NULL, col = colors()[31], lty = "dotted",
          lwd = NULL, equilogs = TRUE)
     
     dev.off()
     return(NULL)
}

