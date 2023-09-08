files = list.files(getwd())

for(theta in c(0.1,0.07, 0.03, 0, -0.03, -0.07, -0.1)) {
    for (propInvalid in c(0.3, 0.5, 0.7)) {
        ind = paste0("theta",theta,"thetaU0.3prop_invalid",propInvalid)
        file = files[grepl(ind,files)]

        
        cat("Finish theta", theta, " Prop",propInvalid, "  N ",length(file),"\n")
    }
}



