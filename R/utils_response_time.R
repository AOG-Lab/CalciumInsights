response_time <- function(data, peak, Puntos_medios){
  filtered_firts <- dplyr::filter(data, data[,1] <= data[,1][peak[1,2]], data[,2]>= Puntos_medios[1,2], peak[1,1] >= data[,2])

  ultimo_peak <- length(peak[,1])

  filtered_second <- dplyr::filter(data, data[,1] >= data[,1][peak[ultimo_peak,2]], data[,2] >= Puntos_medios[ultimo_peak,2], peak[ultimo_peak,1] >= data[,2])


  MCO <- function(M1,M2){
    W <- 1000
    delta_1 <- W*matrix(1, length(M1), length(M1))
    LM1 <- 1:length(M1)
    for(i in LM1){
      for(j in LM1){
        if(M1[i] < M1[j]) delta_1[i,j] <- -1
        if(M1[i] == M1[j]) delta_1[i,j] <- 0
      }
    }
    delta_2 = W*matrix(1, length(M2), length(M2))
    LM2 <- 1:length(M2)
    for(i in LM2){
      for(j in LM2){
        if(M2[i] < M2[j]) delta_2[i,j] <- -1
        if(M2[i] == M2[j]) delta_2[i,j] <- 0
      }
    }
    alpha <- delta_1+delta_2
    Gamma_m1 <- replace(alpha,alpha==0 | alpha==W,W)
    Gamma_m2 <- replace(Gamma_m1,Gamma_m1==2*W ,2*W)
    Gamma_m <- replace(Gamma_m2,Gamma_m2 != 2*W & Gamma_m2 != W,0)
    Beta <- apply(Gamma_m, 1, sum)
    B_selec <- Beta<2*W
    grep('TRUE', B_selec, value=FALSE)
  }

  filtered_second_min = MCO(M1 = filtered_second[,1] ,M2 = filtered_second[,2] )

  filtered_second1 <- filtered_second[filtered_second_min,]


  # Encontrar los valores que minimizan ambas variables simultáneamente
  first_time <- filtered_firts[which.min(filtered_firts[,1] + filtered_firts[,2]), ]
  second_time <- filtered_second1[which.min(filtered_second1[,2]), ]

  Tiempo_respose <- second_time[,1] - first_time[,1]

  return(list(first_time = first_time, second_time = second_time,  Tiempo_respose = Tiempo_respose))
}
