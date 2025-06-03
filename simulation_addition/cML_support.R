
cML_estimate <- function(b_exp,b_out,
                         se_exp,se_out,
                         K,initial_theta = 0,
                         initial_mu = rep(0,length(b_exp)),
                         maxit = 100)
{
  p = length(b_exp)
  ### initialize
  theta = initial_theta
  theta_old = theta-1
  mu_vec = initial_mu

  ###
  ite_ind = 0
  while( (abs(theta_old - theta) > 1e-7) & (ite_ind<maxit))
  {
    theta_old = theta
    ite_ind = ite_ind + 1

    ### first, update v_bg
    if(K>0)
    {
      v_importance = (b_out - theta*mu_vec)^2 / se_out^2
      nonzero_bg_ind = sort((order(v_importance,decreasing = T))[1:K])
      v_bg = rep(0,p)
      v_bg[nonzero_bg_ind] = (b_out - theta*mu_vec)[nonzero_bg_ind]
    } else{
      v_bg = rep(0,p)
    }


    ### second, update mu_vec
    mu_vec =
      (b_exp / se_exp^2 + theta*(b_out - v_bg) / se_out^2) /
      (1 / se_exp^2 + theta^2 / se_out^2)

    ### third, update theta
    theta =
      sum((b_out - v_bg)*mu_vec / se_out^2) /
      sum(mu_vec^2 / se_out^2)

  }
  ### one more step for v_bg and mu
  if(K>0)
  {
    nonzero_ind = which(v_bg!=0)
    mu_vec[nonzero_ind] = b_exp[nonzero_ind]
    v_bg[nonzero_ind] = (b_out - theta*mu_vec)[nonzero_ind]
  }

  return(list(theta = theta,
              b_vec = mu_vec,
              r_vec = v_bg,ite_ind= ite_ind))

}




cML_SdTheta <- function(b_exp,b_out,
                        se_exp,se_out,
                        theta,b_vec,r_vec)
{
  nonzero_ind = which(r_vec!=0)
  zero_ind = which(r_vec==0)

  VarTheta =
    1/(sum((b_vec^2/se_out^2)[zero_ind])
       -sum(
         (
           (2*b_vec*theta - b_out)^2/se_out^4*
             1/(1/se_exp^2 + theta^2/se_out^2)
         )[zero_ind]
       ) )

  if(VarTheta<=0)
  {
    return(NaN)
  } else {
    return(sqrt(VarTheta))
  }

}



cML_estimate_random <- function(b_exp, b_out,
                                se_exp, se_out,
                                K,random_start = 0,
                                maxit = 100)
{
  p = length(b_exp)
  min_theta_range = min(b_out/b_exp)
  max_theta_range = max(b_out/b_exp)

  theta_v_RandomCandidate = NULL
  sd_v_RandomCandidate = NULL
  l_v_RandomCandidate = NULL
  invalid_RandomCandidate = NULL

  for(random_ind in 1:(1+random_start))
  {
    if(random_ind == 1)
    {
      initial_theta = 0
      initial_mu = rep(0,p)
    } else {
      initial_theta = runif(1,min = min_theta_range,max = max_theta_range)
      initial_mu = rnorm(p,mean = b_exp,sd = se_exp)
    }


    MLE_result =
      cML_estimate(b_exp,b_out,
                   se_exp,se_out,
                   K = K,initial_theta = initial_theta,
                   initial_mu = initial_mu,
                   maxit = maxit)

    Neg_l =
      sum( (b_exp - MLE_result$b_vec)^2 / (2*se_exp^2) ) +
      sum( (b_out - MLE_result$theta * MLE_result$b_vec - MLE_result$r_vec)^2 /
             (2*se_out^2))

    sd_theta = cML_SdTheta(b_exp,b_out,
                           se_exp,se_out,
                           MLE_result$theta,
                           MLE_result$b_vec,
                           MLE_result$r_vec)

    theta_v_RandomCandidate = c(theta_v_RandomCandidate,MLE_result$theta)
    sd_v_RandomCandidate = c(sd_v_RandomCandidate,sd_theta)
    l_v_RandomCandidate = c(l_v_RandomCandidate,Neg_l)
    invalid_RandomCandidate = rbind(invalid_RandomCandidate,
                                    as.numeric(MLE_result$r_vec))

  }
  if(sum(is.na(sd_v_RandomCandidate)))
  {
    warning(paste("May not converge to minimums with some given",
                  "start points and maximum number of iteration,",
                  "lead to Fisher Information matrices not positive definite.",
                  "Could try increasing number of iterations (maxit)",
                  "or try different start points. Note: If",
                  "multiple random start points are used,",
                  "this warning does not likely affect result."))
  }
  min_neg_l = which.min(l_v_RandomCandidate)

  theta_est = theta_v_RandomCandidate[min_neg_l]
  sd_est = sd_v_RandomCandidate[min_neg_l]
  l_est = l_v_RandomCandidate[min_neg_l]
  r_est = invalid_RandomCandidate[min_neg_l,]

  return(list(theta = theta_est,
              se = sd_est,
              l = l_est,
              r_est = r_est
  )
  )
}


mr_cML <- function(b_exp,b_out,
                   se_exp,se_out,
                   K_vec = 0:(length(b_exp) - 2),
                   random_start = 0,
                   maxit = 100,
                   random_seed = 0,
                   n)
{
  if(random_seed)
  {
    set.seed(random_seed)
  }

  rand_theta = NULL
  rand_sd = NULL
  rand_l = NULL
  invalid_mat = NULL
  for(K_value in K_vec)
  {
    rand_res = cML_estimate_randomC(b_exp = b_exp,
                                   b_out = b_out,
                                   se_exp = se_exp,
                                   se_out = se_out,
                                   K = K_value,
                                   random_start = random_start,
                                   maxit = maxit)
    rand_theta = c(rand_theta,rand_res$theta)
    rand_sd = c(rand_sd,rand_res$se)
    rand_l = c(rand_l,rand_res$l)
    invalid_mat = rbind(invalid_mat,rand_res$r_est)
  }

  ### get result
  theta_v = rand_theta
  sd_v = rand_sd
  l_v = rand_l

  # cML-MA-BIC
  BIC_vec = log(n) * K_vec + 2 * l_v
  BIC_vec = BIC_vec - min(BIC_vec)
  weight_vec = exp(-1/2 * BIC_vec)
  weight_vec = weight_vec/sum(weight_vec)
  MA_BIC_theta = sum(theta_v * weight_vec)
  MA_BIC_se = sum(weight_vec * sqrt(sd_v^2 + (theta_v - MA_BIC_theta)^2),
                  na.rm = TRUE)
  MA_BIC_p = pnorm(-abs(MA_BIC_theta/MA_BIC_se))*2

  # cML-BIC
  BIC_vec = log(n) * K_vec + 2 * l_v
  BIC_vec = BIC_vec - min(BIC_vec)
  min_ind = which.min(BIC_vec)
  BIC_theta = theta_v[min_ind]
  BIC_se = sd_v[min_ind]
  BIC_p = pnorm(-abs(BIC_theta/BIC_se))*2
  BIC_invalid = which(invalid_mat[min_ind,]!=0)

  # cML-MA-AIC
  AIC_vec = 2*K_vec + 2*l_v
  AIC_vec = AIC_vec - min(AIC_vec)
  weight_vec = exp(-1/2 * AIC_vec)
  weight_vec = weight_vec/sum(weight_vec)
  MA_AIC_theta = sum(theta_v * weight_vec)
  MA_AIC_se = sum(weight_vec * sqrt(sd_v^2 + (theta_v - MA_AIC_theta)^2),
                  na.rm = TRUE)
  MA_AIC_p = pnorm(-abs(MA_AIC_theta/MA_AIC_se))*2

  # cML-AIC
  AIC_vec = 2*K_vec + 2*l_v
  AIC_vec = AIC_vec - min(AIC_vec)
  min_ind = which.min(AIC_vec)
  AIC_theta = theta_v[min_ind]
  AIC_se = sd_v[min_ind]
  AIC_p = pnorm(-abs(AIC_theta/AIC_se))*2
  AIC_invalid = which(invalid_mat[min_ind,]!=0)


  return(list(MA_BIC_theta = MA_BIC_theta,
              MA_BIC_se = MA_BIC_se,
              MA_BIC_p = MA_BIC_p,
              BIC_theta = BIC_theta,
              BIC_se = BIC_se,
              BIC_p = BIC_p,
              BIC_invalid = BIC_invalid,
              BIC_vec = log(n) * K_vec + 2 * l_v,
              MA_AIC_theta = MA_AIC_theta,
              MA_AIC_se = MA_AIC_se,
              MA_AIC_p = MA_AIC_p,
              AIC_theta = AIC_theta,
              AIC_se = AIC_se,
              AIC_p = AIC_p,
              AIC_invalid = AIC_invalid)
         )
}


cML_estimate_DP <- function(b_exp,b_out,
                            se_exp,se_out,
                            K,num_pert = 200,
                            random_start = 0,
                            maxit = 100)
{
  p = length(b_exp)
  theta_v = NULL
  sd_v = NULL
  l_v = NULL
  for(pt_ind in 1:num_pert)
  {
    b_exp_new = b_exp + rnorm(p,0,1)*se_exp
    b_out_new = b_out + rnorm(p,0,1)*se_out
    MLE_result = cML_estimate_randomC(b_exp = b_exp_new,
                                     b_out = b_out_new,
                                     se_exp = se_exp,
                                     se_out = se_out,
                                     K = K,
                                     random_start = random_start,
                                     maxit = maxit)

    theta_v = c(theta_v,MLE_result$theta)
    sd_v = c(sd_v,MLE_result$se)
    l_v = c(l_v,MLE_result$l)
  }

  return(list(theta_v = theta_v,
              se_v = sd_v,
              l_v = l_v))
}


mr_cML_DP <- function(b_exp,b_out,
                      se_exp,se_out,
                      K_vec = 0:(length(b_exp) - 2),
                      random_start = 0,
                      random_start_pert = 0,
                      maxit = 100,
                      num_pert = 200,
                      random_seed = 0,
                      n)
{
  if(random_seed)
  {
    set.seed(random_seed)
  }

  rand_theta = NULL
  rand_sd = NULL
  rand_l = NULL
  invalid_mat = NULL
  rand_pert_theta = NULL
  rand_pert_sd = NULL
  rand_pert_l = NULL
  for(K_value in K_vec)
  {
    rand_res = cML_estimate_randomC(b_exp = b_exp,
                                   b_out = b_out,
                                   se_exp = se_exp,
                                   se_out = se_out,
                                   K = K_value,
                                   random_start = random_start,
                                   maxit = maxit)
    rand_pert_res = cML_estimate_DP(b_exp = b_exp,
                                    b_out = b_out,
                                    se_exp = se_exp,
                                    se_out = se_out,
                                    K = K_value,
                                    num_pert = num_pert,
                                    random_start = random_start_pert,
                                    maxit = maxit)
    rand_theta = c(rand_theta,rand_res$theta)
    rand_sd = c(rand_sd,rand_res$se)
    rand_l = c(rand_l,rand_res$l)
    invalid_mat = rbind(invalid_mat,rand_res$r_est)
    rand_pert_theta = rbind(rand_pert_theta,rand_pert_res$theta_v)
    rand_pert_sd = rbind(rand_pert_sd,rand_pert_res$se_v)
    rand_pert_l = rbind(rand_pert_l,rand_pert_res$l_v)
  }

  ### get result
  theta_v = rand_theta
  sd_v = rand_sd
  l_v = rand_l
  theta_pt_v = rand_pert_theta
  sd_pt_v = rand_pert_sd
  l_pt_v = rand_pert_l

  var_mat = sd_pt_v^2
  K_vec = 0:(nrow(theta_pt_v) - 1)
  numer_perturb = ncol(theta_pt_v)
  l_pt_v = rowMeans(l_pt_v)
  sd_pt_v = sqrt(diag(var(t(theta_pt_v))))
  theta_pt_mat = theta_pt_v
  theta_pt_v = rowMeans(theta_pt_v)

  # cML-MA-BIC
  BIC_vec = log(n) * K_vec + 2 * l_v
  BIC_vec = BIC_vec - min(BIC_vec)
  weight_vec = exp(-1/2 * BIC_vec)
  weight_vec = weight_vec/sum(weight_vec)
  MA_BIC_theta = sum(theta_v * weight_vec)
  MA_BIC_se = sum(weight_vec * sqrt(sd_v^2 + (theta_v - MA_BIC_theta)^2),
                  na.rm = TRUE)
  MA_BIC_p = pnorm(-abs(MA_BIC_theta/MA_BIC_se))*2

  # cML-BIC
  BIC_vec = log(n) * K_vec + 2 * l_v
  BIC_vec = BIC_vec - min(BIC_vec)
  min_ind = which.min(BIC_vec)
  BIC_theta = theta_v[min_ind]
  BIC_se = sd_v[min_ind]
  BIC_p = pnorm(-abs(BIC_theta/BIC_se))*2
  BIC_invalid = which(invalid_mat[min_ind,]!=0)

  # cML-MA-BIC-DP
  BIC_vec = log(n) * K_vec + 2 * l_pt_v
  BIC_vec = BIC_vec - min(BIC_vec)
  weight_vec = exp(-1/2 * BIC_vec)
  weight_vec = weight_vec/sum(weight_vec)
  MA_BIC_DP_theta = sum(theta_pt_v * weight_vec)
  MA_BIC_DP_se = sum(weight_vec * sqrt(sd_pt_v^2 + (theta_pt_v - MA_BIC_DP_theta)^2),
                     na.rm = TRUE)
  MA_BIC_DP_p = pnorm(-abs(MA_BIC_DP_theta/MA_BIC_DP_se))*2

  # cML-BIC-DP
  BIC_vec = log(n) * K_vec + 2 * l_pt_v
  BIC_vec = BIC_vec - min(BIC_vec)
  min_ind = which.min(BIC_vec)
  BIC_DP_theta = theta_pt_v[min_ind]
  BIC_DP_se = sd_pt_v[min_ind]
  BIC_DP_p = pnorm(-abs(BIC_DP_theta/BIC_DP_se))*2

  # GOF 1
  BIC_vec = log(n) * K_vec + 2 * l_v
  BIC_vec = BIC_vec - min(BIC_vec)
  min_ind = which.min(BIC_vec)
  pt_sd_v = sd_pt_v[min_ind]
  origin_sd_v = sd_v[min_ind]
  more_var = var(var_mat[min_ind,])
  #n = numer_perturb
  x = theta_pt_mat[min_ind,]
  sd_x = sqrt((mean((x - mean(x))^4) -
                 (numer_perturb-3)/(numer_perturb-1)*var(x)^2)/numer_perturb +
                more_var)
  Tstat =
    (origin_sd_v^2 - pt_sd_v^2)/(sd_x)
  GOF1_p = pnorm(-abs(Tstat))*2

  # GOF 2
  BIC_vec = log(n) * K_vec + 2 * l_v
  BIC_vec = BIC_vec - min(BIC_vec)
  min_ind = which.min(BIC_vec)
  pt_sd_v = sd_pt_v[min_ind]
  origin_sd_v = sd_v[min_ind]
  more_var = var(var_mat[min_ind,])

  Tstat =
    (origin_sd_v^2 - pt_sd_v^2)/(sqrt(2/(numer_perturb-1)*pt_sd_v^4 +
                                        more_var))
  GOF2_p = pnorm(-abs(Tstat))*2

  return(list(MA_BIC_theta = MA_BIC_theta,
              MA_BIC_se = MA_BIC_se,
              MA_BIC_p = MA_BIC_p,
              BIC_theta = BIC_theta,
              BIC_se = BIC_se,
              BIC_p = BIC_p,
              BIC_invalid = BIC_invalid,
              MA_BIC_DP_theta = MA_BIC_DP_theta,
              MA_BIC_DP_se = MA_BIC_DP_se,
              MA_BIC_DP_p = MA_BIC_DP_p,
              BIC_DP_theta = BIC_DP_theta,
              BIC_DP_se = BIC_DP_se,
              BIC_DP_p = BIC_DP_p,
              GOF1_p = GOF1_p,
              GOF2_p = GOF2_p
              ))
}
